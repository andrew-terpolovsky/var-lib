import dataclasses
import logging
from uuid import UUID
from typing import (
    Optional,
    List,
    Tuple,
    Dict,
    Union,
    TypedDict,
    Callable,
)
import os

import boto3
from variants_lib import Variant, Locus
from variants_lib.merges import canonical_rsids
from variants_lib.variants import get_variants
from variants_lib.format_variants import decode_indel
import pyarrow.dataset as ds
import pyarrow
from pyarrow import fs

TABLE_COLUMNS = ["chrom", "pos", "rsid", "ref", "alt", "gt1", "gt2"]
VariantWithGTDict = TypedDict(
    "VariantWithGTDict",
    {
        "chrom": str,
        "pos": int,
        "ref": str,
        "alt": str,
        "gt1": Optional[int],
        "gt2": Optional[int],
        "rsid": Optional[str],
    },
)


def _get_genome_file_etl_metadata_table():
    genome_file_ddb = os.environ["USER_GENOME_FILE_ETL_DDB"]
    table = boto3.resource("dynamodb").Table(genome_file_ddb)
    return table


def get_parquet_path(file_id: Union[UUID, str]) -> Optional[str]:
    user_genome_file_etl_bucket = os.environ["USER_GENOME_FILE_ETL_BUCKET"]
    genome_file_ddb = _get_genome_file_etl_metadata_table()
    resp = genome_file_ddb.get_item(Key={"file_id": str(file_id)})
    item = resp["Item"]
    if item.get("fetchable", False):
        # return pyarrow backend parquet path
        parquet_path = f'{user_genome_file_etl_bucket}/user_genome_files/parquets/file_id={item["file_id"]}'  # pylint: disable=line-too-long
    else:
        parquet_path = None
    return parquet_path


#
#  Get data from parquet using pyarrow
#


def get_filter(variants: List[Variant]) -> pyarrow.compute.Expression:
    filter_expr = None
    for variant in variants:
        # we filter only using chrom and pos. this is more efficient than
        # filtering on chrom, pos, ref, alt, by about 20%.
        # filtering w.r.t ref/alt happens later
        filter_ = (ds.field("chrom") == variant.chrom) & (
            ds.field("pos") == variant.pos
        )
        if filter_expr is None:
            filter_expr = filter_
        else:
            filter_expr = filter_expr | filter_
    return filter_expr


def read_dataset(base_path: str) -> ds.Dataset:
    s3_fs = fs.S3FileSystem(region="us-east-1")
    return ds.dataset(
        base_path,
        format="parquet",
        partitioning=ds.partitioning(flavor="hive"),
        filesystem=s3_fs,
        partition_base_dir=base_path,
    )


def _validate_sites(sites: List[Union[str, Variant]]) -> None:
    invalid_rsids = [
        site for site in sites if isinstance(site, str) and not site.startswith("rs")
    ]
    if invalid_rsids:
        raise ValueError(f'Invalid rsids: {", ".join(invalid_rsids)}')
    invalid_sites = [
        str(site) for site in sites if not isinstance(site, (str, Variant))
    ]
    if invalid_sites:
        raise ValueError(f'Invalid sites: {", ".join(invalid_sites)}')


def endow_rsid_with_variant(sites: List[Union[str, Variant]]) -> List[Variant]:
    """Return a list computed by replacing each rsid with (rsid, variant) and each
    locus with (None, locus)."""
    rsids = [site for site in sites if isinstance(site, str)]
    result: List[Variant] = []
    if rsids:
        provided_rsid_to_canonical_rsid: Dict[str, str] = {}  # map
        # the given rsids to the corresponding canonical rsids.

        # We deal with slices of 100 rsids because that's the limit for
        # ddb:BatchGetItem.
        for i in range(0, len(rsids), 100):
            rsids_slice = rsids[i : i + 100]
            provided_rsid_to_canonical_rsid.update(**canonical_rsids(rsids_slice))

        canonical_rsid_list = list(provided_rsid_to_canonical_rsid.values())

        variants: Dict[str, List[Variant]] = {}  # map canonical rsids to variants
        # Likewise we slice again for the same reason.
        for i in range(0, len(canonical_rsid_list), 100):
            canonical_rsid_slice = canonical_rsid_list[i : i + 100]
            variants.update(**get_variants(canonical_rsid_slice))

    else:  # Not used, we define them anyway to appease the static checkers.
        provided_rsid_to_canonical_rsid = {}
        variants = {}
    for site in sites:
        if isinstance(site, Variant):
            result.append(site)
        elif isinstance(site, str):
            canonical_rsid = provided_rsid_to_canonical_rsid[site]
            if canonical_rsid in variants:
                result += [
                    dataclasses.replace(variant, rsid=site)
                    for variant in variants[canonical_rsid]
                ]
        else:
            raise ValueError(f"Invalid site: {site}")
    return result


def get_raw_gt(
    variants: List[Variant],
    dataset: ds.Dataset,
) -> Dict[Variant, Tuple[Optional[int], Optional[int]]]:
    """Get the raw genotypes for the given variants. Maps the variant to (gt1, gt2).
    For multiallelic variants, further processing is usually desirable.
    """
    variants_with_gt_dict: List[VariantWithGTDict] = dataset.to_table(
        columns=TABLE_COLUMNS, filter=get_filter(variants)
    ).to_pylist()
    # filter w.r.t ref/alt. more efficient to do it here than in the pyarrow dataset filtering
    def filter_over_ref_alt(variant: VariantWithGTDict) -> bool:
        return any(
            variant_.chrom == variant["chrom"]
            and variant_.pos == variant["pos"]
            and variant_.ref == variant["ref"]
            and variant_.alt == variant["alt"]
            for variant_ in variants
        )

    variants_with_gt_dict = list(filter(filter_over_ref_alt, variants_with_gt_dict))

    return extract_gt(variants_with_gt_dict, variants)


def variant_with_client_rsid(variant: Variant, client_rsid: Optional[str]) -> Variant:
    if client_rsid is None:
        return variant
    return dataclasses.replace(variant, rsid=client_rsid)


def forget_rsid(variant: Variant) -> Variant:
    return dataclasses.replace(variant, rsid=None)


def extract_gt(
    variants_with_gt_dict: List[VariantWithGTDict],
    rsid_with_variants: List[Variant],
) -> Dict[Variant, Tuple[Optional[int], Optional[int]]]:
    """Extract the genotypes from the variants_with_gt_dict."""

    variant_to_rsid = {
        forget_rsid(variant): variant.rsid for variant in rsid_with_variants
    }

    # We map the list of dict to a mapping {variant -> (gt1, gt2)}
    variants_with_gt: Dict[Variant, Tuple[Optional[int], Optional[int]]] = {
        Variant(
            chrom=variant_dict["chrom"],
            pos=variant_dict["pos"],
            ref=variant_dict["ref"],
            alt=variant_dict["alt"],
            rsid=variant_dict.get("rsid"),
        ): (
            variant_dict["gt1"],
            variant_dict["gt2"],
        )
        for variant_dict in variants_with_gt_dict
    }

    # We go through the variants for which an rsid was provided by the client, and replace the
    # rsid from the parquet with the client-provided rsid.

    variants_with_gt = {
        variant_with_client_rsid(variant, variant_to_rsid.get(forget_rsid(variant))): gt
        for variant, gt in variants_with_gt.items()
    }

    return variants_with_gt


def gt_getter(ref: str, alt: str) -> Callable[[Optional[int]], str]:
    def getter(gt: Optional[int]) -> str:
        if gt is None:
            return "."
        return [ref, alt][gt]

    return getter


def format_genotype(
    variants: List[Tuple[str, str, Optional[int], Optional[int]]]
) -> Tuple[Tuple[str, str, str], Tuple[str, str, str]]:
    """Takes a list of (ref, alt, gt1, gt2) and return a pair of properly formatted
    (ref, alt, gt). This is non-trivial only for multiallelic indels."""
    if len(variants) == 1:
        (variant,) = variants
        if is_snp(variant):
            ref, alt, *gt = variant
            gt1, gt2 = map(gt_getter(ref, alt), gt)
            return ((ref, alt, gt1), (ref, alt, gt2))

    variants_gt1 = [(variant[0], variant[1], variant[2]) for variant in variants]
    variants_gt2 = [(variant[0], variant[1], variant[3]) for variant in variants]
    ref_alt_gt1, ref_alt_gt2 = tuple(map(decode_indel, (variants_gt1, variants_gt2)))
    return ref_alt_gt1, ref_alt_gt2


def is_snp(variant: Tuple[str, str, Optional[int], Optional[int]]) -> bool:
    return len(variant[0]) == 1 and len(variant[1]) == 1


def format_raw_gt(
    raw_gt: Dict[Variant, Tuple[Optional[int], Optional[int]]]
) -> Dict[Locus, Tuple[Variant, Variant]]:
    # We deal with the variants separately. Now a variant is a (chrom, pos, rsid). There are theoretically
    # some situations where a given locus has two different variants which cannot be discriminated
    # based on the rsid, because there could be no rsid. We do not deal with these situations yet.
    # One way to do so would be to create a function which, given a number of variations (ref, alt)
    # at a given locus, would split these variations into subset corresponding to different logical variants.
    # This is likely overkill for the time being.
    logical_variants = {
        (variant.chrom, variant.pos, variant.rsid) for variant in raw_gt
    }
    result: Dict[Locus, Tuple[Variant, Variant]] = {}

    logging.info("Found %d logical variants.", len(logical_variants))
    for logical_variant in logical_variants:
        # keep the variants matching the given logical variant
        variants = [
            (variant.ref, variant.alt, gt1, gt2)
            for variant, (gt1, gt2) in raw_gt.items()
            if (variant.chrom, variant.pos, variant.rsid) == logical_variant
        ]

        (ref1, alt1, gt1), (ref2, alt2, gt2) = format_genotype(variants)
        v1 = Variant(
            chrom=logical_variant[0],
            pos=logical_variant[1],
            ref=ref1,
            alt=alt1,
            rsid=logical_variant[2],
            genotype=gt1,
        )
        v2 = Variant(
            chrom=logical_variant[0],
            pos=logical_variant[1],
            ref=ref2,
            alt=alt2,
            rsid=logical_variant[2],
            genotype=gt2,
        )
        result[Locus(*logical_variant)] = (v1, v2)
    return result


def get_genotypes_raw(
    sites: List[Union[str, Variant]],
    file_id: Union[UUID, str],
) -> Optional[Dict[Variant, Tuple[Optional[int], Optional[int]]]]:
    if not sites:
        logging.warning("No sites specified.")
        return {}  # type: ignore
    s3_path = get_parquet_path(file_id)

    if s3_path is None:
        logging.info("Genome file has not been ingested for file_id %s", file_id)
        return None
    _validate_sites(sites)

    # Get the variants (= chrom, pos, ref, alt) for the given rsids.
    rsid_with_variant = endow_rsid_with_variant(
        sites
    )  # list of (rsid, variant) when the rsid
    # was provided by the client, otherwise (None, variant | locus).
    if not rsid_with_variant:
        logging.info(
            "No site (%s) could be matched to a variant.",
            ", ".join(str(site) for site in sites),
        )
        return {}  # type: ignore

    dataset = read_dataset(s3_path)  # lazy read
    return get_raw_gt(rsid_with_variant, dataset)


def get_genotypes(
    sites: List[Union[str, Variant]], file_id: Union[UUID, str]
) -> Optional[Dict[Locus, Tuple[Variant, Variant]]]:
    raw_genotypes = get_genotypes_raw(sites, file_id)
    if not raw_genotypes:
        return raw_genotypes  # type: ignore
    return format_raw_gt(raw_genotypes)
