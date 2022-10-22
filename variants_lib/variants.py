from collections import defaultdict
from typing import Iterable, List, Dict, Tuple
import os
from variants_lib import Variant, Locus
from variants_lib.utils import get_ddb_client

### Get build38 coordinates and variants from rsid

TABLE_NAME = os.environ.get("BUILD38_TABLE", "genome_reference_build38")


def get_loci(rsids: List[str]) -> Dict[str, Locus]:
    """Return a dictionary mapping rsid to Locus for those rsids that
    are found in the build38 table. If no rsid is found, return and empty dict."""
    client = get_ddb_client("build38")
    resp = client.batch_get_item(
        RequestItems={TABLE_NAME: {"Keys": [{"rsid": {"S": rsid}} for rsid in rsids]}}
    )
    try:
        items = resp["Responses"][TABLE_NAME]
    except KeyError:
        return {}
    return {
        item["rsid"]["S"]: Locus(item["chrom"]["S"], int(item["pos"]["N"]))
        for item in items
    }


def get_variants(rsids: List[str], normalized=True) -> Dict[str, List[Variant]]:
    client = get_ddb_client("build38")
    resp = client.batch_get_item(
        RequestItems={TABLE_NAME: {"Keys": [{"rsid": {"S": rsid}} for rsid in rsids]}}
    )
    try:
        items = resp["Responses"][TABLE_NAME]
    except KeyError:
        return {}
    result = defaultdict(list)
    for item in items:
        rsid, chrom, pos, ref, alts = (
            item["rsid"]["S"],
            item["chrom"]["S"],
            int(item["pos"]["N"]),
            item["ref"]["S"],
            [alt["S"] for alt in item["alt"]["L"]],
        )
        if not normalized:
            for alt in alts:
                result[rsid].append(Variant(chrom=chrom, pos=pos, ref=ref, alt=alt))
        else:
            for normalized_ref, normalized_alt in normalize_variant(ref, alts):
                result[rsid].append(
                    Variant(
                        chrom=chrom, pos=pos, ref=normalized_ref, alt=normalized_alt
                    )
                )
    return result


# TODO: Implement (chrom, pos) -> List[Variant]


def normalize_variant(ref: str, alts: List[str]) -> Iterable[Tuple[str, str]]:
    "Yield parsimonious (ref, alt) pairs. The input is expected to be left aligned."
    # if not all(ref.startswith(alt) or alt.startswith(ref) for alt in alts):
    #     raise ValueError()
    for alt in alts:
        k = 0
        while ref[-(1 + k) :] == alt[-(1 + k) :]:
            k += 1
        yield (ref[: len(ref) - k], alt[: len(alt) - k])
