from pytest import fixture
import pytest
from variants_lib import athena, Variant


@fixture
def canonical_rsids_mock(monkeypatch):
    def mock(rsids):
        ddb = {"rs123": "rs456", "rs321": "rs654"}
        result = {rsid: new_rsid for rsid, new_rsid in ddb.items() if rsid in rsids}
        for rsid in rsids:
            if rsid not in result:
                result[rsid] = rsid
        return result

    monkeypatch.setattr(athena, "canonical_rsids", mock)
    yield


@fixture
def get_variants_mock(monkeypatch):
    def mock(rsid_list):
        variants = {
            "rs456": [Variant(chrom="chr1", pos=1001, ref="A", alt="G")],
            "rs999": [Variant(chrom="chr1", pos=1002, ref="C", alt="T")],
        }
        return {rsid: variants[rsid] for rsid in rsid_list if rsid in variants}

    monkeypatch.setattr(athena, "get_variants", mock)
    yield


def test_format_genotype():
    variants_snp = [("A", "G", 0, 1)]
    assert athena.format_genotype(variants_snp) == (("A", "G", "A"), ("A", "G", "G"))
    variants_multiallelic_snp = [("A", "G", 0, 0), ("A", "T", 0, 1)]
    assert athena.format_genotype(variants_multiallelic_snp) == (
        ("A", "", "A"),
        ("A", "T", "T"),
    )
    variants_insertion = [("A", "AG", 0, 0), ("A", "AGG", 1, 0)]
    assert athena.format_genotype(variants_insertion) == (
        ("A", "AGG", "I"),
        ("A", "", "A"),
    )
    variants_deletion = [("AGG", "A", 0, 1), ("AG", "A", 0, 0)]
    assert athena.format_genotype(variants_deletion) == (
        ("A", "", "A"),
        ("AGG", "A", "I"),
    )
    variants_indel = [
        ("AGTGT", "A", 0, 0),
        ("AGT", "A", 1, 0),
        ("A", "AGTGTGT", 0, 1),
    ]
    assert athena.format_genotype(variants_indel) == (
        ("AGT", "A", "I"),
        ("A", "AGTGTGT", "I"),
    )


@pytest.mark.usefixtures("canonical_rsids_mock")
class TestEndowRsidWithVariant:
    def test_no_rsid(self):
        sites = [
            Variant(chrom="chr2", pos=2002, ref="A", alt="G"),
        ]
        assert athena.endow_rsid_with_variant(sites) == [
            Variant(chrom="chr2", pos=2002, ref="A", alt="G", rsid=None),
        ]

    def test_rsid_only(self, get_variants_mock):
        sites = ["rs123", "rs999"]
        assert athena.endow_rsid_with_variant(sites) == [
            Variant(chrom="chr1", pos=1001, ref="A", alt="G", rsid="rs123"),
            Variant(chrom="chr1", pos=1002, ref="C", alt="T", rsid="rs999"),
        ]

    def test_rsid_and_variant(self, get_variants_mock):
        sites = ["rs123", Variant(chrom="chr2", pos=2002, ref="A", alt="G")]
        assert athena.endow_rsid_with_variant(sites) == [
            Variant(chrom="chr1", pos=1001, ref="A", alt="G", rsid="rs123"),
            Variant(chrom="chr2", pos=2002, ref="A", alt="G", rsid=None),
        ]

    def test_missing_rsid(self, get_variants_mock):
        sites = ["rs321"]
        assert athena.endow_rsid_with_variant(sites) == []


class TestExtractGt:
    def test_biallelic_variant(self):
        assert athena.extract_gt(
            [
                {
                    "chrom": "chr1",
                    "pos": 1001,
                    "ref": "A",
                    "alt": "G",
                    "rsid": "rs456",
                    "gt1": 0,
                    "gt2": 1,
                }
            ],
            [Variant(chrom="chr1", pos=1001, ref="A", alt="G", rsid="rs123")],
        ) == {Variant(chrom="chr1", pos=1001, ref="A", alt="G", rsid="rs123"): (0, 1)}

    def test_multiallelic_variant(self):
        variants_with_gt_dict = [
            {
                "chrom": "chr1",
                "pos": 1002,
                "ref": "C",
                "alt": "T",
                "rsid": "rs999",
                "gt1": 1,
                "gt2": 0,
            },
            {
                "chrom": "chr1",
                "pos": 1002,
                "ref": "C",
                "alt": "G",
                "rsid": "rs999",
                "gt1": 0,
                "gt2": 1,
            },
        ]
        rsid_with_variants = [
            Variant(chrom="chr1", pos=1002, ref="C", alt="T", rsid="rs999"),
            Variant(chrom="chr1", pos=1002, ref="C", alt="G", rsid="rs999"),
        ]
        assert athena.extract_gt(variants_with_gt_dict, rsid_with_variants) == {
            Variant(chrom="chr1", pos=1002, ref="C", alt="T", rsid="rs999"): (1, 0),
            Variant(chrom="chr1", pos=1002, ref="C", alt="G", rsid="rs999"): (0, 1),
        }
