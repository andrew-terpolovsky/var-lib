from pytest import fixture
from unittest.mock import Mock
from variants_lib.variants import get_loci, get_variants, normalize_variant
from variants_lib import variants as variants_module
from variants_lib import Locus, Variant
from .utils import MockDdbClient


def test_variants_equality():
    v1 = Variant("chr1", 1, "A", "T", "rs123")
    v2 = Variant("chr1", 1, "A", "T", "rs456")
    assert v1 != v2


def test_variant_id():
    v1 = Variant("chr1", 1000, "A", "T", "rs123")
    assert v1.id == "1_1000_A_T"


@fixture(autouse=True)
def mock_build38_table(monkeypatch):
    monkeypatch.setattr(
        variants_module,
        "get_ddb_client",
        lambda role: MockDdbClient(
            table_name="genome_reference_build38",
            items=[
                {
                    "rsid": {"S": "rs123"},
                    "chrom": {"S": "chr1"},
                    "pos": {"N": "10000"},
                    "ref": {"S": "A"},
                    "alt": {"L": [{"S": "G"}]},
                },
                {
                    "rsid": {"S": "rs456"},
                    "chrom": {"S": "chr2"},
                    "pos": {"N": "15000"},
                    "ref": {"S": "AGTGT"},
                    "alt": {"L": [{"S": "A"}, {"S": "AGT"}, {"S": "AGTGTGT"}]},
                },
            ],
            key_name="rsid",
        )
        if role == "build38"
        else None,
    )


def test_get_locus():
    assert get_loci(["rs123"]) == {"rs123": Locus(chrom="chr1", pos=10000)}
    assert get_loci(["rs780"]) == {}


class TestGetVariants:
    def test_get_variants_none(self):
        assert get_variants(["rs000"]) == {}

    def test_get_variants_snp(self):
        assert get_variants(["rs123"]) == {
            "rs123": [Variant(chrom="chr1", pos=10000, ref="A", alt="G")]
        }

    def test_get_variants_indel_normalized(self):
        assert get_variants(["rs456"], normalized=True) == {
            "rs456": [
                Variant(chrom="chr2", pos=15000, ref="AGTGT", alt="A"),
                Variant(chrom="chr2", pos=15000, ref="AGT", alt="A"),
                Variant(chrom="chr2", pos=15000, ref="A", alt="AGT"),
            ]
        }

    def test_get_variants_indel_not_normalized(self):
        assert get_variants(["rs456"], normalized=False) == {
            "rs456": [
                Variant(chrom="chr2", pos=15000, ref="AGTGT", alt="A"),
                Variant(chrom="chr2", pos=15000, ref="AGTGT", alt="AGT"),
                Variant(chrom="chr2", pos=15000, ref="AGTGT", alt="AGTGTGT"),
            ]
        }


def test_normalize_variant():
    assert list(normalize_variant("A", ["AG", "AGG"])) == [("A", "AG"), ("A", "AGG")]
    assert list(normalize_variant("AGT", ["A", "AGTGT"])) == [
        ("AGT", "A"),
        ("A", "AGT"),
    ]
