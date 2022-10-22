import pytest
from variants_lib.format_variants import (
    detect_variation_type,
    VariationType,
    factor_sequence,
    describe_alt,
    decode_indel,
)


DELETION_EXAMPLES_MONOALLELIC = [("TCTAC", "T"), ("TA", "T"), ("CTG", "C")]

DELETION_EXAMPLES_MULTIALLELIC = [
    [("TA", "T"), ("TAA", "T")],
    [("GAC", "G"), ("GACAC", "G")],
]

INSERTION_EXAMPLES_MONOALLELIC = [("T", "TCTAC"), ("T", "TA"), ("C", "CTG")]

INSERTION_EXAMPLES_MULTIALLELIC = [
    [("T", "TA"), ("T", "TAA")],
    [("G", "GAC"), ("G", "GACAC")],
]

INDEL_EXAMPLES = [
    [("T", "TTTTGTTTG"), ("TTTTGTTTG", "T"), ("T", "TTTTG"), ("TTTTG", "T")],
    [("G", "GAA"), ("GAA", "G"), ("G", "GAAA"), ("GA", "G")],
]


def test_detect_snp():
    monoallelic_snp = [("A", "G")]
    assert detect_variation_type(monoallelic_snp) == VariationType.SNP
    multiallelic_snp = [("A", "G"), ("A", "T")]
    assert detect_variation_type(multiallelic_snp) == VariationType.SNP


def test_detect_mnp():
    monoallelic_mnp = [("CAT", "TGC")]
    assert detect_variation_type(monoallelic_mnp) == VariationType.MNP
    multiallelic_mnp = [("CAT", "TGC"), ("CAT", "CGT")]
    assert detect_variation_type(multiallelic_mnp) == VariationType.MNP


@pytest.mark.parametrize("ref,alt", DELETION_EXAMPLES_MONOALLELIC)
def test_detect_deletion_monoallelic(ref, alt):
    assert detect_variation_type([(ref, alt)]) == VariationType.DELETION


@pytest.mark.parametrize("alleles", DELETION_EXAMPLES_MULTIALLELIC)
def test_detect_deletion_multiallelic(alleles):
    assert detect_variation_type(alleles) == VariationType.DELETION


@pytest.mark.parametrize("ref,alt", INSERTION_EXAMPLES_MONOALLELIC)
def test_detect_insertion_monoallelic(ref, alt):
    assert detect_variation_type([(ref, alt)]) == VariationType.INSERTION


@pytest.mark.parametrize("alleles", INSERTION_EXAMPLES_MULTIALLELIC)
def test_detect_insertion_multiallelic(alleles):
    assert detect_variation_type(alleles) == VariationType.INSERTION


@pytest.mark.parametrize("alleles", INDEL_EXAMPLES)
def test_detect_indels(alleles):
    assert detect_variation_type(alleles) == VariationType.INDEL


SEQUENCE_EXAMPLES = [
    ("A", ("A", 1)),
    ("AA", ("A", 2)),
    ("AAA", ("A", 3)),
    ("AB", ("AB", 1)),
    ("ABB", ("ABB", 1)),
    ("ABAB", ("AB", 2)),
    ("ABCBC", ("ABCBC", 1)),
]


@pytest.mark.parametrize("sequence,expected", SEQUENCE_EXAMPLES)
def test_factor_sequence(sequence, expected):
    infix, count = factor_sequence(sequence)
    assert expected == (infix, count)
    assert infix * count == sequence


INSDEL_EXAMPLES = [
    ("A", "AA", "INS(A)_1"),
    ("A", "ABC", "INS(BC)_1"),
    ("A", "ABABA", "INS(BA)_2"),
    ("AB", "A", "DEL(B)_1"),
    ("ABB", "AB", "DEL(B)_1"),
    ("ABCABABAB", "ABC", "DEL(AB)_3"),
]


@pytest.mark.parametrize("ref,alt,expected", INSDEL_EXAMPLES)
def test_describe_alt(ref, alt, expected):
    assert describe_alt(ref, alt) == expected


INDEL_DECODE_EXAMPLES = [
    (
        [
            ("A", "AGT", 0),
            ("AGT", "A", 0),
            ("AGTGT", "A", 0),
            ("A", "AGTGTGT", 0),
            ("A", "AGTGTGTGT", 0),
            ("A", "AGTGT", 0),
        ],
        ("A", "", "A"),
    ),
    (
        [
            ("A", "AGT", 1),
            ("AGT", "A", 0),
            ("AGTGT", "A", 0),
            ("A", "AGTGTGT", 0),
            ("A", "AGTGTGTGT", 0),
            ("A", "AGTGT", 0),
        ],
        ("A", "AGT", "I"),
    ),
    (
        [
            ("A", "AGT", 0),
            ("AGT", "A", 1),
            ("AGTGT", "A", 0),
            ("A", "AGTGTGT", 0),
            ("A", "AGTGTGTGT", 0),
            ("A", "AGTGT", 0),
        ],
        ("AGT", "A", "I"),
    ),
    (
        [
            ("A", "AGT", 0),
            ("AGT", "A", 0),
            ("AGTGT", "A", 1),
            ("A", "AGTGTGT", 0),
            ("A", "AGTGTGTGT", 0),
            ("A", "AGTGT", 0),
        ],
        ("AGTGT", "A", "I"),
    ),
    (
        [
            ("A", "AGT", 0),
            ("AGT", "A", 0),
            ("AGTGT", "A", 0),
            ("A", "AGTGTGT", 1),
            ("A", "AGTGTGTGT", 0),
            ("A", "AGTGT", 0),
        ],
        ("A", "AGTGTGT", "I"),
    ),
]


@pytest.mark.parametrize("rows,expected", INDEL_DECODE_EXAMPLES)
def test_decode_indel(rows, expected):
    assert decode_indel(rows) == expected


DELETION_DECODE_EXAMPLES = [
    (
        [["TACATATATACACATGTATATACAC", "T", 1], ["TAC", "T", 0]],
        ("TACATATATACACATGTATATACAC", "T", "I"),
    ),
    ([["TACATATATACACATGTATATACAC", "T", 0], ["TAC", "T", 1]], ("TAC", "T", "I")),
    (
        [
            ["GTGT", "G", 1],
            ["GTGTGTGTGT", "G", 0],
            ["GT", "G", 0],
            ["GTGTGT", "G", 0],
        ],
        ("GTGT", "G", "I"),
    ),
    ([["GT", "G", 0], ["GTT", "G", 1]], ("GTT", "G", "I")),
]


@pytest.mark.parametrize("rows,expected", DELETION_DECODE_EXAMPLES)
def test_decode_deletion(rows, expected):
    assert decode_indel(rows) == expected


INSERTION_DECODE_EXAMPLES = [
    (
        [
            ["T", "TATATATATATATATATA", 1],
            ["T", "TATATATATATATATA", 0],
            ["T", "TATATATATATATA", 0],
            ["T", "TATATATATATATATATATA", 0],
            ["T", "TATATATATATATATATATATA", 0],
        ],
        ("T", "TATATATATATATATATA", "I"),
    ),
    (
        [["C", "CTGTG", 1], ["C", "CTGTGTGTGTG", 0], ["C", "CTG", 0]],
        ("C", "CTGTG", "I"),
    ),
]


@pytest.mark.parametrize("rows,expected", INSERTION_DECODE_EXAMPLES)
def test_decode_insertion(rows, expected):
    assert decode_indel(rows) == expected
