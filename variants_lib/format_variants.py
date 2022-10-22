from typing import List, Tuple, Optional
from variants_lib import VariationType

### Format variants for user-friendly display


def detect_variation_type(refalt_alleles: List[Tuple[str, str]]) -> VariationType:
    ref_alleles, alt_alleles = zip(*refalt_alleles)
    ref_lengths = set(len(ref) for ref in ref_alleles)
    alt_lengths = set(len(alt) for alt in alt_alleles)
    lengths = ref_lengths | alt_lengths
    if len(lengths) == 1:
        if lengths == {1}:
            return VariationType.SNP
        else:
            return VariationType.MNP

    # For the two below conditions we use the fact that variants are normalized
    if all(len(alt) == 1 and ref.startswith(alt) for ref, alt in refalt_alleles):
        return VariationType.DELETION

    if all(len(ref) == 1 and alt.startswith(ref) for ref, alt in refalt_alleles):
        return VariationType.INSERTION

    return VariationType.INDEL


def factor_sequence(sequence: str) -> Tuple[str, int]:
    for length in range(1, len(sequence) // 2 + 1):
        i = len(sequence) // length
        if sequence == sequence[:length] * i:
            return (sequence[:length], i)
    return (sequence, 1)


def describe_alt(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return alt
    if ref.startswith(alt):
        verb = "DEL"
        sequence = ref[len(alt) :]
    else:
        verb = "INS"
        sequence = alt[len(ref) :]
    infix, count = factor_sequence(sequence)
    return f"{verb}({infix})_{count}"


def decode_indel(
    rows: List[Tuple[str, str, Optional[int]]], factor=False
) -> Tuple[str, str, str]:
    if all(call is None for ref, alt, call in rows):  # Will most likely never happen
        return ("", "", "")
    variation_type = detect_variation_type([(ref, alt) for ref, alt, _ in rows])
    for row in rows:
        ref, alt, call = row
        if call == 1:
            if variation_type == VariationType.SNP:
                return ref, alt, alt
            return ref, alt, "I"  # stands for "INDEL"
    ref = rows[0][0][0]
    alt = ""
    return ref, alt, ref
