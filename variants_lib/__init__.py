from enum import Enum
from dataclasses import dataclass
from typing import Optional

VARIANT_ID_DELIMITER = "_"


class VariationType(Enum):
    INDEL = "INDEL"
    INSERTION = "INSERTION"
    DELETION = "DELETION"
    SNP = "SNP"
    MNP = "MNP"  # Multi-Nucleotide Polymorphism (we haven't any yet.)


@dataclass(frozen=True)
class Locus:
    chrom: str
    pos: int
    rsid: Optional[str] = None


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    rsid: Optional[str] = None
    genotype: Optional[str] = None

    def __eq__(self, other) -> bool:
        if not isinstance(other, Variant):
            return False
        if self.rsid != other.rsid:
            return False
        return self.id == other.id

    def __hash__(self) -> int:
        return hash((self.chrom, self.pos, self.ref, self.alt))

    @property
    def id(self) -> str:
        if not self.alt:
            return ""
        chrom = self.chrom.lstrip("chr")  # remove leading 'chr' if it exists
        return VARIANT_ID_DELIMITER.join([chrom, str(self.pos), self.ref, self.alt])


# Use cases:
# I have rsid, I want genotype for that rsid. EASY.
# I have a locus, I want genotype for that locus. HARD.
# 1. There could be more than one rsid for that locus. Example: rs1294270260 and rs150461910 at 8_1672228. The first
#    an SNP, the second an INDEL.
# 2. In that case we could give the user both results. The rsids would just come from the data.
# I have a variant (chrom pos ref alt), what do I want? Two distinct cases:
# 1. If the variant is a SNP, or is monoallelic we just want the genotype.
# 2. If the variant is more complex, e.g. an INDEL, we want to compute the genotype for this INDEL variant.
# 3. Or maybe we just want the corresponding gt without any context?

# The opentargets use case is different than the recommendation use case.


# If I query using an rsid, I want the same rsid back. If I query using a locus, I want the rsid
# that is in the data (if there is any).
# If I want to know the canonical rsid, it is my responsibility to fetch it using the appropriate
# function.
