"""
Utility script to map allele sequences to variant types.
The script expectes to find a `data` directory containing gzipped JSON files.
The JSON schema for the input data is
 - `chrom`: string
 - `pos`: int
 - `ref`: string
 - `alt`: string
The script loads all the loci in memory before unstacking the multiallelic sites,
 because of that it deals with the chromosomes separately (in order to avoid OOM issues).
The output is a single json.gz file containing all the input data with the alleles unstacked
(for multilallelic sites) and a `type` field containing the variant type.
"""
import json
from collections import defaultdict
import gzip
import os

from variants_lib.format_variants import detect_variation_type

BASE_DIR = os.path.dirname(__file__)
DATA_DIR = os.path.join(BASE_DIR, "data")

CHROM_LIST = [f"chr{i}" for i in range(1, 23)] + ["chrX"]


def load_chrom_data(
    file_list: list[str], chrom: str
) -> dict[int, list[tuple[str, str]]]:
    """Load the data for a given chromosome and return the result as a dictionary
    (each locus is a key, and each value is a list of ref/alt alleles)."""
    d: dict[int, list[tuple[str, str]]] = defaultdict(list)
    for file in file_list:
        print(f"Loading data from {file}")
        for row in map(json.loads, gzip.open(os.path.join(BASE_DIR, DATA_DIR, file))):
            if row["chrom"] != chrom:
                continue
            d[row["pos"]].append((row["ref"], row["alt"]))
    return d


def main():
    files = os.listdir(DATA_DIR)
    print(f"Found {len(files)} files.")
    for chrom in CHROM_LIST:
        print("Loading data for chrom {0}".format(chrom))
        data = load_chrom_data(files, chrom)
        print("All data loaded.")
        with gzip.open("output.json.gz", "at", newline="\n") as output:
            for locus in data:
                type = detect_variation_type(data[locus]).name
                output.write(
                    json.dumps(
                        {
                            "chrom": chrom,
                            "pos": locus,
                            "type": type,
                            "alleles": data[locus],
                        }
                    )
                    + "\n"
                )
    print("All done")


if __name__ == "__main__":
    main()
