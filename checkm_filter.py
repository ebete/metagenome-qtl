#!/usr/bin/env python3

import sys
import csv
from pathlib import Path
import re


_tax_ranks = {
    "k": "kingdom",
    "p": "phylum",
    "c": "class",
    "o": "order",
    "f": "family",
    "g": "genus",
    "s": "species"
}


def main(checkm_dir, metabat_dir, out_dir):
    for path in checkm_dir.glob("*/"):
        if not path.is_dir():
            continue
        csv_file = path.joinpath("checkm_lineagewf.tsv")
        if not csv_file.is_file():
            continue
        taxa_bins = get_good_bins(csv_file)
        print(path, taxa_bins, sep="\t")

        for taxon, bins in taxa_bins.items():
            for bin_name in bins:
                metabat_bin = metabat_dir.joinpath(path.stem, bin_name + ".fa").resolve()
                if not metabat_bin.is_file():
                    continue
                bindir_out = out_dir.joinpath("_".join(taxon))
                binfile_out = bindir_out.joinpath(bin_name + ".fa")
                bindir_out.mkdir(exist_ok=True)
                binfile_out.symlink_to(metabat_bin)
                print(binfile_out, "->", metabat_bin)


def get_good_bins(qa_file):
    taxa_bins = dict()
    with qa_file.open("r", newline="") as f:
        handle = csv.reader(f, delimiter="\t")
        next(handle)
        for record in handle:
            bin_id = record[0]
            taxon = record[1].split(" ")[0].split("__")
            taxa_rank = _tax_ranks.get(taxon[0], "unknown")
            taxa_name = taxon[-1]
            completeness = float(record[11])
            contamination = float(record[12])
            heterogeneity = float(record[13])

            if completeness < 90:
                continue
            if contamination > 5:
                continue

            taxa_bins.setdefault((taxa_rank, taxa_name), list()).append(bin_id)
    return taxa_bins


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))
