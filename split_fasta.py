#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path

from Bio import SeqIO


def split_fasta(fasta_file, output_path, split_count):
    logging.info("Reading FASTA records from %s ...", fasta_file)
    seqs = dict()
    seqcount = 0
    with fasta_file.open("r") as fasta:
        for record in SeqIO.parse(fasta, "fasta"):
            seqcount += 1
            fasta_bin = seqcount % split_count + 1
            seqs.setdefault(fasta_bin, list()).append(record)

    for fasta_bin, fasta_records in seqs.items():
        new_fasta = output_path.joinpath(f"{fasta_bin}_of_{split_count}_{fasta_file.name}")
        logging.info("Writing bin %d to %s [%d records] ...", fasta_bin, new_fasta, len(fasta_records))
        with new_fasta.open("w") as fasta_out:
            SeqIO.write(fasta_records, fasta_out, "fasta")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="[%(asctime)s] %(message)s")
    # get command-line arguments
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="FASTA file to split.", metavar="FASTA")
    parser.add_argument("-n", "--splits", help="Number of files to spread the FASTA records over.", metavar="N",
                        type=int, default=2)
    parser.add_argument("-o", "--output", help="Output directory of the split FASTA files.", metavar="DIR",
                        default="./")
    args = parser.parse_args()

    fasta_path = Path(args.input_fasta).expanduser().resolve().absolute()
    output_path = Path(args.output).expanduser().resolve().absolute()
    split_count = args.splits

    split_fasta(fasta_path, output_path, split_count)
