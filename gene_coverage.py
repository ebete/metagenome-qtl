#!/usr/bin/env python3

import argparse
import logging
from pathlib import Path


def get_orfs(gff_fname):
    logging.info("Reading features from %s ...", gff_fname.name)
    d = dict()
    with gff_fname.open("rt") as f:
        for line in f:
            if not line or line[0] == '#':
                continue
            cols = line.strip().split("\t")
            if cols[2] != "CDS":
                continue

            contig_id = cols[0]
            start = int(cols[3])
            end = int(cols[4])
            orf_id = cols[8].split(";")[0][3:]

            d.setdefault(contig_id, list()).append((start, end, orf_id))
    return d


def get_depths(bed_fname, contig_ids):
    logging.info("Reading alignment depths from %s ...", bed_fname.name)
    d = dict()
    with bed_fname.open("rt") as f:
        for line in f:
            if not line:
                continue
            cols = line.strip().split("\t")

            contig_id = cols[0]
            start = int(cols[1])
            end = int(cols[2])
            depth = int(cols[3])

            if not contig_id in contig_ids:
                continue

            d.setdefault(contig_id, list()).append((start, end, depth))
    return d


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_file", metavar="GFF", help="GFF of reference.")
    parser.add_argument("bed_file", metavar="BED", help="Reference alignment depth BED file.")
    args = parser.parse_args()

    logging.basicConfig(format="[%(asctime)s] %(message)s", level=logging.INFO)

    gff_file = Path(args.gff_file)
    bed_file = Path(args.bed_file)

    orf_regions = get_orfs(gff_file)
    depths = get_depths(bed_file, orf_regions.keys())

    logging.info("Writing average depth per ORF ...")
    print("sample_id", "orf_id", "avg_depth", sep="\t")
    for contig_id in orf_regions:
        depth = depths[contig_id]
        depth_idx = 0
        depth_start, depth_end, depth_val = depth[depth_idx]
        for orf_start, orf_end, orf_id in orf_regions[contig_id]:
            avg_depth = 0
            while orf_start > depth_end and depth_idx < len(depth):
                depth_idx += 1
                depth_start, depth_end, depth_val = depth[depth_idx]
            while orf_end > depth_start and depth_idx < len(depth) - 1:
                overlap = min(orf_end, depth_end) - max(orf_start, depth_start)
                avg_depth += depth_val * overlap
                depth_idx += 1
                depth_start, depth_end, depth_val = depth[depth_idx]
            if avg_depth != 0:
                avg_depth /= (orf_end - orf_start)
            print(bed_file.name, orf_id, f"{avg_depth:.2f}", sep="\t")
