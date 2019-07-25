#!/usr/bin/env python3

import sys
from pathlib import Path
import yaml


def main(args):
    data_loc = Path(args[0])
    structure = {
        "data": dict()
    }

    for f in data_loc.glob("*.fq.gz"):
        if not f.is_file():
            continue
        parts = f.name.split(".")[0].split("_")
        sample_id = parts[0]
        direction = parts[-1]
        lane_nr = parts[-2]

        structure["data"]\
            .setdefault(sample_id, dict())\
            .setdefault(direction, list())\
            .append(f.absolute().as_posix())
    return structure


if __name__ == "__main__":
    d = main(sys.argv[1:])
    sys.stdout.write(yaml.dump(d, default_flow_style=False))
