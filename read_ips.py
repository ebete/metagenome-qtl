#!/usr/bin/env python3

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

if __name__ == "__main__":
    ips_xml = Path(sys.argv[1])

    ns = {"ips": "http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5"}
    tree = ET.parse(ips_xml)
    root = tree.getroot()
    for protein in root:
        qname = protein.find("ips:xref", ns).attrib["id"]
        for matches in protein.find("ips:matches", ns):
            program = re.match(r"^\{.+\}([\w\-]+)\-match$", matches.tag, re.IGNORECASE)
            if not program:
                continue
            program = str(program.group(1))
            print(program)
