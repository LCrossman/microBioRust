#!/usr/bin/python

import microbiorust as mb
collection = mb.parse_gbk("Rhiz3841.gbk.gb")

with open("rhiz_through_parse.faa", 'w') as f:
    for record in collection.values():
        for info in record.faa().values():
            f.write(f">{info.locus_tag}\n{info.faa}\n")
