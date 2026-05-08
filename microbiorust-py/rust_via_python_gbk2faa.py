#!/usr/bin/python

import microbiorust
result = microbiorust.gbk_to_faa("Rhiz3841.gbk.gb")
outfile = open("rhiz.faa", 'w')
for r in result:
    outfile.write(">{}\n{}\n".format(r, result[r].faa))
