#!/usr/bin/python

import os
import microbiorust as mb

bench_dir = os.path.dirname(os.path.realpath(__file__))
filepath = os.path.join(bench_dir, "Rhiz3841.gbk.gb") 

result = mb.gbk_to_faa_count(filepath)
outfile = open("rhiz.txt", 'w')
outfile.write("{}\n".format(result))
