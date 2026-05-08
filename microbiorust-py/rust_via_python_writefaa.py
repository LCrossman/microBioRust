#!/usr/bin/python

import microbiorust
result = microbiorust.parse_gbk("Rhiz3841.gbk.gb")
result.write_faa("rhiz.faa")
