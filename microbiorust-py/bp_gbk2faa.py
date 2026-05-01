#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq


genbank_file = open("Rhiz3841.gbk.gb",'r')

record = SeqIO.read(genbank_file, "genbank")

proteins = []

for feature in record.features:
    if feature.type == "CDS":
        # The location object handles:
        # - complement()
        # - join()
        # - frames
        # - partials
        nuc_seq = feature.extract(record.seq)

        # Determine translation table (default = 1)
        transl_table = feature.qualifiers.get("transl_table", ["1"])[0]
        transl_table = int(transl_table)

        # If GenBank already provides the translation, use it
        provided_translation = feature.qualifiers.get("translation", [None])[0]

        # Otherwise translate ourselves
        translated = nuc_seq.translate(table=transl_table, to_stop=True)

        # Write a fasta header name
        protein_id = feature.qualifiers.get("protein_id", ["unknown"])[0]
        locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
        
        #print(f">{protein_id}  {locus_tag}")
        #print(translated)
        #print()

        proteins.append((protein_id, locus_tag, str(translated)))

# If you want to save to a FASTA file:
with open("predicted_proteins.faa", "w") as out:
    for pid, tag, aa in proteins:
        out.write(f">{pid} {tag}\n{aa}\n")

print("Saved translations to predicted_proteins.faa")

