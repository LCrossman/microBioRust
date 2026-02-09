import sys
import os
from Bio import SeqIO

def main():
    if len(sys.argv) < 2:
        sys.exit(1)
    
    filepath = sys.argv[1]
    count = 0
    
    # Matching the exact logic of your Rust tool
    for record in SeqIO.parse(filepath, "genbank"):
        genome_seq = record.seq
        for feature in record.features:
            if feature.type == "CDS":
                parts = getattr(feature.location, 'parts', [feature.location])
                for part in parts:
                    _ = str(part.extract(genome_seq).translate(table=11)).split('*')[0]
                    count += 1
    # Minimal output to match Rust tool behavior
    output = open("rhiz_bp.txt", 'w')
    output.write("{}".format(count))

if __name__ == "__main__":
    main()
