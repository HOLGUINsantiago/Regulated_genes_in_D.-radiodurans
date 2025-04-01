import os
import pandas as pd
import re

def parse_fasta(fasta_file):
    peak_to_chr = {}
    with open(fasta_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                match = re.match(r">(\S+) (\S+) from:", line)
                if match:
                    peak_id, chromosome = match.groups()
                    peak_to_chr[peak_id] = chromosome
    return peak_to_chr

def update_csv_files(directory, peak_to_chr):
    for filename in os.listdir(directory):
        if filename.endswith(".tsv"):
            csv_path = os.path.join(directory, filename)
            df = pd.read_csv(csv_path, sep="\t")  
            
            df["sequence_name"] = df["sequence_name"].map(peak_to_chr).fillna(df["sequence_name"])
            
            df.to_csv(csv_path, sep="\t", index=False)

fasta_file = "Fimo136/FinalPeaks_DdrO_Paper.fasta"
directory = "Fimo136"  

peak_to_chr = parse_fasta(fasta_file)
update_csv_files(directory, peak_to_chr)
