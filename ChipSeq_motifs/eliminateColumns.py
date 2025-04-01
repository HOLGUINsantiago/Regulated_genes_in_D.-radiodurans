import pandas as pd
from intervaltree import Interval, IntervalTree
import re
import os

def load_gff3(gff3_file):
    gff3_df = pd.read_csv(gff3_file, sep='\t', comment='#', header=None)
    gff3_df.columns = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes']
    genes_df = gff3_df[gff3_df['type'] == 'gene'].sort_values(by=['seqid', 'start'])

    interval_trees = {}
    promoter_regions = {}
    promoter_minus500_regions = {}

    prev_start, prev_end, prev_strand, prev_attributes = None, None, None, None
    
    for seqid, group in genes_df.groupby('seqid'):
        interval_trees[seqid] = IntervalTree()
        promoter_regions[seqid] = IntervalTree()
        promoter_minus500_regions[seqid] = IntervalTree()

        for _, row in group.iterrows():
            start, end, strand, attributes = row['start'], row['end'], row['strand'], row['attributes']
            interval_trees[seqid].add(Interval(start, end, (attributes, strand)))

            if prev_start is not None:
                if prev_strand == '+':
                    promoter_start, promoter_end = prev_end, start
                else:
                    promoter_start, promoter_end = end, prev_start

                if promoter_start < promoter_end:
                    promoter_regions[seqid].add(Interval(promoter_start, promoter_end, prev_attributes))
            
            tss = start if strand == '+' else end  
            promoter_start = max(0, tss - 500) if strand == '+' else tss
            promoter_end = tss if strand == '+' else tss + 500
            promoter_minus500_regions[seqid].add(Interval(promoter_start, promoter_end, attributes))

            prev_start, prev_end, prev_strand, prev_attributes = start, end, strand, attributes
    
    return interval_trees, promoter_regions, promoter_minus500_regions

def find_associated_genes(peak_file, interval_trees, promoter_regions, promoter_minus500_regions):
    peaks_df = pd.read_csv(peak_file, sep='\t', header=0)
    
    associated_genes = []
    for _, peak in peaks_df.iterrows():
        seqid = peak['sequence_name']
        if seqid in interval_trees:
            peak_start, peak_end = peak['start'], peak["start"] + 100 if  (peak['stop'] - peak['start']) > 100 else peak['stop']

            gene_matches = interval_trees[seqid][peak_start:peak_end]
            promoter_matches = promoter_regions[seqid][peak_start:peak_end]
            promoter_minus500_matches = promoter_minus500_regions[seqid][peak_start:peak_end]

            for interval in gene_matches:
                associated_genes.append((seqid, peak_start, peak_end, interval.data[0], "Gene body"))

            for interval in promoter_matches:
                associated_genes.append((seqid, peak_start, peak_end, interval.data, "Promoter region"))

            for interval in promoter_minus500_matches:
                associated_genes.append((seqid, peak_start, peak_end, interval.data, "Proximal promoter (-500)"))
    
    return associated_genes

def extract_gene_name(attributes):
    match = re.search(r'Name=([^;]+)', attributes)
    return match.group(1) if match else "Unknown"

def write_results_to_tsv(results, output_file):
    unique_results = list(set(results))  
    with open(output_file, 'w') as f:
        f.write("Gene_Name\tChromosome\tStart\tEnd\tRegion\n") 
        for gene in unique_results:
            gene_name = extract_gene_name(gene[3])
            f.write(f"{gene_name}\t{gene[0]}\t{gene[1]}\t{gene[2]}\t{gene[4]}\n")

def process_all_files_in_folder(folder, gff3_file, output_file):
    interval_trees, promoter_regions, promoter_minus500_regions = load_gff3(gff3_file)
    all_associated_genes = set()  

    for filename in os.listdir(folder):
        if filename.endswith(".tsv"):  
            peak_file = os.path.join(folder, filename)
            associated_genes = find_associated_genes(peak_file, interval_trees, promoter_regions, promoter_minus500_regions)
            all_associated_genes.update(associated_genes)  

    write_results_to_tsv(all_associated_genes, output_file)
    print(f"Results saved in {output_file}")

gff3_file = "ChipSeq_motifs/seq.gff3"
fimo_folder = "ChipSeq_motifs/Fimo136"
output_file = "ChipSeq_motifs/associated_genes136.tsv"

process_all_files_in_folder(fimo_folder, gff3_file, output_file)

fimo_folder = "ChipSeq_motifs/Fimo452"
output_file = "ChipSeq_motifs/associated_genes452.tsv"

process_all_files_in_folder(fimo_folder, gff3_file, output_file)
