# extract_and_blast.py
# Extracts regions from bins based on hit_list.tsv, writes to multi-FASTA.
# Runs blastn on local core_nt DB.
# Filters blast output to max bitscore per query, unique sscinames.

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from subprocess import run
import sys

# Read hit_list.tsv (skip if needed, but headers start with #FILE)
hits = pd.read_csv('hit_list.tsv', sep='\t')

# Prepare multi-FASTA
extracted = []
cgs = []
for _, row in hits.iterrows():
    bin_file = row['#FILE']  # e.g., SemiBin_149.fa
    sequence = row['SEQUENCE']  # e.g., contig_7264
    start = row['START'] - 1  # 0-based for slicing
    end = row['END']
    strand = row['STRAND']
    
    # Parse bin FASTA, find contig
#    with open(bin_file, 'r') as f:
    with open(f"{sys.argv[1]}/{bin_file}", 'r') as f:                # for reprocess.sh
        records = list(SeqIO.parse(f, 'fasta'))
        contig = next(r for r in records if r.id == sequence)
    
    # Extract region
    seq = contig.seq[start:end]
    if strand == '-':
        seq = seq.reverse_complement()  # Reverse complement for - strand

    # Extract region
    seq2 = contig.seq
    if strand == '-':
        seq2 = seq.reverse_complement()  # Reverse complement for - strand
    
    # Add to FASTA with meaningful ID
    id = f"{row['#FILE']}_{sequence}_{row['START']}-{row['END']}"
    extracted.append(SeqRecord(seq, id=id, description=''))

    cgs_id = f"{row['#FILE']}_{sequence}"
    cgs.append(SeqRecord(seq2, id=cgs_id, description=''))

# Write multi-FASTA
with open('extracted_regions.fasta', 'w') as f:
    SeqIO.write(extracted, f, 'fasta')

with open('extracted_contigs.fasta', 'w') as f:
    SeqIO.write(cgs, f, 'fasta')

# Run blastn (assume $BLASTDB set)
run([
    'blastn', '-query', 'extracted_regions.fasta', '-db', 'core_nt', #'"core_nt ref_prok_rep_genomes"',
    '-outfmt', '6 qseqid sseqid staxids sscinames pident qcovs length mismatch bitscore',
    '-out', 'blast_out.tsv', '-num_threads', '32'
])

# Process blast output
blast = pd.read_csv('blast_out.tsv', sep='\t', header=None, names=[
    'qseqid', 'sseqid', 'staxids', 'sscinames', 'pident', 'qcovs', 'length', 'mismatch', 'bitscore'
])

# Group by qseqid, find max bitscore
max_bits = blast.groupby('qseqid')['bitscore'].max().reset_index()

# Filter rows with max bitscore
filtered = pd.merge(blast, max_bits, on=['qseqid', 'bitscore'])

filtered['sscinames'] = filtered['sscinames'].str.split().str[:2].str.join(' ')

# Unique by sscinames per qseqid (drop duplicates)
unique = filtered.drop_duplicates(subset=['qseqid', 'sscinames'])

# Select columns: qseqid, sscinames, pident, qcovs
report = unique[['qseqid', 'sscinames', 'pident', 'qcovs', 'length']]

# Save to TSV
report.to_csv('blast_report.tsv', sep='\t', index=False)
