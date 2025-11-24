#!/usr/bin/env python3
import sys
from Bio import SeqIO
from Bio import AlignIO
from collections import Counter
import matplotlib.pyplot as plt

usage = 'Usage: ' + sys.argv[0] + ' <stk file> <output pdf file>'
if len(sys.argv) != 3:
    print(usage)
    exit()

input_stk_file = sys.argv[1]
output_file = sys.argv[2]
ref_id = 'Y09347.1'

align = AlignIO.read(input_stk_file,'stockholm')
dot_bracket = align.column_annotations["secondary_structure"]
print(dot_bracket)

seq_data = []
for record in align:
    if record.id == ref_id:
        seq_data.insert(0,(record.id,record.seq))        
    else:
        seq_data.append((record.id,record.seq))

id_list,seq_list = zip(*seq_data)
seq_list = list(seq_list)
n_seqs = len(seq_list)
    
positions = []
mut_rates = []
gap_rates = []
#set ref_seq to our Y09352.1 seq (enter where that falls in the list)
ref_seq   = seq_list[0]
#loop through each column in filem then create percent identity for each position column
pos_counter = 1
for index, column in enumerate(zip(*seq_list), start=1):
    ref_base    = ref_seq[index-1]
    counts      = Counter(column)
    match_count = counts.get(ref_base, 0)
    gap_count   = counts.get('-',0)
    percent_id  = match_count / n_seqs * 100
    mut_rate    = 100.0 - percent_id
    # parsimonious assumption
    gap_rate    = min(100*(gap_count / n_seqs), 100*(1 - gap_count / n_seqs))
    var_counts  = {b: c for b, c in counts.items() if b != ref_base}
    positions.append(pos_counter)
    mut_rates.append(mut_rate)
    gap_rates.append(gap_rate)
    if ref_base != '-':
        pos_counter += 1
    # print(f"Pos {pos}: {percent_id:.2f}% identity "
    #      f"({match_count}/{n_seqs}), variants: {var_counts}")

plot_positions = []
structure_flags = []

pos_counter = 1
for i, (char, base) in enumerate(zip(dot_bracket, ref_seq)):
    if base != "-":
        plot_positions.append(pos_counter)
        structure_flags.append(char)
        pos_counter += 1
#highlight bases that have ( or ) w/ blue                                                     
for pos, sflag in zip(plot_positions, structure_flags):
    if sflag in '()':
        plt.axvspan(pos-0.5, pos + 0.5, color='lightblue', alpha=0.3)
    
plt.plot(positions, mut_rates, label='mutation rate')
plt.plot(positions, gap_rates, label='insertion/deletion rate')
plt.xlabel('Position (nt)')
plt.ylabel('Fraction')
plt.legend()
plt.savefig('output_file')
