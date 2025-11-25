import sys
from Bio import SeqIO
from pycirclize import Circos
import numpy as np
import pandas as pd
import math

usage = 'Usage: ' + sys.argv[0] + ' <fasta file> ' + ' <mfe file> '

if len(sys.argv) < 2:
    print(usage)
    exit()

def get_genome_size(fasta_file):
    data = SeqIO.parse(fasta_file, "fasta")
    for record in data:
        #print(record.id)
        return len(record.seq)

def extract_bps_from_mfe(mfe_file):
    with open(mfe_file, 'r') as mfe:
        dot_bracket_str = ''
        for line in mfe:
            if line.startswith('.') or line.startswith('('):
                dot_bracket_str += line.split()[0]
    stack = []
    base_pairs = []
    for idx in range(len(dot_bracket_str)):
        if dot_bracket_str[idx] == '(':
            stack.append(idx)
        if dot_bracket_str[idx] == ')':
            bp_first_half = stack.pop()
            base_pairs.append((bp_first_half, idx))
    return base_pairs

def main():
    
    fasta_file = sys.argv[1]
    genome_length = get_genome_size(fasta_file)

    mfe_file = sys.argv[2]
    bps = extract_bps_from_mfe(mfe_file)

    # -----------------------------
    # CONFIGURATION
    # -----------------------------
    # Map files into groups: time_slot -> {"forward": [...], "reverse": [...]}
    depth_files = {
            "T1": {
                "forward": ['Cusa_sRNAseq_T1_HSVd_R1_HSVd_forward.depth', 'Cusa_sRNAseq_T1_HSVd_R2_HSVd_forward.depth', 'Cusa_sRNAseq_T1_HSVd_R3_HSVd_forward.depth'],
                "reverse": ['Cusa_sRNAseq_T1_HSVd_R1_HSVd_reverse.depth', 'Cusa_sRNAseq_T1_HSVd_R2_HSVd_reverse.depth', 'Cusa_sRNAseq_T1_HSVd_R3_HSVd_reverse.depth']
                },
            "T2": {
                "forward": ['Cusa_sRNAseq_T2_HSVd_R1_HSVd_forward.depth', 'Cusa_sRNAseq_T2_HSVd_R2_HSVd_forward.depth', 'Cusa_sRNAseq_T2_HSVd_R3_HSVd_forward.depth'],
                "reverse": ['Cusa_sRNAseq_T2_HSVd_R1_HSVd_reverse.depth', 'Cusa_sRNAseq_T2_HSVd_R2_HSVd_reverse.depth', 'Cusa_sRNAseq_T2_HSVd_R3_HSVd_reverse.depth']
                },
            "T3": {
                "forward": ['Cusa_sRNAseq_T3_HSVd_R1_HSVd_forward.depth', 'Cusa_sRNAseq_T3_HSVd_R2_HSVd_forward.depth', 'Cusa_sRNAseq_T3_HSVd_R3_HSVd_forward.depth'],
                "reverse": ['Cusa_sRNAseq_T3_HSVd_R1_HSVd_reverse.depth', 'Cusa_sRNAseq_T3_HSVd_R1_HSVd_reverse.depth', 'Cusa_sRNAseq_T3_HSVd_R1_HSVd_reverse.depth']
                }
            }
    #print(depth_files)
    

    # -----------------------------
    # BUILD CIRCOS
    # -----------------------------
    circos = Circos({"HSVd": genome_length}, space=5)
    
    colors = {"forward":"steelblue", "reverse":"tomato"}
    
    sector = circos.sectors[0]

    # --- Compute max depth across forward & reverse for this time_slot ---
    max_depth = 0
    for idx, (time_slot, strands) in enumerate(depth_files.items()):
        for strand, files in strands.items():   
            for file in files:
                df_tmp = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                file_max = df_tmp["depth"].max()
                if file_max > max_depth:
                    max_depth = file_max

    track_vmin, track_vmax = -max_depth, max_depth

    # --- Add tracks for each time slot ---
    min_r = 50
    max_r = 100
    gap = 0
    n_tracks = 3
    height = (max_r - min_r)/ n_tracks
    times = ["10 dpi","17 dpi", "24 dpi"] 
    for idx, (time_slot, strands) in enumerate(depth_files.items()):
        # Set the inner and outer radius for the track
        inner_r = min_r + idx * (height + gap)
        outer_r = inner_r + height
        outer_r = int(outer_r)
        
        track = sector.add_track((inner_r, outer_r))
        
        if idx == 2:
            track.xticks_by_interval(interval=20)

        if idx == 0:
            circos.rect(r_lim=(49, 50), deg_lim=(0, 360), fc="black", zorder=10)

        circos.text(times[idx], r = (inner_r + outer_r) / 2, deg = 90, orientation = "horizontal", size = 12, color = "black")
        
        # --- Plot forward & reverse replicates ---
        for strand, files in strands.items():
            for file in files:
                df = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                if strand == "reverse":
                    df["depth"] = df["depth"] * -1

                df_sub = df[(df["pos"] >= 1) & (df["pos"] <= genome_length)]

                track.line(
                    df_sub["pos"].values,
                    df_sub["depth"].values,
                    color=colors[strand],
                    vmin=track_vmin,
                    vmax=track_vmax,
                )
        # --- Plot links ---
        for bp in bps:
            bp_first_half = bp[0]
            bp_second_half = bp[1]
            circos.link(("HSVd", bp_first_half, bp_first_half), ("HSVd", bp_second_half, bp_second_half))
        
    circos.savefig("HSVd_timeslots_with_bps_7.pdf")
main()


