import sys
from Bio import SeqIO
from pycirclize import Circos
import numpy as np
import pandas as pd
import math
import matplotlib
import matplotlib.pyplot as plt

usage = 'Usage: ' + sys.argv[0] + ' <fasta file> ' + ' <mfe file> '

if len(sys.argv) < 3:
    print(usage)
    exit()

def main():
    fasta_file = sys.argv[1]
    genome_length = get_genome_size(fasta_file)

    mfe_file = sys.argv[2]
    bps = extract_bps_from_mfe(mfe_file)

    # -----------------------------
    # CONFIGURATION
    # -----------------------------
    depth_files = {
        "T1": {
            "forward": [
                'Cusa_sRNAseq_T1_HSVd_R1_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T1_HSVd_R2_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T1_HSVd_R3_HSVd_shift_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T1_HSVd_R1_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T1_HSVd_R2_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T1_HSVd_R3_HSVd_shift_reverse.depth'
            ],
        },
        "T2": {
            "forward": [
                'Cusa_sRNAseq_T2_HSVd_R1_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T2_HSVd_R2_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T2_HSVd_R3_HSVd_shift_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T2_HSVd_R1_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T2_HSVd_R2_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T2_HSVd_R3_HSVd_shift_reverse.depth'
            ],
        },
        "T3": {
            "forward": [
                'Cusa_sRNAseq_T3_HSVd_R1_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T3_HSVd_R2_HSVd_shift_forward.depth',
                'Cusa_sRNAseq_T3_HSVd_R3_HSVd_shift_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T3_HSVd_R1_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T3_HSVd_R2_HSVd_shift_reverse.depth',
                'Cusa_sRNAseq_T3_HSVd_R3_HSVd_shift_reverse.depth'
            ],
        }
    }
    
    # -----------------------------
    # BUILD CIRCOS
    # -----------------------------
    # Define one continuous genome sector
    sectors = {"HSVd": (0, genome_length)}
    circos = Circos(sectors, space=0)
    #print(circos.sectors)
    #for s in circos.sectors:
        #print(s)
    colors = {"forward": "steelblue", "reverse": "tomato"}

    # -----------------------------
    # Annotate labeled subregions
    # -----------------------------
    region_labels = {
        "TL1": (1, 31),
        "P1": (32, 62),
        "C1": (63, 100),
        "V1": (101, 124),
        "TR": (125, 175),
        "V2": (176, 198),
        "C2": (198, 238),
        "P2": (239, 269),
        "TL2": (270, 297),
    }

    rect_list = [
            (1,31, "#66c2a5"), (32, 62, "#fc8d62"), (63,100, "#8da0cb"), (101,124,"#e78ac3"), (125,175,"#a6d854"), (176,198,"#e78ac3"), (198,238,"#8da0cb"),(239,269,"#fc8d62"),(270,297            ,"#66c2a5")
    ]

    region_colors = [
        "#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854",
        "#e78ac3", "#8da0cb", "#fc8d62", "#66c2a5"
    ]

    hsvd_sector = circos.sectors[0]
    print(f"Sector start={hsvd_sector.start}, end={hsvd_sector.end}, size={hsvd_sector.deg_size}")

    # -----------------------------
    # Compute global max depth
    # -----------------------------
    max_depth = 0
    for strands in depth_files.values():
        for files in strands.values():
            for file in files:
                df_tmp = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                file_max = df_tmp["depth"].max()
                if file_max > max_depth:
                    max_depth = file_max
    track_vmin, track_vmax = -max_depth, max_depth

    # -----------------------------
    # Plot base-pairing links
    # -----------------------------
    min_r = 50
    max_r = 120
    gap = 0
    n_tracks = 4
    height = (max_r - min_r) / n_tracks

    for bp in bps:
        circos.link(("HSVd", bp[0], bp[0]), ("HSVd", bp[1], bp[1]), min_r + height/2 - 5, min_r + height/2 - 5)

    # -----------------------------
    # Plot coverage tracks
    # -----------------------------
    times = ["10 dpi", "17 dpi", "24 dpi"]
    for idx, (time_slot, strands) in enumerate(depth_files.items()):
        inner_r = min_r + idx * (height + gap)
        outer_r = inner_r + height

        track = hsvd_sector.add_track((inner_r, outer_r))
     # Force full circle and prevent clipping
        """
        circos.start = 0
        circos.end = 360
        circos.space = 0
        """

        if idx == 0:
            circos.rect(r_lim=(min_r + height/2 - 5, min_r + height/2 - 4), deg_lim=(0, 360), fc="black", zorder=10)

        # Add a dedicated link track very close to the ideogram
        #link_track = hsvd_sector.add_track((0, 50))

        circos.text(times[idx], r=(inner_r + outer_r) / 2 + 5, deg=90,
                    orientation="horizontal", size=8, color="black")

        # Plot forward and reverse replicates
        for strand, files in strands.items():
            for file in files:
                df = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                if strand == "reverse":
                    df["depth"] = df["depth"] * -1
                df_sub = df[(df["pos"] >= 0) & (df["pos"] <= genome_length)]

                track.line(
                    df_sub["pos"].values,
                    df_sub["depth"].values,
                    color=colors[strand],
                    vmin=track_vmin,
                    vmax=track_vmax,
                )
        # Plot sector rectangles
        inner_r = min_r + (idx + 1) * (height + gap)
        outer_r = inner_r + height
        track = hsvd_sector.add_track((inner_r, outer_r))
        for start,end,color in rect_list:
            track.rect(start, end, r_lim = (max_r - 5, max_r), facecolor=color, edgecolor="none")
            # Add x ticks
            track.xticks_by_interval(interval=10)
            #track.rect((50,100), None, fc=color, ec="None")

        # Plot labels
        for (label, (start, end)), color in zip(region_labels.items(), region_colors):
        mid = (start + end) / 2
        for pos in (start, end):
            # add the text at the correct genome position (mid) and radius (r_label)
            hsvd_sector.text(
                label,
                x=mid,       # genomic coordinate (within sector range)
                r=outer_r + 7,   # distance from center
                size=10,
                color=color,
            )


        fig = circos.plotfig(figsize=(7, 7))
        circos.savefig("HSVd_timeslots_with_bps_and_sectors.pdf")

def get_genome_size(fasta_file):
    data = SeqIO.parse(fasta_file, "fasta")
    for record in data:
        print(len(record.seq))
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
        elif dot_bracket_str[idx] == ')':
            bp_first_half = stack.pop()
            base_pairs.append((bp_first_half, idx))
    return base_pairs

