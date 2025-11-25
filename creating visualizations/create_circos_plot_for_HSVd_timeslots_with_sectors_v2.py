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
                'Cusa_sRNAseq_T1_HSVd_R1_HSVd_forward.depth',
                'Cusa_sRNAseq_T1_HSVd_R2_HSVd_forward.depth',
                'Cusa_sRNAseq_T1_HSVd_R3_HSVd_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T1_HSVd_R1_HSVd_reverse.depth',
                'Cusa_sRNAseq_T1_HSVd_R2_HSVd_reverse.depth',
                'Cusa_sRNAseq_T1_HSVd_R3_HSVd_reverse.depth'
            ],
        },
        "T2": {
            "forward": [
                'Cusa_sRNAseq_T2_HSVd_R1_HSVd_forward.depth',
                'Cusa_sRNAseq_T2_HSVd_R2_HSVd_forward.depth',
                'Cusa_sRNAseq_T2_HSVd_R3_HSVd_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T2_HSVd_R1_HSVd_reverse.depth',
                'Cusa_sRNAseq_T2_HSVd_R2_HSVd_reverse.depth',
                'Cusa_sRNAseq_T2_HSVd_R3_HSVd_reverse.depth'
            ],
        },
        "T3": {
            "forward": [
                'Cusa_sRNAseq_T3_HSVd_R1_HSVd_forward.depth',
                'Cusa_sRNAseq_T3_HSVd_R2_HSVd_forward.depth',
                'Cusa_sRNAseq_T3_HSVd_R3_HSVd_forward.depth'
            ],
            "reverse": [
                'Cusa_sRNAseq_T3_HSVd_R1_HSVd_reverse.depth',
                'Cusa_sRNAseq_T3_HSVd_R2_HSVd_reverse.depth',
                'Cusa_sRNAseq_T3_HSVd_R3_HSVd_reverse.depth'
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
    
    
    r_label = 112
    for (label, (start, end)), color in zip(region_labels.items(), region_colors):
        mid = (start + end) / 2
        for pos in (start, end):
            # add the text at the correct genome position (mid) and radius (r_label)
            hsvd_sector.text(
                label,
                x=mid,       # genomic coordinate (within sector range)
                r=r_label,   # distance from center
                size=10,
                color=color,
            )
            """
            hsvd_sector.text(
                f"({start}-{end})",
                x=mid,
                r=105,         # place just outside the outermost track
                size=9,
                color=color,
                ha="center"
            )
            
            hsvd_sector.line(
                rad_lim=[pos,pos],
                r=[95, 100],  # draws a full radial line
                color=color,
                lw=0.5,
                zorder=10,
            )

            hsvd_sector.line(
                rad_lim=[pos,pos],
                r=[95, 100],
                color=color,
                lw=0.5,
                zorder=10,
            )
            """

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
    # Plot coverage tracks
    # -----------------------------
    min_r = 50
    max_r = 100
    gap = 0
    n_tracks = 3
    height = (max_r - min_r) / n_tracks
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

        if idx == 2:
            track.xticks_by_interval(interval=10)
        
        if idx == 0:
            circos.rect(r_lim=(50, 51), deg_lim=(0, 360), fc="black", zorder=10)

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
        if idx == 2:
            for start,end,color in rect_list:
                track.rect(start, end, r_lim = (95,100), facecolor=color, edgecolor="none")
                #track.rect((50,100), None, fc=color, ec="None")

    # -----------------------------
    # Plot base-pairing links
    # -----------------------------
    for bp in bps:
        circos.link(("HSVd", bp[0], bp[0]), ("HSVd", bp[1], bp[1]))
    
    fig = circos.plotfig(figsize=(7, 7))
    circos.savefig("HSVd_timeslots_with_bps_and_sectors.pdf")

    """
    # --- Debug / force-full-circle block (put at end of main(), before exiting) ---

    # Force circos-level settings (best-effort for older/newer pycirclize)
    try:
        circos.start = 0
        circos.end = 360
        circos.space = 0
    except Exception:
        pass

    # Force sector to cover full 360 and radial domain to 0-100
    for s in circos.sectors:
        try:
            s.start = 0
            s.end = 360
        except Exception:
            pass
        # degree_size attr may be read-only in some versions, ignore if not settable
        if hasattr(s, "degree_size"):
            try:
                s.degree_size = 360.0
            except Exception:
                pass
        # ensure internal radial bounds (if present)
        if hasattr(s, "r_min"):
            try:
                s.r_min = 0
                s.r_max = 100
            except Exception:
                pass

    # Make sure all tracks / rects are inside 0..100
    # (If you used variables min_r/max_r earlier, consider reducing max_r slightly)
    # e.g., max_r = min(max_r, 98)
    # Create the figure (this is the important step)
    fig = circos.plotfig(figsize=(8, 8))

    # Force equal aspect and wide limits so nothing is clipped
    ax = fig.axes[0]
    ax.set_aspect("equal", adjustable="box")

    # set generous x/y limits (these numbers should contain the circle)
    ax.set_xlim(-120, 120)
    ax.set_ylim(-120, 120)

    # Debug prints
    print("MATPLOTLIB BACKEND:", matplotlib.get_backend(), file=sys.stderr)
    try:
        import pycirclize
        print("pycirclize version:", getattr(pycirclize, "__version__", "unknown"), file=sys.stderr)
    except Exception:
        print("pycirclize import failed", file=sys.stderr)

    for s in circos.sectors:
        print(f"Sector='{s.name}' start={getattr(s,'start',None)} end={getattr(s,'end',None)} size={getattr(s,'degree_size',None)}", file=sys.stderr)

    print("Axes xlim, ylim:", ax.get_xlim(), ax.get_ylim(), file=sys.stderr)

    # Save both PNG and PDF (PNG is easier to preview)
    png_out = "HSVd_full_circle_debug.png"
    pdf_out = "HSVd_full_circle_debug.pdf"
    fig.savefig(png_out, bbox_inches="tight", dpi=300)
    fig.savefig(pdf_out, bbox_inches="tight")
    print("Saved:", png_out, pdf_out, file=sys.stderr)

    # check genome length
    print("genome_length:", genome_length, file=sys.stderr)

    # depth files min/max
    global_min = None
    global_max = None
    for time_slot, strands in depth_files.items():
        for strand, files in strands.items():
            for file in files:
                try:
                    df_tmp = pd.read_csv(file, sep="\t", header=None, names=["chrom", "pos", "depth"])
                except Exception as e:
                    print(f"Failed reading {file}: {e}", file=sys.stderr)
                    continue
                if df_tmp.empty:
                    print(f"{file} is empty!", file=sys.stderr)
                    continue
                pmin = int(df_tmp["pos"].min())
                pmax = int(df_tmp["pos"].max())
                print(f"{file}: pos_min={pmin}, pos_max={pmax}", file=sys.stderr)
                if global_min is None or pmin < global_min:
                    global_min = pmin
                if global_max is None or pmax > global_max:
                    global_max = pmax

    print("global depth pos range:", global_min, global_max, file=sys.stderr)

    # mfe base-pairs min/max
    if len(bps) > 0:
        bp_positions = [p for pair in bps for p in pair]
        print("bps min, max:", np.min(bp_positions), np.max(bp_positions), file=sys.stderr)
    else:
        print("No base-pairs parsed from MFE.", file=sys.stderr)
    """
if __name__ == "__main__":
    main()



