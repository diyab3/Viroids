import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

usage = 'Usage: ' + sys.argv[0] + ' < clustering output files > '

if len(sys.argv) < 1:
    print(usage)
    exit()

cluster_files = sys.argv[1:]

fig,ax = plt.subplots()
plt.xticks([10, 20, 30], ["T1", "T2", "T3"])
plt.xlim(0, 40)

for cluster_file in cluster_files:
    df = pd.read_csv(cluster_file, sep="\t", header=None, names=['pos', 'strand', 'num reads', "5' het", 'seq', 'RPM'])
    for idx, abundance in enumerate(df['RPM']):
        if df['strand'].iloc[idx] == '+':
            color = 'steelblue'
        else:
            color = 'tomato'
        if "T1" in cluster_file:
            ax.scatter(10, abundance, c=color)
        elif "T2" in cluster_file:
            ax.scatter(20, abundance, c=color)
        else:
            ax.scatter(30, abundance, c=color)


plt.xlabel('Time')
plt.ylabel('Abundance (in rpm)')

"""
labels_added = {}

for cluster_file in cluster_files:
    df = pd.read_csv(cluster_file, sep="\t", header=None, names=['pos', 'strand', 'num reads', "5' het", 'seq'])
    avg_abundance = df['num reads'].mean()
    if "R1" in cluster_file:
        color = 'tomato'
        x_pos = 9.75
        replicate = "R1"
    elif "R2" in cluster_file:
        color = 'steelblue'
        x_pos = 10
        replicate = "R2"
    else:
        color = 'green'
        x_pos = 10
        replicate = "R3"
    
    # only add label the first time a replicate appears to avoid duplicates in the legend
    if replicate not in labels_added:
        label = replicate
        labels_added[replicate] = True
    else:
        label = None

    if "T1" in cluster_file:
        ax.scatter(x_pos, avg_abundance, c=color, label=label)
        if "R1" in cluster_file:
            R1_y_pos = avg_abundance - 5
            R1_x_pos = 8.5
        elif "R2" in cluster_file:
            R2_y_pos = avg_abundance + 5
            R2_x_pos = 10.5
        else:
            R3_y_pos = avg_abundance

    elif "T2" in cluster_file:
        ax.scatter(20, avg_abundance, c=color, label=label)
        line1_x = np.linspace(10, 20)
    else:
        ax.scatter(30, avg_abundance, c=color, label=label)


plt.text(9, R1_y_pos, "R1", fontsize = 8)
plt.text(9, R2_y_pos, "R2", fontsize = 8)
plt.text(9, R3_y_pos, "R3", fontsize = 8)


ax.legend()

line1_x = np.linspace(

    for idx, pos in enumerate(df['pos']):
        if df['strand'].iloc[idx] == '+':
            ax.scatter(df["5' het"].iloc[idx], df['num reads'].iloc[idx], c='steelblue')
        else:
            ax.scatter(df["5' het"].iloc[idx], df['num reads'].iloc[idx], c='tomato')

plt.xlabel('Five Prime Heterogeneity')
plt.ylabel('Abundance (hundreds of thousands)')
plt.yticks([100000, 200000, 300000, 400000], [1, 2, 3, 4])

for idx, pos in enumerate(df['pos']):
    if df['num reads'].iloc[idx] > 19000 and (df["5' het"].iloc[idx] < 0.78 or df["5' het"].iloc[idx] > 0.79):
        plt.text(df["5' het"].iloc[idx], df['num reads'].iloc[idx], pos, fontsize=7)
"""
plt.savefig('abundance_vs_time_scatter_26.pdf')
