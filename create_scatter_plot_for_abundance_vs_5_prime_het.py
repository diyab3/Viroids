import sys
import pandas as pd
import matplotlib.pyplot as plt

usage = 'Usage: ' + sys.argv[0] + ' < clustering output file > '

if len(sys.argv) < 1:
    print(usage)
    exit()

clusters = sys.argv[1]

df = pd.read_csv(clusters, sep="\t", header=None, names=['pos', 'strand', 'num reads', "5' het", 'seq'])

fig,ax = plt.subplots()

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

plt.savefig('abundance_vs_5_prime_het_scatter_12.pdf')

