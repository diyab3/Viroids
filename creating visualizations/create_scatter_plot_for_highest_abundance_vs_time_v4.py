import sys, re
#import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

usage = 'Usage: ' + sys.argv[0] + ' < clustering output files > '

if len(sys.argv) < 1:
    print(usage)
    exit()

cluster_files = sys.argv[1:]

data = {}
highest_six_pos = []
for cluster_file in cluster_files:
    # normalized_Cusa_sRNAseq_T3_HSVd_R3_HSVd_clusters.sh
    match = re.search('normalized_Cusa_sRNAseq_(.*)_HSVd_(.*)_HSVd_clusters',cluster_file)
    time = match.group(1)
    rep = match.group(2)

    # find loci to track from T3
    if 'T3' in time:
        with open(cluster_file) as F:
            #print(rep)
            for line in F:
                if not line.startswith('pos'):
                    pos,strand,num_reads,het5,seq,rpm = line.strip().split('\t')
                    #print(pos,strand,num_reads,het5,seq,rpm)
                    rpm = float(rpm)
                    pos = int(pos)
                    key = (pos,strand)

                    existing_key = None
                    for k in data:
                        existing_pos, existing_strand = k
                        if existing_strand == strand and abs(existing_pos - pos) < 5:
                            existing_key = k
                            break

                    if existing_key is None:
                        if len(highest_six_pos) < 6:
                            if (pos,strand) not in data:
                                data[(pos,strand)] = {}
                            if time not in data[(pos,strand)]: 
                                data[(pos,strand)][time] = []
                            data[(pos,strand)][time].append(rpm)
                            highest_six_pos.append((pos,strand))

                        else:
                            lowest_key = None
                            lowest_rpm = float('inf')
                            
                            for k in data:
                                if time in data[k]:
                                    current_rpm = data[k][time][-1]
                                else:
                                    current_rpm = 0

                                if current_rpm < lowest_rpm:
                                    lowest_rpm = current_rpm
                                    lowest_key = k

                            if rpm > lowest_rpm:
                                del data[lowest_key]
                                highest_six_pos.remove((lowest_key[0], lowest_key[1]))

                                data[key] = {time : [rpm]}
                                highest_six_pos.append((pos,strand))

                    else:
                        if existing_key not in data:
                            data[existing_key] = {}
                        if time not in data[existing_key]:
                            data[existing_key][time] = []
                        data[existing_key][time].append(rpm)

                    #print(data)
                    #print(highest_six_pos)

#print(data)
#print(highest_six_pos)
for cluster_file in cluster_files:
    match = re.search('normalized_Cusa_sRNAseq_(.*)_HSVd_(.*)_HSVd_clusters',cluster_file)
    time = match.group(1)
    rep = match.group(2)

    if 'T1' in time or 'T2' in time:
        with open(cluster_file) as F:
            #print(rep)
            for line in F:
                if not line.startswith('pos'):
                    pos,strand,num_reads,het5,seq,rpm = line.strip().split('\t')
                    pos = int(pos)
                    rpm = float(rpm)
                    for posn,forward_or_reverse in highest_six_pos:
                        if pos == posn and forward_or_reverse == strand:
                            if time not in data[(pos,strand)]:
                                data[(pos,strand)][time] = []
                            data[(pos,strand)][time].append(rpm)
            #print(data)
            #print(highest_six_pos)
"""                                     
print(data)
print(highest_six_pos)
for cluster_file in cluster_files:
    with open(cluster_file) as F:
        for line in F:
            if not line.startswith('pos'):
                pos,strand,num_reads,het5,seq,rpm = line.strip().split('\t')
                if (pos,strand) not in data:
                    data[(pos,strand)] = {}
                if time not in data[(pos,strand)]:
                    data[(pos,strand)][time] = []
                data[(pos,strand)][time].append(float(rpm))
"""
time_map = {'T1':10, 'T2':17, 'T3':24}
times = [10,17,24]

n_plus = 0
n_minus = 0
colors_plus = ['#c51b7d', '#de77ae', '#f1b6da']
colors_minus = ['#b8e186', '#7fbc41', '#4d9221']
for n,(pos,strand) in enumerate(data):
    plot_data = []
    if strand == '+':
        c = colors_plus[n_plus]
        n_plus += 1
    else:
        c = colors_minus[n_minus]
        n_minus += 1

    for time in data[(pos,strand)]:

        #print(pos,strand,time,np.mean(data[(pos,strand)][time]))
        plot_data.append((time_map[time], np.mean(data[(pos,strand)][time])))
        
    plot_data.sort(key=lambda x:x[0])
    x_vals,y_vals = zip(*plot_data)
    x_values = list(x_vals)
    y_values = list(y_vals)
    
    plt.plot(x_values,y_values,color=c)
    plt.scatter(x_values,y_values,color=c, zorder = 2, label=f'{pos} {strand}')
        
        
plt.legend()




    
"""

exit()


fig,ax = plt.subplots()
plt.xticks([10, 20, 30], ["T1", "T2", "T3"])
plt.xlim(0, 40)

labels_added = {}

T1_df = pd.DataFrame()
T2_df = pd.DataFrame()
T3_df = pd.DataFrame()

for cluster_file in cluster_files:
    replicate_df = pd.read_csv(cluster_file, sep="\t", header=None, names=['pos', 'strand', 'num reads', "5' het", 'seq'])
    if "T1" in cluster_file:
        T1_df = pd.concat([T1_df, replicate_df], ignore_index=True)
    elif "T2" in cluster_file:
        T2_df = pd.concat([T2_df, replicate_df], ignore_index=True)
    else:
        T3_df = pd.concat([T3_df, replicate_df], ignore_index=True)


T1_sorted_abundance = T1_df.sort_values(by='num reads', ascending=False)
T1_top_six = T1_sorted_abundance.head(6)
for idx,num_reads in enumerate(T1_top_six['num reads']):
    ax.scatter(10, num_reads)


T2_sorted_abundance = T2_df.sort_values(by='num reads', ascending=False)
T2_top_six = T2_sorted_abundance.head(6)
for idx,num_reads in enumerate(T2_top_six['num reads']):
    ax.scatter(20, num_reads)

T3_sorted_abundance = T3_df.sort_values(by='num reads', ascending=False)
T3_top_six = T3_sorted_abundance.head(6)
for idx,num_reads in enumerate(T3_top_six['num reads']):
    ax.scatter(30, num_reads)


"""
plt.savefig("highest_abundance_vs_time_scatter_HSVd.pdf")

