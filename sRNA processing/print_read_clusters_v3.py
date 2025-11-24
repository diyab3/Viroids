# with pysam

#!/usr/bin/python3.9
import sys
# import pysam


def main():
    usage = 'Usage:' + sys.argv[0] + ' <sam files>'

    if len(sys.argv) < 2:
        print(usage)
        exit()


    sam_files = sys.argv[1:]

    #print("reading sam files, building distinct reads")
    dist_reads = {}
    for sam_file in sam_files:
        with open(sam_file) as F:
            for line in F:
                if not line.startswith('@'):
                    terms = line.strip().split('\t')
                    start = int(terms[3])
                    read_seq = terms[9]
                    length = len(read_seq)
                    end = start + length
                    sam_flag = int(terms[1])
                    if sam_flag >= 16:
                        strand = '-'
                    else:
                        strand = '+'
                    if terms[2] != '*':
                        if (start,end,strand,read_seq) in dist_reads.keys():
                            dist_reads[(start,end,strand,read_seq)] += 1
                        else:
                            dist_reads[(start,end,strand,read_seq)] = 1


    # Sort by value in descending order
    sorted_items_desc = sorted(dist_reads.items(), key=lambda item: item[1], reverse=True)

    # Keep track of clusters
    product_cluster = {}
    for read_data in sorted_items_desc:
        read_key,read_count = read_data
        start,end,strand,read_seq = read_key
        if strand == '-':
            rel_start = end
        else:
            rel_start = start
        updated_read = (rel_start,start,end,strand,read_count,read_seq)

        min_dist, min_cluster_key = find_min_dist(product_cluster,rel_start,strand)
        
        if min_dist <= 5:
            product_cluster[min_cluster_key].append(updated_read)
        else:
            cluster_key = (rel_start,strand)
            product_cluster[cluster_key] = [updated_read]

    for cluster_key in product_cluster:
        cluster_rel_start,cluster_strand = cluster_key
        top_read_seq = get_top_read_seq(product_cluster[cluster_key])
        total_reads = get_total_reads(product_cluster[cluster_key])
        five_prime_het = get_five_prime_het(product_cluster[cluster_key],cluster_rel_start)
        print(f"{cluster_rel_start}\t{cluster_strand}\t{total_reads}\t{five_prime_het}\t{top_read_seq}")

###############
# SUBROUTINES #
###############

def get_top_read_seq(product_cluster):
    top_read_data = product_cluster[0]
    rel_start,start,end,strand,read_count,read_seq = top_read_data
    return read_seq


def find_min_dist(product_cluster,rel_start,strand):

    min_dist = 1000
    min_cluster_key = ''
    for cluster_key in product_cluster:
        cluster_rel_start,cluster_strand = cluster_key
        if cluster_strand == strand:
            dist = abs(rel_start - cluster_rel_start)
            if dist < min_dist:
                min_dist = dist
                min_cluster_key = cluster_key
    '''
    if len(dist_list) > 1 and dist_list[0][0] == dist_list[1][0]:
        print("Error: read is equidistant to two clusters")
        exit()
    '''
    return min_dist,min_cluster_key


def get_total_reads(product_cluster):
    total_reads = 0
    for read in product_cluster:
        read_rel_start,read_start,read_end,read_strand,read_count,read_seq = read
        total_reads += read_count
    return total_reads

def get_five_prime_het(product_cluster,cluster_rel_start):
    total_reads = 0
    offset_reads = 0
    for read in product_cluster:
        read_rel_start,read_start,read_end,read_strand,read_count,read_seq = read
        if read_rel_start != cluster_rel_start:
            offset_reads += read_count
        total_reads += read_count
    return offset_reads / total_reads



main()


"""
samfile = pysam.AlignmentFile(sys.argv[1], "r")
count_dict = {}

for read in samfile.fetch():
    if read.is_mapped():
        start = read.reference_start
        end = read.reference_end
        if (start,end) in count_dict.keys():
            count_dict[(start,end)] += 1
        else:
            count_dict[(start,end)] = 1
samfile.close()

print(count_dict)
"""
