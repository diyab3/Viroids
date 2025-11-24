import sys

# Usage statement
usage = 'Usage: ' + sys.argv[0] + ' <fasta file of genome we mapped to>' + '<SAM files of reads>'

if len(sys.argv) < 3:
    print(usage)
    exit()

fasta = sys.argv[1]
sam = sys.argv[2:]

seq = ""
# get the length of the viroid genome
with open(fasta, 'r') as genome:
    for line in genome:
        if not line.startswith(">"):
            seq += line.strip()

seq_len = len(seq)

print("read_length,strand,positions,counts")  # print the column labels

# access the information we need from the sam files
for sam_file in sam:
    # intialize dictionaries for the forward and reverse strands
    plus_d = {}
    minus_d = {}

    min_length = 20
    max_length = 26
    # for each number in the range between minimum and maximum possible read length, create an entry in each dictionary
    for L in range(min_length, max_length):
        # the keys of the dictionary will be the length of the read
        plus_d[L] = [0] * seq_len
        # the values will be a list of the same length as the viroid genome, initialized to zeros
        minus_d[L] = [0] * seq_len

    with open(sam_file, 'r') as SAM:
        for line in SAM:  # for each line in the sam file
            terms = line.strip().split('\t')  # remove whitespace at the beginning and end of the line and split on tabs (columns are separated by tabs)
            if len(terms) >= 11:
                flag = terms[1]  # the flag will tell us which strand we are on
                pos = int(terms[3]) - 1  # SAM is 1-based so we need to subtract by 1 to get the 0-based start position of the read
                read_seq = terms[9] # the actual read sequence
                length = len(read_seq)  # the length of the read
                if 20 <= length < 26:
                    if terms[2] != '*':  # if the RNAME column is not '*', so the read was actually aligned to the genome
                        if int(flag) == 0:  # a SAM flag of 0 means the read was on the forward strand
                            plus_d[length][pos] += 1
                        elif int(flag) == 16:  # a SAM flag of 16 means the read was on the reverse strand
                            minus_d[length][pos] += 1
                        else:
                            print("Error: read not mapped to forward or reverse strand")
                            exit()


    # print the information out
    positions = list(range(seq_len))  # generate a list with all of the possible positions
    for key in plus_d:  # iterate over the forward strand dictionary
        for pos,count in zip(positions, plus_d[key]):  # create tuples of each posiiton and the count at that position, and iterate over them
            print(f"{key},+,{pos},{count}")

    for key in minus_d:  # iterate over the reverse strand dictionary
        for pos,count in zip(positions, minus_d[key]):
            print(f"{key},-,{pos},{count}")



