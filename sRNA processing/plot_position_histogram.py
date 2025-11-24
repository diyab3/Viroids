import sys
import numpy as np
import matplotlib.pyplot as plt


usage = 'Usage:' + sys.argv[0] + ' <sam file>'

if len(sys.argv) != 2:
    print(usage)
    exit()

sam_file = sys.argv[1]
seqname = 'Y09352.1'

data = []
with open(sam_file) as F:
    for line in F:
        if not line.startswith('@'):
            terms = line.strip().split('\t')
            if terms[2] == seqname:
                read = terms[9]
                length = len(read)
                data.append(int(terms[3])+length)

bin_size = 1
bins = np.arange(min(data),max(data)+bin_size,bin_size)
plt.hist(data, bins=bins)
plt.savefig('test.png')
