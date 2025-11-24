import re
import sys

usage = 'Usage: ' + sys.argv[0] + ' <flagstat output files> '

if len(sys.argv) < 1:
    print(usage)
    exit()

text_files = sys.argv[1:]
pattern = "+ 0 mapped"
total_reads = 0

for text_file in text_files:
    with open(text_file, 'r') as text:
        for line in text:
            if pattern in line:
                values = line.split()
                total_reads += int(values[0])

print(total_reads)
