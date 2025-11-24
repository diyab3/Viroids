import sys

usage = 'Usage: ' + sys.argv[0] + ' <lane info file>'

if len(sys.argv) != 2:
    print(usage)
    exit()

lane_info_file = sys.argv[1]

with open(lane_info_file) as F:
    for line in F:
        sra_acc, name = line.strip().split()
        command = f'hqsub "fasterq-dump {sra_acc} -v" -q "*" -r {sra_acc}.log'
        print(command)

