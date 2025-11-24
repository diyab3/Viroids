import sys
import os

usage = 'Usage:' + sys.argv[0] + ' <sam files>'

if len(sys.argv) < 2:
    print(usage)
    exit()

for input_file in sys.argv[1:]:
    root,ext = os.path.splitext(input_file)
    forward_outfile_name = root + "_forward" + ext
    reverse_outfile_name = root + "_reverse" + ext
    with open(forward_outfile_name, 'w') as forward_outfile:
        with open(reverse_outfile_name, 'w') as reverse_outfile:
            with open(input_file, 'r') as infile:
                for line in infile:
                    if line.startswith('@'):
                        forward_outfile.write(line)
                        reverse_outfile.write(line)
                    else:
                        terms = line.strip().split('\t')
                        if terms[2] != '*':
                            sam_flag = int(terms[1])
                            if sam_flag == 16:
                                reverse_outfile.write(line)
                            elif sam_flag == 0:
                                forward_outfile.write(line)
                            else:
                                print("Sam flag indicates read is not strictly aligned to forward or reverse strand")
                                exit()

