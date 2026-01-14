import sys
import os

### MAIN FUNCTION ###
def main():
    usage = 'Usage:' + sys.argv[0] + ' <sam file(s)>'

    if len(sys.argv) < 2:
        print(usage)
        exit()

    sam_files = sys.argv[1:]
    
    for sam in sam_files:
        fixed = os.path.splitext(sam)[0] + "_fix.sam"
        #print(fixed)
        
        with open(sam, 'r') as infile:
            with open(fixed, 'w') as outfile:
                for line in infile:
                    if line.startswith('@'):
                        outfile.write(line)
                    else:
                        terms = line.strip().split('\t')
                        pos = int(terms[3])
                        if terms[2] != '*':
                           pos += 50
                           terms[3] = str(pos)
                        outfile.write('\t'.join(terms) + '\n')

### SUBROUTINES ###

if __name__ == "__main__":
    main()
