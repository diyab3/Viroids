import sys
import os

### MAIN FUNCTION ###
def main():
    usage = 'Usage:' + sys.argv[0] + '<HSVd sam files>' + ' -- ' + '<shift sam files>'

    if '--' not in sys.argv:
        print(usage)
        exit()

    split_idx = sys.argv.index('--')
    HSVd_sam_files = sys.argv[1:split_idx]
    shift_sam_files = sys.argv[split_idx + 1:]

    for HSVd_sam, shift_sam in zip(HSVd_sam_files, shift_sam_files):
        output = os.path.splitext(HSVd_sam)[0] + "_merge.sam"
        #print(output)
        print(f'samtools merge -o {output} {HSVd_sam} {shift_sam}')


### SUBROUTINES ###

if __name__ == "__main__":
    main()
