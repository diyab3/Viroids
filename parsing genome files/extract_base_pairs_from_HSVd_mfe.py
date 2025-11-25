import sys

usage = 'Usage: ' + sys.argv[0] + ' <mfe (dot-bracket) file> '

if len(sys.argv) < 2:
    print(usage)
    exit()

mfe_file = sys.argv[1]
with open(mfe_file, 'r') as mfe:
    dot_bracket_str = ''
    for line in mfe:
        if line.startswith('.') or line.startswith('('):
            dot_bracket_str += line.split()[0]

stack = []
base_pairs = []
print(len(dot_bracket_str))
for idx in range(len(dot_bracket_str)):
    if dot_bracket_str[idx] == '(':
        stack.append(idx)
    if dot_bracket_str[idx] == ')':
        bp_first_half = stack.pop()
        base_pairs.append((bp_first_half, idx))

#print(base_pairs)

