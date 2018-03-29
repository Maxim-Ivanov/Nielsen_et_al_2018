# This script modified bedGraph files, so that a record covering more than one basepair is converted to multiple records each covering a single basepair;
# E.g. input record:
# chr1 1000 1003 5
# Output records:
# chr1 1000 1001 5
# chr1 1001 1002 5
# chr1 1002 1003 5
# This increases the size of bedGraph files but prevents CAGEfightR from crashing on variableStep BigWig files;

import os, sys

input_dir = sys.argv[1]
filenames = [f for f in os.listdir(input_dir) if os.path.isfile(os.path.join(input_dir, f)) and (f.endswith('bg') or f.endswith('bedgraph'))]
for filename in filenames:
    print(filename)
    output_filename = filename[:filename.rfind('.')] + '_expanded.bg'
    with open(os.path.join(input_dir, filename), 'r') as input_file, open(os.path.join(input_dir, output_filename), 'w') as output_file:
        count = 0
        for line in input_file:
            if line.startswith('track'):
                output_file.write(line)
                continue
            fields = line.rstrip('\n').split('\t')
            start, end = int(fields[1]), int(fields[2])
            width = end - start
            if width > 1:
                count += 1
                for i in range(width):
                    new_start = start + i
                    new_end = new_start + 1
                    output_line = '\t'.join([fields[0], str(new_start), str(new_end), fields[3]])
                    output_file.write(output_line + '\n')
            else:
                output_file.write(line)
    print(count, 'lines were expanded\n;')
print('Done!')