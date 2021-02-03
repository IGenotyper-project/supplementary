#!/bin/env python
import sys

# File should be from IMGT_V-QUEST_reference_directory
infile = sys.argv[1]

seqs = []
with open(infile, 'r') as fh:
    for line in fh:
        line = line.rstrip()
        if ">" in line:
            if len(seqs) != 0:
                print(">%s" % name)
                print("%s" % "".join(seqs))
                seqs = []
            line = line.split('|')
            gene = line[1].split("*")[0]
            allele = line[1].split("*")[1]
            function = line[3]
            name = "gene=%s_allele=%s_function=%s" % (gene,allele,function)
        else:
            seqs.append(line.replace(".",''))
            
print(">%s" % name)
print("%s" % "".join(seqs))


