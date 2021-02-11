#!/bin/env python
import sys

variants_infn = sys.argv[1]
mapping_infn = sys.argv[2]
lifted_variants_outfn = sys.argv[3]
imprecise_lifted_variants_outfn = sys.argv[4]
not_lifted_variants_outfn = sys.argv[5]

def read_in_liftOver_mapping(input_liftOver_file):
    mapping = {}
    with open(input_liftOver_file,'r') as input_liftOver_fh:
        for line in input_liftOver_fh:
            line = line.rstrip().split('\t')
            to_chrom = line[0]
            to_pos = line[1]
            to_base = line[2]
            from_pos = int(line[3]) - 1
            from_base = line[4]
            if len(from_base) != 1:
                continue
            if len(to_base) != 1:
                continue
            mapping[from_pos] = [from_base,to_chrom,to_pos,to_base]
    return mapping

def write_lifted_indel(mapping,start,end,stats,lifted_fh):
    out_chrom = mapping[start][1]
    out_start = mapping[start][2]
    out_end = mapping[end][2]
    out_line = [out_chrom,out_start,out_end] + stats
    lifted_fh.write("%s\n" % "\t".join(map(str,out_line)))

def write_imprecise_lifted_indel(mapping,start,end,stats,imprecise_lifted_fh):
    if end == None:
        out_chrom = mapping[start][1]
        out_start = mapping[start][2]
        out_end = int(out_start) + 1
    elif start == None:
        out_chrom = mapping[end][1]
        out_end = mapping[end][2]
        out_start = int(out_end) - 1
    out_line = [out_chrom,out_start,out_end] + stats
    imprecise_lifted_fh.write("%s\n" % "\t".join(map(str,out_line)))
    
    
liftOver_mapping = read_in_liftOver_mapping(mapping_infn)

lifted_variants_fh = open(lifted_variants_outfn,'w')
imprecise_lifted_variants_fh = open(imprecise_lifted_variants_outfn,'w') 
not_lifted_variants_fh = open(not_lifted_variants_outfn, 'w')

with open(variants_infn,'r') as variants_infh:
    for line in variants_infh:
        line = line.rstrip().split('\t')
        start = int(line[1])
        end = int(line[2])
        stats = line[3:]
        if start in liftOver_mapping and end in liftOver_mapping:
            write_lifted_indel(liftOver_mapping,start,end,stats,lifted_variants_fh)
        elif start in liftOver_mapping:
            write_imprecise_lifted_indel(liftOver_mapping,start,None,stats,imprecise_lifted_variants_fh)
        elif end in liftOver_mapping:
            write_imprecise_lifted_indel(liftOver_mapping,None,end,stats,imprecise_lifted_variants_fh)
        else:
            not_lifted_variants_fh.write("%s\n" % "\t".join(map(str,line)))

