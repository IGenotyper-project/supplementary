#!/bin/env python
import sys

input_indel = sys.argv[1]
input_liftOver_file = sys.argv[2]
lifted_indels = sys.argv[3]
imprecise_lifted_indels = sys.argv[4]
not_lifted_indels = sys.argv[5]

def read_in_liftOver_mapping(input_liftOver_file):
    mapping = {}
    with open(input_liftOver_file,'r') as input_liftOver_fh:
        for line in input_liftOver_fh:
            line = line.rstrip().split('\t')
            to_lift_to_chrom = line[0]
            to_lift_to_pos = line[1]
            to_lift_to_base = line[2]
            lift_from_pos = int(line[3]) - 1
            lift_from_base = line[4]
            if len(lift_from_base) != 1:
                continue
            if len(to_lift_to_base) != 1:
                continue
            mapping[lift_from_pos] = [lift_from_base,to_lift_to_chrom,to_lift_to_pos,to_lift_to_base]
    return mapping

def write_lifted_indel(liftOver_mapping,in_start_pos,in_end_pos,indel_stats,lifted_indels_fn):
    out_chrom_pos = liftOver_mapping[in_start_pos][1]
    out_start_pos = liftOver_mapping[in_start_pos][2]
    out_end_pos = liftOver_mapping[in_end_pos][2]
    out_line = [out_chrom_pos,out_start_pos,out_end_pos] + indel_stats
    output_indel_fw.write("%s\n" % "\t".join(map(str,out_line)))

def write_imprecise_lifted_indel(liftOver_mapping,start,end,indel_stats,imprecise_lifted_indels_fh):
    if end == None:
        out_chrom_pos = liftOver_mapping[start][1]
        out_start_pos = liftOver_mapping[start][2]
        out_end_pos = int(out_start_pos) + 1
    elif start == None:
        out_chrom_pos = liftOver_mapping[end][1]
        out_end_pos = liftOver_mapping[end][2]
        out_start_pos = int(out_end_pos) - 1
    out_line = [out_chrom_pos,out_start_pos,out_end_pos] + indel_stats
    imprecise_lifted_indels_fh.write("%s\n" % "\t".join(map(str,out_line)))
    
    
liftOver_mapping = read_in_liftOver_mapping(input_liftOver_file)

lifted_indels_fn = open(lifted_indels,'w')
imprecise_lifted_indels_fh = open(imprecise_lifted_indels,'w') 
not_lifted_indels_fh = open(not_lifted_indels, 'w')

with open(input_indel,'r') as output_indel_fh:
    for line in output_indel_fh:
        line = line.rstrip().split('\t')
        in_start_pos = int(line[1])
        in_end_pos = int(line[2])
        indel_stats = line[3:]
        if in_start_pos in liftOver_mapping and in_end_pos in liftOver_mapping:
            write_lifted_indel(liftOver_mapping,in_start_pos,in_end_pos,indel_stats,lifted_indels_fn)                
        elif in_start_pos in liftOver_mapping:
            write_imprecise_lifted_indel(liftOver_mapping,in_start_pos,None,indel_stats,imprecise_lifted_indels_fh)
        elif in_end_pos in liftOver_mapping:
            write_imprecise_lifted_indel(liftOver_mapping,None,in_end_pos,indel_stats,imprecise_lifted_indels_fh)
        else:
            not_lifted_indels_fh.write("%s\n" % "\t".join(map(str,line)))

