#!/bin/env python
import sys
import vcf
import pysam
from collections import namedtuple

input_vcf = sys.argv[1]
input_lift_over_file = sys.argv[2]
output_prefix = sys.argv[3]
bam_file = sys.argv[4]

# To do
# 1. Take care of input multiallelic sites


def read_in_lift_over_mapping(input_lift_over_file):
    mapping = {}
    Position = namedtuple('Position',['new_chrom','new_pos','new_ref_base','og_pos','og_ref_base'])
    with open(input_lift_over_file,'r') as input_liftOver_fh:
        for line in input_liftOver_fh:
            line = line.rstrip().split('\t')
            line[3] = int(line[3]) -1
            position = Position._make(line)
            #position.og_pos = int(position.og_pos) - 1
            if len(position.og_ref_base) != 1:
                continue
            if len(position.new_ref_base) != 1:
                continue
            mapping[position.og_pos] = position
    return mapping

def get_ref_mismatches(lift_over_mapping):
    ref_mismatches = []
    for pos in lift_over_mapping:
        position = lift_over_mapping[pos]        
        if position.new_ref_base != position.og_ref_base:
            ref_mismatches.append(pos)
    return ref_mismatches

def lift_ref_base_match(record,lift_over_mapping):
    info = record.INFO
    #info["igh_pos"] = record.POS
    record.INFO = info
    position = lift_over_mapping[record.POS]
    record.CHROM = position.new_chrom
    record.POS = int(position.new_pos)
    return record

def lift_ref_base_mismatch(record,lift_over_mapping):    
    position = lift_over_mapping[record.POS]
    if len(record.ALT) == 1:
        if record.ALT[0] == position.new_ref_base: # IGH ALT base is equal to new reference base
            het_sites = record.get_hets()
            if len(het_sites) != 0: # position is heterozygous
                record.ALT = [vcf.model._Substitution(position.og_ref_base)] #position.og_ref_base
            else:
                # In IGH, SNP is homozygous ALT. However, the ALT is the reference base in the
                # new reference. Therefore when lifted over, the SNV would be homozygous REF
                # and not a SNV
                return None
        else:
            # In IGH, the REF and ALT is not equal to the REF or ALT in the new reference
            # Therefore, the ALT is still a SNV in the new reference
            het_sites = record.get_hets()
            if len(het_sites) != 0: # SNV is multiallelic in new reference
                record.ALT = [vcf.model._Substitution(record.REF),vcf.model._Substitution(record.ALT[0])]
            else:
                # Site is homozygous
                # Unnecessary but complete
                #record.ALT = record.ALT            
                pass
    else:
        if position.new_ref_base in record.ALT:
            if record.ALT[0] == position.new_ref_base:
                alt_base = record.ALT[1]
            else:
                alt_base = record.ALT[0]
            record.ALT = [vcf.model._Substitution(alt_base)] 
            CallData = namedtuple('CallData',['GT','GQ'])
            call = vcf.model._Call(record,"sample",CallData(GT="0/1",GQ="60"))
            record.samples = [call]
        else:
            pass
    info = record.INFO
    #info["igh_pos"] = record.POS
    record.INFO = info    
    record.CHROM = position.new_chrom
    record.POS = int(position.new_pos)
    record.REF = position.new_ref_base
    return record

def assembly_location(read_name):
    read_origin = read_name.split("_")[0].split('=')[1]
    chrom = read_origin.split(":")[0]
    start = int(read_origin.split(":")[1].split("-")[0])
    end = int(read_origin.split(":")[1].split("-")[1])
    return [chrom,start,end]

def is_overlapping(a, b):
    if a[0] != b[0]:
        return False
    overlapping = False
    num_overlapping = max(0, min(a[2], b[2]) - max(a[1], b[1]))
    if num_overlapping > 0:
        overlapping = True
    return overlapping

def get_haplotype(read_name):
    return read_name.split("_")[1].split('=')[1]
        
def snps_per_hap(bamfile,position,ref_base):
    samfile = pysam.AlignmentFile(bamfile)
    regions = {}
    for read in samfile.fetch('igh',position,position + 1):
        if read.is_unmapped:
            continue
        if read.is_supplementary:
            continue
        if read.is_secondary:
            continue
        assembled_region = assembly_location(read.query_name)
        mapped_chrom = samfile.get_reference_name(read.reference_id)
        mapped_start = int(read.reference_start)
        mapped_end = int(read.reference_end)
        mapped_region = [mapped_chrom,mapped_start,mapped_end]
        if not is_overlapping(assembled_region,mapped_region):
            continue
        haplotype = get_haplotype(read.query_name)
        for read_pos, ref_pos in read.get_aligned_pairs():
            if read_pos == None:
                continue
            if ref_pos == None:
                continue
            if ref_pos != position:
                continue
            read_base = read.query_sequence[read_pos].upper()
            if ref_base != read_base:
                read_qual = read.query_qualities[read_pos]
                if mapped_chrom not in regions:
                    regions[mapped_chrom] = {}
                if ref_pos not in regions[mapped_chrom]:
                    regions[mapped_chrom][ref_pos] = {}
                regions[mapped_chrom][ref_pos][haplotype] = (ref_base,read_base,read_qual,read.query_name)
    return regions

def get_new_snp(bamfile,position,lift_over_mapping):
    ref_base = lift_over_mapping[position].new_ref_base
    ref_chrom = lift_over_mapping[position].new_chrom
    ref_pos = int(lift_over_mapping[position].new_pos)
    haplotype_snps = snps_per_hap(bamfile,position,ref_base)
    if len(haplotype_snps) == 0:
        return None
    if "0" in haplotype_snps["igh"][position]:
        genotype = "1/1"
        alt_bases = [haplotype_snps["igh"][position]["0"][1]]
    else:
        if len(haplotype_snps["igh"][position]) == 1:
            genotype = "0/1"
            if "1" in haplotype_snps["igh"][position]:
                alt_bases = [haplotype_snps["igh"][position]["1"][1]]
            else:
                alt_bases = [haplotype_snps["igh"][position]["2"][1]]                
        else:
            if haplotype_snps["igh"][position]["1"][1] == haplotype_snps["igh"][position]["2"][1]:
                genotype = "1/1"
                alt_bases = [haplotype_snps["igh"][position]["2"][1]]
            else:
                genotype = "1/2"
                alt_bases = [haplotype_snps["igh"][position]["1"][1],haplotype_snps["igh"][position]["2"][1]]
    contig_names = []
    for hap in haplotype_snps["igh"][position]:
        contig_names.append(haplotype_snps["igh"][position][hap][-1])
    # info = { "contig": contig_names,
    #          "region": ".",
    #          "read_support": "Yes",
    #          "intronic": ".",
    #          "LP1": ".",
    #          "RSS": ".",
    #          "gene": ".",
    #          "igh_region": ".",
    #          "phased_genotype": ".",
    #          "haplotype_block": ".",
    #          "sv_filter": ["No"],
    #          "igh_pos": position}
    info = {"read_support": "Yes",
            "contig": contig_names,
            "SV": ["No"],
            "intronic": ".",
            "LP1": ".",
            "RSS": ".",
            "gene": ".",
            "VDJ": ".",
            "read_genotype": "." }
    substitutions = []
    for alt_base in alt_bases:        
        substitution = vcf.model._Substitution(alt_base)
        substitutions.append(substitution)
    new_record = vcf.model._Record(ref_chrom,ref_pos,".",ref_base,substitutions,"60","PASS",info,"GT:GQ",[])
    CallData = namedtuple('CallData',['GT','GQ'])
    call = vcf.model._Call(new_record,"sample",CallData(GT=genotype,GQ="60"))
    new_record.samples.append(call)
    return new_record

def do_lift_over(input_vcf,lift_over_mapping,ref_mismatches,bamfile):
    vcf_reader = vcf.Reader(filename=input_vcf)
    lifted_records = []
    unlifted_records = []
    record_positions = []
    for record in vcf_reader:
        if record.POS in lift_over_mapping:
            record_positions.append(record.POS)
            if record.POS not in ref_mismatches:
                changed_record = lift_ref_base_match(record,lift_over_mapping)
            else:
                changed_record = lift_ref_base_mismatch(record,lift_over_mapping)
            if changed_record != None:
                lifted_records.append(changed_record)
        else:
            unlifted_records.append(record)
    for position in ref_mismatches:
        if position in record_positions: # Already taken care of
            continue
        record = get_new_snp(bamfile,position,lift_over_mapping)
        if record != None:
            lifted_records.append(record)
    return lifted_records,unlifted_records

def write_to_vcf(positions,output_vcf,input_vcf):
    vcf_reader = vcf.Reader(filename=input_vcf)
    vcf_writer = vcf.Writer(open(output_vcf, 'w'), vcf_reader)
    records = sorted(positions,key=lambda x: (x.CHROM,x.POS))
    for record in records:
        record.POS = int(record.POS) - 1 # added
        vcf_writer.write_record(record)
    vcf_writer.close()    

lift_over_mapping = read_in_lift_over_mapping(input_lift_over_file)
ref_mismatches = get_ref_mismatches(lift_over_mapping)
lifted_position,unlifted_records = do_lift_over(input_vcf,lift_over_mapping,ref_mismatches,bam_file)
lifted_vcf = "%s_lifted.vcf" % output_prefix
not_lifted_vcf = "%s_not_lifted.vcf" % output_prefix
write_to_vcf(lifted_position,lifted_vcf,input_vcf)
write_to_vcf(unlifted_records,not_lifted_vcf,input_vcf)



