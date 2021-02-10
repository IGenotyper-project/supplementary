#!/bin/env python
import pysam
from Bio import SeqIO
from collections import Counter

def get_matches(query_entries,database_entries):
    matches = {}
    for query_entry in query_entries:
        query_seq = query_entries[query_entry].seq.upper()
        matches[(query_entry,query_seq)] = set()
        for database_entry in database_entries:
            database_seq = database_entries[database_entry].seq.upper()
            if query_seq.count(database_seq) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
            if query_seq.count(database_seq.reverse_complement()) > 0:
                matches[(query_entry,query_seq)].add(database_entry)
    return matches

def query_gene_hap(entry):
    gene_name = entry.split("_")[0].split("=")[1]
    hap = entry.split("_")[1].split("=")[1]
    return (gene_name,hap)

def hit_gene_allele(entry):
    hit_gene = entry.split("_")[0].split("=")[1]
    hit_allele = int(entry.split("_")[1].split("=")[1])
    return (hit_gene,hit_allele)

def genotype_genes(seq_fa,db):
    query_entries = SeqIO.to_dict(SeqIO.parse(seq_fa,"fasta"))
    db_entries = SeqIO.to_dict(SeqIO.parse(db,"fasta"))
    matches = get_matches(query_entries,db_entries) 
    alleles = {}
    for query_gene, query_seq in matches:
        gene_name,hap = query_gene_hap(query_gene)
        added = False
        alleles[(gene_name,hap)] = []
        for hit in matches[(query_gene,query_seq)]:
            hit_gene, hit_allele = hit_gene_allele(hit)
            if hit_gene != gene_name:
                continue
            alleles[(gene_name,hap)].append(("MATCH",hit_allele))
            added = True
        if added:
            continue
        alleles[(gene_name,hap)].append(("NOVEL",query_seq))
    return alleles
            
def write_genotypes(alleles,fn):
    header = ["gene","hap","category","allele","obs_count"]    
    with open(fn,'w') as fh:
        fh.write("%s\n" % "\t".join(header))
        for gene_name,hap in alleles:
            all_alleles = alleles[(gene_name,hap)]
            allele_counts = Counter(all_alleles)
            for allele_count in allele_counts:                
                cat,allele = allele_count
                output = [gene_name,hap,cat,allele,allele_counts[allele_count]]
                fh.write("%s\n" % "\t".join(map(str,output)))

def extract_sequence_from(read,chrom,start,end):
    read_start = None
    read_end = None
    aligned_pairs = read.get_aligned_pairs()
    for query_pos, ref_pos in aligned_pairs:
        if query_pos == None:
            continue
        if ref_pos == None:
            continue
        if ref_pos <= start:
            read_start = query_pos
        read_end = query_pos
        if ref_pos > end:
            break
    if read_start == None:
        return ""
    if read_end == None:
        return ""
    return read.query_sequence[read_start:read_end].upper()

def load_bed_regions(bedfile,add_fourth=False):
    bed_regions = []
    with open(bedfile,'r') as bedfh:
        for line in bedfh:
            line = line.rstrip().split('\t')
            chrom = line[0]
            start = int(line[1])
            end = int(line[2])
            if add_fourth:
                annotation = True
                if len(line) == 4:
                    annotation = line[3]
                bed_regions.append([chrom,start,end,annotation])
            else:
                bed_regions.append([chrom,start,end])
    return bed_regions

def extract_sequence(phased_bam,bed,fasta):
    records = []
    coords = load_bed_regions(bed,True)
    samfile = pysam.AlignmentFile(phased_bam)
    for chrom,start,end,feat in coords:
        i = 0
        for read in samfile.fetch(chrom,start,end):
            if skip_read(read):
                continue
            if read.has_tag("RG"):
                hap = read.get_tag("RG",True)[0]
            else:
                hap = "None"
            name = "feat=%s_hap=%s_i=%s" % (feat,hap,i)
            seq = extract_sequence_from(read,chrom,start,end)
            if len(seq) == 0:
                continue
            record = SeqRecord(Seq(seq),id=name,name=name,description="")
            records.append(record)
            i += 1
    SeqIO.write(records,fasta,"fasta")


bam = sys.argv[1]
bed = sys.argv[2]
db = sys.argv[3]
gene_outfn = sys.argv[4]
assignment_outfn = sys.argv[5]

extract_sequence(bam,bed,gene_outfn)
alleles = genotype_genes(gene_outfn,db)
write_genotypes(alleles,assignment_outfn)
