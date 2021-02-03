#!/bin/env python
import sys
import pandas as pd

assembly = sys.argv[1]
ccs = sys.argv[2]
outdir = sys.argv[3]

common = "%s/common.txt" % outdir
assembly_uniq = "%s/assembly_uniq.txt" % outdir
ccs_uniq = "%s/ccs_uniq.txt" % outdir

assembly = pd.read_csv(assembly,sep="\t")
ccs = pd.read_csv(ccs,sep="\t")

merged = assembly.merge(ccs,how='inner',on=["gene","hap","category","allele"],suffixes=('_assembly', '_ccs'))
all_ = assembly.merge(ccs,how='outer',on=["gene","hap","category","allele"],suffixes=('_assembly', '_ccs'))


ccs_unique = all_[all_["obs_count_assembly"].isnull()]
assembly_unique = all_[all_["obs_count_ccs"].isnull()]

merged.to_csv(common,sep="\t",index=False)
ccs_unique.to_csv(ccs_uniq,sep="\t",index=False)
assembly_unique.to_csv(assembly_uniq,sep="\t",index=False)
