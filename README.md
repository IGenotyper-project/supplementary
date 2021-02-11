# supplementary
## parse IMGT database
```
python parse_imgt.py <IMGT_DB_FASTA> > parsed_imgt.fasta
```

## Compare IGenotyper assembly and CCS allele
```
python compare_assembly_ccs_alleles.py <assembly_alleles.txt> <ccs_alleles.txt> <outdir>
```
  
## Assign overlapping gene sequence to allele
```
python assign_alleles.py  <bam alignment> <gene coords bed> <parsed imgt db> <extracted gene seq outfn> <assignment outfn>
```

## LiftOver igh indels and SVs to hg38
```
python lift_over/lift_over_non_snps.py \
  igh_variants.bed \
  lift_over/igh_to_hg38_liftOver.txt \
  hg38_lifted.bed \
  hg38_imprecise.bed \
  igh_to_hg38_not_lifted.bed
```
