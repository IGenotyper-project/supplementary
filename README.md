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
