
VarCal is python script to call SNPs. It takes BAM and reference fasta as input to call variants, which are reported in a tabula format.

It has the following columns:
1. Chr - chromosome
2. POS - position of the variant
3. REF - Reference base
4. ALT - ALT base
5. Avg_BQ	- Average base quality
6. Avg_MQ	- Average mapping quality
7. Coverage	- Read coverage
8. Avg_PiR - Average position of the base in the read
9. ALT_DP	- Depth of alterate allele
10. AF - Allele frequency
11. Entropy	- Entropy based on number of +ve and -ve strands
12. Total_Pos_Strand - Number of total reads mapped to +ve strand
13. Total_Neg_Strand - Number of total reads mapped to -ve strand

The following filters are applied by default:
1. Uniquely mapped reads only
2. Minimum total coverage of 5
3. Minimum variant-supporting coverage of 2
4. Minimum average base quality of q15
5. Minimum variant allele frequency of 0.01

```
$ conda env create -f env.yml 
$ conda activate varcall
$ python3 varcall.py -h
usage: varcall.py [-h] -b BAM -r REF -o OUT
Generate variant calls from a BAM file using ref fasta
optional arguments:
  -h, --help         show this help message and exit
  -b BAM, --bam BAM  bam file
  -r REF, --ref REF  reference fasta
  -o OUT, --out OUT  file name to output variant calls
```

