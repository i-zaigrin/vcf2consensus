# vcf2consensus
Create consensus sequences from reference genome and VCF.

```bash
usage: vcf2consensus.py [-h] -v VCF -g GENOME -l LENGTH -c COUNT -f FREQ [-a]
                        [-d] [-n N] [-o OUTPUT]

optional arguments:
  -h, --help  show this help message and exit
  -v VCF      VCF file [required]
  -g GENOME   Genome in fasta format [required]
  -l LENGTH   Consensus length [required, >=10]
  -c COUNT    Consensus count [required, >=1]
  -f FREQ     Min allele frequence [required, 0<x<=1]
  -a          Consensus must contain at least 1 snp from vcf (can be reference
              for some samples) [default: on]
  -d          Use if vcf-header doesnt contain list of chromosomes and their
              length, or if chromosomes order is different from fasta file
              [default: off]
  -n N        Split consensuses into rows with length N [default: off, >=60]
  -o OUTPUT   Output name [default: output.fa]
```
