[![Review Assignment Due Date](https://classroom.github.com/assets/deadline-readme-button-22041afd0340ce965d47ae6ef1cefeee28c7c493a6346c4f15d667ab976d596c.svg)](https://classroom.github.com/a/Fx_NP1KI)
# MicroDNA
Using k-mers to search strings


## Smith-Waterman
### Usage
```
usage: kmer_idx.py [-h] [--kmer KMER] --reference REFERENCE [--read_size READ SIZE] --num_reads NUMBER OF READS [--error_rate ERROR RATE] [--max_mismatches MAX MISMATCHES] [--experiment EXPERIMENT]

optional arguments:
  -h, --help                        show this help message and exit
  --kmer KMER                       Kmer size (default=3)
  --error_rate ERROR_RATE           Error rate (default=0.01)
  --max_mismatches MAX_MISMATCHES   Maximum number of mismatches (default=3)
  --experiment EXPERIMENT           Run k-mer size and error rate experiments
```
### Example
```
$ python3 src/kmer_idx.py 
    --reference data/chr22.fa.gz 
    --num_reads 5 
    --experiment

```
<center><img src="alignment_score_vs_error.png" width="600"/></center>
<center><img src="alignments_vs_kmer.png" width="600"/></center>
<center><img src="runtime_vs_kmer.png" width="600"/></center>

