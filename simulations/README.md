## Extension of partial transcripts using simulated reads
### Data sets
1. Reads with constant coverage
2. Reads with variable coverage

### Reference used:
CYP12C1 (NM_001300180.1)
CYP4G1 (NM_080292.4)
CYP9B2 (NM_078922.4)

### Reads simulator
https://github.com/zstephens/neat-genreads

### Dataset 1
```
python ./neat-genreads/genReads.py -r CYP12C1.fasta -R 100 -o constant_cov --pe 250 30 -c 50
python ./neat-genreads/genReads.py -r CYP4G1.fasta -R 100 -o constant_cov --pe 250 30 -c 50
python ./neat-genreads/genReads.py -r CYP9B2.fasta -R 100 -o constant_cov --pe 250 30 -c 50

```
### Dataset2
```
python ./neat-genreads/genReads.py -r CYP12C1.fasta -R 100 -o one --pe 250 30 -c 10
python ./neat-genreads/genReads.py -r CYP12C1.fasta -R 100 -o two --pe 250 30 -c 100 -t cyp12c1.bed
cat two.read1.fastq >> one.read1.fastq
cat two.read2.fastq >> two.read2.fastq
```
```
python ./neat-genreads/genReads.py -r CYP4G1.fasta -R 100 -o one --pe 250 30 -c 10
python ./neat-genreads/genReads.py -r CYP4G1.fasta -R 100 -o two --pe 250 30 -c 100 -t cyp4g1.bed
cat two.read1.fastq >> one.read1.fastq
cat two.read2.fastq >> two.read2.fastq
```
```
python ./neat-genreads/genReads.py -r CYP9B2.fasta -R 100 -o one --pe 250 30 -c 10
python ./neat-genreads/genReads.py -r CYP9B2.fasta -R 100 -o two --pe 250 30 -c 100 -t cyp9b2.bed
cat two.read1.fastq >> one.read1.fastq
cat two.read2.fastq >> two.read2.fastq
```
