# RPF_simulator
RPF_simulator is a user-friendly tool to simulate ribosome protected fragments(ribo-seq data) with different periodicity in ORFs.
RPF_simulator is based on python3
## Dependencies
* Numpy
* Biopython
## Installation
Install from source:
```
git clone https://www.github.com/GennongChen/RPF_simulator
```
## Tutorial
1. Input files  
  * Genome (fasta)  
  * Gencode annotation (gtf)  
  * Configration file  
2. Parameters  
  * "-m": minimum length of AA length  
  * "-s": substitution rate  
  * "-s": indel rate  
  * "-e": sequencing error rate  
3. Output files
  * Out file (contain ORF information, raw simulated reads and reads with snv)  
  * Fasta file (simulated reads with snv and sequencing error)  
