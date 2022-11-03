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
### 1. Prepare transcript and ORF index  
This step is used to generate ORF information and is used to subsequently generate simulation data.  
```
mkdir -p $out_dir && cd $out_dir
python /home/chengennong/code-manual/vscode/ORFfinding/simulator/prepare_transcripts.py \
    -g ${gtf} \
    -f ${fa} \
    -o ${out_dir}
```
#### Parameters  
  **-f/--fasta**  
    Human/Mouse genome sequences file. (fasta format)  
  **-g/--gtf**  
    Human/Mous gencode annotation file with CCDS tag. (gtf format)  
  **-o/--out_dir**  
    Annotation directory name. This directory is necessary for subsequent step.  
  
### 2. Fill in the configration file

### 3. Simulate RPF reads with SNV and sequencing error
This step mainly contains three parts:  
  Step1: Produce P-site tracks form the beginning to the distant position (determined by configration) and generate their corresponding coverage area anchored by a specific read length and offset.  
  Step2: Give each P-site track a CCDS ORF (filtered by parameters and transcript attribute), simualte substitution and indel (determined by parameters) to the corresponding transcript sequence and extract reads by corrected position with or without 5'utr indel.  
  Step3: Simualte sequencing error (determined by parameters) to simulated RPF reads generated in Step2 and output the information of simualtion data to a txt file and simulated reads to a fasta file.
```
python /home/chengennong/code-manual/vscode/ORFfinding/simulator/orf_feature.py \
    -c ${confi} \
    -a ${out_dir} \
    -o ${out_file} \
    -f ${T1_fasta} \
    -s 0.001 \
    -i 0.0001 \
    -e 0.005 \
    -m 7 > $T1_log_file
```

  * Configration file  
2. Parameters  
  * -m: minimum length of AA length  
  * -s: substitution rate  
  * -i: indel rate  
  * -e: sequencing error rate  
3. Output files
  * Out file (contain ORF information, raw simulated reads and reads with snv)  
  * Fasta file (simulated reads with snv and sequencing error)  
