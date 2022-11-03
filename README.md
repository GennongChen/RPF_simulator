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
```
mkdir -p $out_dir && cd $out_dir
python /home/chengennong/code-manual/vscode/ORFfinding/simulator/prepare_transcripts.py \
    -g ${gtf} \
    -f ${fa} \
    -o ${out_dir}
```
#### Parameters  
  **-f/--fasta**  
    &emsp;&emsp;Human/Mouse genome sequences file. (fasta format)  
  **-g/--gtf**  
    &emsp;&emsp;Human/Mous gencode annotation file with CCDS tag. (gtf format)  
  **-o/--out_dir**  
    &emsp;&emsp;Annotation directory name. This directory is necessary for subsequent step.  

This step is used to generate ORF information and is used to subsequently generate simulation data.  

### 2. Fill in the configration file
```
reads_num: 10000000
CCDS_ORF_num: 10000

read_length: [25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35]
offset: [12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]
periodicty: [[3,3,3],[5,3,2],[6,3,1],[7,2,1],[8,1,1],[9,1,0],[8,1,1],[7,2,1],[6,3,1],[5,3,2],[3,3,3]]
predicted_status: [0,0,1,1,1,1,1,1,1,0,0]
read_length_proportion: [2.35, 2.99, 5.16, 8.93, 13.68, 17.8, 15.2, 15.2, 10.1, 4.74, 3.85]
```
Note: The sum of read_length_proportion must be equal to 100.
Note: CCDS_ORF_num * read_length_proportion_k / 100 must be a integer.

### 3. Simulate RPF reads with SNV and sequencing error
```
python /home/chengennong/code-manual/vscode/ORFfinding/simulator/orf_feature.py \
    -c ${config} \
    -a ${out_dir} \
    -o ${output_file} \
    -f ${output_fasta} \
    -s 0.001 \
    -i 0.0001 \
    -e 0.005 \
    -m 7 > ${log}
```
#### Parameters  
  **-c/--config_file**  
    &emsp;&emsp;Configration file.  
  **-a/--annot_dir**  
    &emsp;&emsp;Annotation directory name generated by prepare_transcripts.py.  
  **-m/--min_aa_length**  
    &emsp;&emsp;Minimum length of AA sequence.  
  **-s/--sub_rate**  
    &emsp;&emsp;Substitution rate of transcript sequence.  
  **-i/--ind_rate**  
    &emsp;&emsp;Indel rate of transcript sequence of transcript sequence.  
  **-e/--err_rate**  
    &emsp;&emsp;Sequencing error rate of simulated RPF reads.  
  **-o/--output_file**  
    &emsp;&emsp;The information of simualtion data.(txt format)  
  **-f/--output_fasta**  
    &emsp;&emsp;Final simulated reads. (fasta format)  

#### Output files
  **output_file**  
    &emsp;&emsp;Contain ORF information, raw simulated reads and reads with snv)  
  **output_fastae**  
    &emsp;&emsp;Simulated reads with snv and sequencing error)  

This step mainly contains three parts:  
    &emsp;&emsp;Step1: Produce P-site tracks form the beginning to the distant position (determined by configration) and generate their corresponding coverage area anchored by a specific read length and offset.  
    &emsp;&emsp;Step2: Give each P-site track a CCDS ORF (filtered by parameters and transcript attribute), simualte substitution and indel (determined by parameters) to the corresponding transcript sequence and extract reads by corrected position with or without 5'utr indel.  
    &emsp;&emsp;Step3: Simualte sequencing error (determined by parameters) to simulated RPF reads generated in Step2 and output the information of simualtion data to a txt file and simulated reads to a fasta file.
