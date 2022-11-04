#stat ref
#spe=Human
#gtf=/home/chengennong/ref/Human/gencode.v34.annotation.gtf
#grep ccdsid $gtf|grep _NF|head -n30
#grep ccdsid $gtf|grep -v _NF|cut -f3|head -n30
#grep ccdsid $gtf|grep mRNA_start_NF|cut -f3|head -n30
#grep ccdsid $gtf|grep cds_start_NF|awk '$3!="transcript"'|cut -f1-6
###Tx num
#cat $gtf|awk '$3=="transcript"'|wc -l
###ccds Tx num
#grep ccdsid $gtf|awk '$3=="transcript"'|wc -l
###complete ccds Tx num 
#grep ccdsid $gtf|grep -v "mRNA_start_NF\|cds_start_NF\|mRNA_stop_NF\|cds_stop_NF" |awk '$3=="transcript"'|wc -l

#annotation
##Human
out_dir=/home/chengennong/test/benchmark_simulator/annotation/Human
gtf=/home/chengennong/ref/Human/gencode.v34.annotation.gtf
fa=/home/chengennong/ref/Human/GRCh38.genome.fa
mkdir -p $out_dir && cd $out_dir
#python /home/chengennong/code-manual/vscode/ORFfinding/simulator/prepare_transcripts.py \
#    -g $gtf \
#    -f $fa \
#    -o $out_dir
##Mouse
out_dir=/home/chengennong/test/benchmark_simulator/annotation/Mouse
gtf=/home/chengennong/ref/Mouse/gencode.vM25.annotation.gtf
fa=/home/chengennong/ref/Mouse/GRCm38.genome.fa
#mkdir -p $out_dir && cd $out_dir
#python /home/chengennong/code-manual/vscode/ORFfinding/simulator/prepare_transcripts.py \
#    -g $gtf \
#    -f $fa \
#    -o $out_dir

for spe in Human Mouse
do
    out_dir=/home/chengennong/test/benchmark_simulator/annotation/${spe}
    T1_out_file=/home/chengennong/test/benchmark_simulator/out/${spe}.T1.txt
    T2_out_file=/home/chengennong/test/benchmark_simulator/out/${spe}.T2.txt
    #out_fasta=/home/chengennong/test/benchmark_simulator/out/${spe}.fasta
    config=/home/chengennong/test/benchmark_simulator/out/${spe}_config.txt
    T1_fasta=/home/chengennong/test/benchmark_simulator/out/${spe}.T1.fasta
    T2_fasta=/home/chengennong/test/benchmark_simulator/out/${spe}.T2.fasta
    T1_log_file=${T1_fasta}.log
    T2_log_file=${T2_fasta}.log
    python /home/chengennong/code-manual/vscode/ORFfinding/simulator/orf_feature.py \
        -c $config \
        -a $out_dir \
        -o ${T1_out_file} \
        -f ${T1_fasta} \
        -s 0.001 \
        -i 0.0001 \
        -e 0.005 \
        -m 7 > $T1_log_file && \
    python /home/chengennong/code-manual/vscode/ORFfinding/simulator/orf_feature.py \
        -c $config \
        -a $out_dir \
        -o ${T2_out_file} \
        -f ${T2_fasta} \
        -s 0.005 \
        -i 0.002 \
        -e 0.02 \
        -m 7 > $T2_log_file
done
#add mutation and error
#T1: 0.001 substitution, 0.0001 indel 0.005, error
#T2: 0.005 substitution, 0.002 indel, 0.01 error
#T3: 0.03 substitution, 0.005 indel, 0.02 error

#cat Human.T1.fasta | sed ":a;N;s/\n//"|grep -e "^[A-Z]"|head
