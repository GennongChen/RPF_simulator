#!/usr/bin/env python

from prepare_transcripts import *
#from ccds_orf_psite import *
from collections import OrderedDict, defaultdict
import numpy as np
import os
import sys
from Bio import SeqIO
from Bio.Seq import MutableSeq, Seq
from Bio.SeqRecord import SeqRecord
import argparse
import random
import yaml

def parsing_args():
    parser = argparse.ArgumentParser(
        description=
        "Designed for extract position feature by transcrip_id with chrom,\
        starnd, t or gcoords."
    )
    parser.add_argument( "-c", "--config_file",dest="config_file",required=True,type=str,help="simulation config file.")
    parser.add_argument( "-a", "--annot_dir", dest="annot_dir", required=True, type=str, help="transcripts annotation directory, generated by prepare_transcripts.")
    #parser.add_argument("-S","--start_codon",default="ATG",type=str,dest="start_codon",help="The canonical start codon. default: ATG")
    #parser.add_argument("-A","--alt_start_codons",default="",type=str,dest="alternative_start_codons",help="The alternative start codon, such as CTG,GTG, default: None. Multiple codons should be separated by comma.")
    parser.add_argument("-m","--min_aa_length",dest="min_aa_length",default="7",required=False,help="The minimal length of predicted peptides,default 7",type=int)
    parser.add_argument("-o","--output_file",dest="outname",default="final_result",help="output file name, default: final_result",type=str)
    parser.add_argument("-f","--output_fasta",dest="outfasta",default="final_fasta",help="output the fasta file of simulated reads",type=str)
    parser.add_argument("-s","--sub_rate",dest="read_sub_rate",required=True,help="read substitution rate",type=float)
    parser.add_argument("-i","--ind_rate",dest="read_indel_rate",required=True,help="read indel rate",type=float)
    parser.add_argument("-e","--err_rate",dest="read_error_rate",required=True,help="read sequencing error rate",type=float)    

    args = parser.parse_args()
    #if args.start_codon:
    #    args.start_codon = args.start_codon.strip().split(",")
    #if args.alternative_start_codons:
    #    args.alternative_start_codons = args.alternative_start_codons.strip(
    #    ).split(",")
    #else:
    #    args.alternative_start_codons = None
    if not os.path.exists(args.annot_dir):
        raise ValueError(
            "Error, the transcript annotation directory not found: {} \n \
                     pls run prepare_transcripts.py first. ".format(
                args.annot_dir))
    if not os.path.exists(args.config_file):
        raise ValueError(
            "Error, the configration file not found")
    return args

def ORF_record(only_ATG=False):
    keys = [
        "ORF_ID", "ORF_type", "transcript_id", "transcript_type", "transcript_length", "gene_id",
        "gene_name", "gene_type", "chrom", "strand", "ORF_length",
        "ORF_tstart", "ORF_tstop", "ORF_gstart", "ORF_gstop",
        "annotated_tstart", "annotated_tstop", "annotated_gstart",
        "annotated_gstop", "transcript_tag", "cutoff", "predicted_status",
        "offset", "read_length", "psite_codon_num", "psite_track", "read_loc", "read_seq", "read_muta_seq"
    ]
    if only_ATG is False:
        keys.insert(9, "start_codon")
    r = OrderedDict(zip(keys, [None] * len(keys)))
    return r

def read_config(config_file):
    config = open(config_file, 'r')
    config_var = yaml.load(config.read(), Loader=yaml.FullLoader)
    reads_num = int(config_var['reads_num']); CCDS_ORF_num = int(config_var['CCDS_ORF_num'])
    read_length = config_var['read_length']; offset = config_var['offset']; read_length_proportion = config_var['read_length_proportion']
    periodicty = config_var['periodicty']; predicted_status = config_var['predicted_status']
    config.close()
    print("Reading cofigration ", config_file); print("read_length offset periodicty read_length_proportion predicted_status reads_num_length ORF_num_length codon_num")
    psite_track_list = []; codon_num_list = []
    for ORF_num in range(len(read_length)):
        read_length_k = read_length[ORF_num];  offset_k = offset[ORF_num]; read_length_proportion_k = round(read_length_proportion[ORF_num],2)
        periodicty_k = periodicty[ORF_num]; predicted_status_k = str(predicted_status[ORF_num])
        reads_num_length_k = int(round(reads_num * read_length_proportion_k / 100, 4))
        ORF_num_length_k = int(round(CCDS_ORF_num * read_length_proportion_k / 100, 4))
        codon_num_k = int(round((reads_num * read_length_proportion_k / 100) / (CCDS_ORF_num * read_length_proportion_k / 100) / sum(periodicty_k),4))
        codon_num_list.append(codon_num_k)
        print(read_length_k,offset_k,periodicty_k,read_length_proportion_k,predicted_status_k,reads_num_length_k,ORF_num_length_k,codon_num_k)
        for n in range(ORF_num_length_k):
            psite_track_k = Psite_track(reads_num,CCDS_ORF_num,read_length_k,periodicty_k,offset_k,read_length_proportion_k,predicted_status_k)
            psite_track_list.append(psite_track_k)
    max_codon_num = max(codon_num_list)
    print("CCDS_ORF_num: ",CCDS_ORF_num,"\tmax_codon_num: ", max_codon_num, "\tpsite_track_num", len(psite_track_list))
    if round(sum(read_length_proportion), 4) != 100:
        raise ValueError("Error: The sum of read_length_proportion must be equal to 100!")
    if CCDS_ORF_num != len(psite_track_list):
        raise ValueError("Error: CCDS_ORF_num * read_length_proportion_k / 100 must be a integer!")
    return psite_track_list, CCDS_ORF_num, max_codon_num

def write_result(orf_results, outname):

    header = "\t".join(list(orf_results[0].keys())[:-1])
    with open(outname, "w") as fout:
        fout.write(header + "\n")
        for v in orf_results:
            #fout.write("\t".join(map(str, list(v.values())[:-1])))
            fout.write("\t".join(map(str, list(v.values())[:-1])) + "\n")

def classfy_orf(orfiv, cdsiv):
    """
    define the orf type according the orf position related to cds region
    """
    if orfiv.end == cdsiv.end or orfiv.end == cdsiv.end - 3:
        if orfiv.start == cdsiv.start:
            orftype = "ORF_annotated"
        else:
            orftype = "Others"
    return orftype

def alt_seq(mutable_seq, read_sub_rate, read_indel_rate):
    read_sub = np.random.binomial(1, read_sub_rate, len(mutable_seq))
    try: 
        mutable_seq_alt = MutableSeq(mutable_seq)
        for i in np.nditer(np.nonzero(read_sub)):
            mutable_seq_alt[int(i)] = random.sample("ATCG".replace(mutable_seq[int(i)],""), 1)[0]
        mutable_seq = mutable_seq_alt
    except: 
        pass
    read_indel = np.random.binomial(1, read_indel_rate, len(mutable_seq))
    try: 
        d = 0; i = 0; n_utr_base_var = 0
        for id in np.nditer(np.nonzero(read_indel)):
            id =int(id)
            if random.sample("01", 1)[0] == "0":
                id = id - d
                mutable_seq_l = list(mutable_seq)
                mutable_seq_l.pop(id)
                mutable_seq = "".join(mutable_seq_l)
                d += 1
                if id < start:
                    n_utr_base_var -= 1
            else:
                id = id + i
                mutable_seq_l = list(mutable_seq)
                mutable_seq_l.insert(id, random.sample("ATCG", 1)[0])
                mutable_seq = "".join(mutable_seq_l)
                i += 1
                if id < start:
                    n_utr_base_var += 1
    except: 
        pass
    return mutable_seq, n_utr_base_var

def add_sequencing_error_seq(mutable_seq, read_error_rate):
    read_error = np.random.binomial(1, read_error_rate, len(mutable_seq))
    try: 
        mutable_seq_alt = MutableSeq(mutable_seq)
        for i in np.nditer(np.nonzero(read_error)):
            mutable_seq_alt[int(i)] = random.sample("ATCGN".replace(mutable_seq[int(i)],""), 1)[0]
        mutable_seq = mutable_seq_alt
    except: 
        pass
    return mutable_seq

def orf_tlevel(annot_dir, transcript_dict, psite_track_list, CCDS_ORF_num, max_codon_num, min_aa_length, read_sub_rate, read_indel_rate):
    tid_list=[]
    candicate_orf_list=[]
    transcript_seq = GenomeSeq(
        os.path.join(annot_dir, "transcripts_sequence.fa"))
    random.shuffle(list(transcript_dict.keys()))
    n_orf = 0
    print("\nRandomly electing CCDS_ORF_num CCDS ORFs with complete mRNA, CDS tag and > 35nt 5'UTR\n")
    for tid in transcript_dict.keys():
        if 'tag' not in transcript_dict[tid].attr: continue #no tag
        if 'CCDS' not in transcript_dict[tid].attr.get('tag',[]): continue #no CCDS tag
        not_com_tag = ['mRNA_start_NF','cds_start_NF','mRNA_stop_NF','cds_stop_NF']
        if len(list(set(transcript_dict[tid].attr.get("tag",[])).intersection(set(not_com_tag)))) > 0: continue #ORF not complete
        tid_list.append(tid)
        start = transcript_dict[tid].cds.start
        end = transcript_dict[tid].cds.end
        strand = transcript_dict[tid].cds.strand
        length = transcript_dict[tid].length
        transcript_id = tid
        if length - start < 35: continue #5' UTR too short
        if length - end < 35: continue #3'UTR too short
        if end - start < min_aa_length * 3:
            print("Filtered: CDS length < min_aa_length * 3", tid, (end - start), " < ", min_aa_length * 3)
            continue #ORF shorter than threshold
        if n_orf < CCDS_ORF_num:
            psite_track_obj = psite_track_list[n_orf]
            codon_num = psite_track_obj.codon_num
            if (end - start) < max_codon_num * 3 + 35:
                print("Filtered: CDS length < max_codon_num * 3 + 35", transcript_id, (end - start), " < ", max_codon_num * 3 + 35)
                continue #ORF short than max codon number 
            offset = psite_track_obj.offset
            read_length = psite_track_obj.read_length
            cutoff = psite_track_obj.periodicty
            predicted_status = psite_track_obj.predicted_status
            psite_track = psite_track_obj.psite_track
            read_start=[[(n+start-offset)]*psite_track[n] for n in range(len(psite_track))]
            read_end=[[(n+start-offset+read_length)]*psite_track[n] for n in range(len(psite_track))]
            read_loc=[[str(n+start-offset + 1) + '-' + str(n+start-offset+read_length)]*psite_track[n] for n in range(len(psite_track))] #convert to 1-base
            read_seq=[[transcript_seq.get_seq(transcript_id, (n+start-offset),(n+start-offset+read_length), '+')]*psite_track[n] for n in range(len(psite_track))]
            #add snv
            mutable_seq, n_utr_base_var = alt_seq(transcript_seq.get_seq(transcript_id), read_sub_rate, read_indel_rate)
            read_muta_seq=[[str(mutable_seq[(n+start-offset+n_utr_base_var):(n+start-offset+read_length+n_utr_base_var)])]
                *psite_track[n] for n in range(len(psite_track))]

            candicate_orf_list.append(
                    Interval_source(int(start), int(end), strand, transcript_id,
                        cutoff, predicted_status, psite_track, offset, read_length, read_loc, codon_num, read_seq, read_muta_seq)
                    )
            n_orf += 1
            if n_orf % 1e3 == 0 or n_orf == CCDS_ORF_num:
                print("\t" + str(int(n_orf/1e3)) + "k CCDS ORF selected\n")
    if CCDS_ORF_num > n_orf:
        print("\n\tAll avaiable CCDS ORF number is  .....\n")
        exit(0)
    else:
        print("\n\tThe total number of complete CCDS ORF (5 UTR > 35nt): " + str(len(candicate_orf_list)) + "\n")
    return candicate_orf_list, tid_list

def orf_feature(candicate_orf_list, gene_dict, transcript_dict, annot_dir, #START_CODON, ALTERNATIVE_START_CODON_LIST,
                outname, outfasta, CCDS_ORF_num, read_error_rate, only_ATG=False):
    transcript_seq = GenomeSeq(
        os.path.join(annot_dir, "transcripts_sequence.fa"))
    orf_results = []
    seq_records = []
    n_orf = 0
    #random.shuffle(candicate_orf_list)
    for orf_iv in candicate_orf_list:
        tid = orf_iv.transcript_id
        tobj = transcript_dict[tid]
        gobj = gene_dict[tobj.gene_id]
        tseq = transcript_seq.get_seq(tobj.transcript_id)
        if n_orf < CCDS_ORF_num:
            n_orf += 1
            orf_dict = ORF_record(only_ATG)
            orf_dict["orf_iv"] = orf_iv
            orf_dict["transcript_id"] = tobj.transcript_id
            orf_dict["transcript_type"] = tobj.transcript_type
            orf_dict["transcript_length"] = tobj.length
            orf_dict["gene_id"] = tobj.gene_id
            orf_dict["gene_name"] = tobj.gene_name
            orf_dict["gene_type"] = gobj.gene_type
            orf_dict["chrom"] = tobj.chrom
            orf_dict["strand"] = tobj.genomic_iv.strand
            orf_dict["ORF_tstart"] = orf_iv.start + 1  # convert to 1-base
            orf_dict["ORF_tstop"] = orf_iv.end + 3  # include stop-codon
            if (orf_iv.end + 3) > tobj.length:
                orf_dict["ORF_tstop"] = tobj.length
            orf_gstart = transcript_pos_transform(
                tobj, orf_iv.start) + 1  # convert to 1-base
            orf_gstop = transcript_pos_transform(
                tobj, orf_iv.end + 2) + 1  # include stop-codon
            orf_dict["ORF_gstart"] = orf_gstart
            orf_dict["ORF_gstop"] = orf_gstop
            orf_dict["ORF_length"] = orf_iv.end - orf_iv.start
            if tobj.startcodon is not None:
                annotated_tstart = tobj.startcodon.start + 1
                annotated_gstart = transcript_pos_transform(
                    tobj, tobj.startcodon.start) + 1
            else:
                annotated_tstart = None
                annotated_gstart = None
            if tobj.stopcodon is not None:
                annotated_tstop = tobj.stopcodon.end
                annotated_gstop = transcript_pos_transform(
                    tobj, tobj.stopcodon.end - 1) + 1
            else:
                annotated_tstop = None
                annotated_gstop = None
            orf_dict["annotated_tstart"] = annotated_tstart
            orf_dict["annotated_tstop"] = annotated_tstop
            if annotated_tstop and annotated_tstop > tobj.length:
                orf_dict["annotated_tstop"] = tobj.length
            orf_dict["annotated_gstart"] = annotated_gstart
            orf_dict["annotated_gstop"] = annotated_gstop
            orf_dict["ORF_ID"] = "%s_%s_%s_%i" % (tobj.gene_id, orf_gstart,
                                                  orf_gstop, orf_iv.length / 3)
            orf_dict["transcript_tag"] = ";".join(tobj.attr["tag"] if "tag" in tobj.attr else ["No_tag_transcript"])
            orf_dict["cutoff"] = orf_iv.cutoff
            orf_dict["predicted_status"] = orf_iv.predicted_status
            orf_dict["offset"] = orf_iv.offset
            orf_dict["read_length"] = orf_iv.read_length
            orf_dict["psite_codon_num"] = orf_iv.codon_num
            orf_dict["psite_track"] = orf_iv.psite_track
            orf_dict["read_loc"] = orf_iv.read_loc
            orf_dict["read_seq"] = orf_iv.read_seq
            orf_dict["read_muta_seq"] = orf_iv.read_muta_seq
            orf_seq = transcript_seq.get_seq(tobj.transcript_id, orf_iv.start,
                                             orf_iv.end, "+")
            orf_dict["ORF_type"] = classfy_orf(orf_iv, tobj.cds)
            if not only_ATG:
                orf_dict["start_codon"] = orf_seq[:3]
            orf_results.append(orf_dict)

            for n_loc in range(len(orf_iv.read_loc)):
                for n_read in range(len(orf_iv.read_muta_seq[n_loc])):
                    if (orf_iv.read_muta_seq[n_loc][n_read]):
                        add_error_seq = Seq(str(add_sequencing_error_seq(orf_iv.read_muta_seq[n_loc][n_read],read_error_rate)))
                        seq_record = SeqRecord(add_error_seq)
                        seq_record.id = orf_dict["ORF_ID"] + "_" + tobj.chrom + "-" + tobj.transcript_id + ":" + orf_iv.read_loc[n_loc][n_read] + "-" + orf_iv.predicted_status + "-" + "_".join(map(str,orf_iv.cutoff)) + "-" + str(orf_iv.read_length) + "-" + str(orf_iv.offset) + "-" + str(n_loc) + "-" + str(n_read)
                        seq_record.description = ",".join(map(str,orf_iv.cutoff))
                        seq_records.append(seq_record)
        if n_orf % 1e3 == 0 or n_orf == CCDS_ORF_num:
            print(str(int(n_orf/1e3)) + "k CCDS ORF processed\n")
    print("\n\tThe number of ORF used in simulated data: " + str(len(orf_results)))

    if len(orf_results) == 0:
        sys.stderr.write("Opps, no translated ORFs were detected.\n")
        orf_dict = ORF_record(only_ATG)
        orf_results.append(orf_dict)
        header = "\t".join(list(orf_results[0].keys())[:-1])
        with open(outname, "w") as fout:
            fout.write(header + "\n")

    else:
        print("\n\tWriting the results to file .....\n")
        write_result(orf_results, outname)
        print("\n\t" + str(len(seq_records)) + " reads was in " + outfasta + "\n")
        SeqIO.write(seq_records, outfasta, "fasta")

def main():
    args = parsing_args()
    psite_track_list, CCDS_ORF_num, max_codon_num = read_config(args.config_file)
    gene_dict, transcript_dict = load_transcripts_pickle(
        os.path.join(args.annot_dir, "transcripts.pickle"))
    candicate_orf_list, tid_list = orf_tlevel(args.annot_dir, transcript_dict,
        psite_track_list = psite_track_list,
        CCDS_ORF_num = CCDS_ORF_num,
        max_codon_num = max_codon_num,
        min_aa_length=args.min_aa_length,
        read_sub_rate=args.read_sub_rate,read_indel_rate=args.read_indel_rate)
    orf_feature(candicate_orf_list=candicate_orf_list,
                gene_dict=gene_dict,
                transcript_dict=transcript_dict,
                annot_dir=args.annot_dir,
                #START_CODON=args.start_codon,
                #ALTERNATIVE_START_CODON_LIST=args.alternative_start_codons,
                outname=args.outname,
                outfasta=args.outfasta,
                CCDS_ORF_num = CCDS_ORF_num,
                read_error_rate=args.read_error_rate,
                only_ATG=False)

if __name__ == "__main__":
    main()