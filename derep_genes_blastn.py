import os, sys, re, subprocess, shlex, argparse
from collections import defaultdict
from Bio import SeqIO

args_parser = argparse.ArgumentParser(description="Script for dereplicating genes", epilog="Virginia Tech Department of Biological Sciences")
args_parser.add_argument('-i','--input_seq', required=True, help='Input sequences to dereplicate')
args_parser.add_argument('-o','--out_prefix', required=True, help='Prefix for output files')
args_parser.add_argument('-s','--seq_type', required=True, help='Type of input sequences, option are "nucl" or "prot", for nucleotide and amino acid sequences, respectively')
args_parser.add_argument('-id','--id_min', required=False, default = str(100) ,help='Minimum percent identity (out of 100) between two aligned sequences')
args_parser.add_argument('-c','--cov_min', required=False, default = str(1),help='Minimum proportion of alignment length of at least one of the sequences. Default is 1. Ex: A sequence of 560 bp may align over its entire length to a seq of 780 bp and its proportion would be 1')

args_parser = args_parser.parse_args()

input_seq = args_parser.input_seq
out_prefix = args_parser.out_prefix
seq_type = args_parser.seq_type
id_min = args_parser.id_min
cov_min = args_parser.cov_min
id_min = int(id_min)
cov_min = float(cov_min)

seq2len = {}
seq_open = open(input_seq,'r')

# Get sequence length table
for record in SeqIO.parse(seq_open,'fasta'):
	name = record.id
	length = len(record.seq)
	seq2len[name] = length

# Makeblastdb
blastdb_name = out_prefix + "_blastdb"
cmd = "makeblastdb -in " + input_seq + " -out " + blastdb_name + " -dbtype " + seq_type + " -parse_seqids"
print(cmd)
cmd2 = shlex.split(cmd)
subprocess.call(cmd2, stdout=open("test_out", "w"), stderr=open("test_err", "w"))

# run blastn
blastout = out_prefix + "_blast_tbl.tsv"
blast_cmd = str()
if seq_type == "nucl":
	blast_cmd = "blastn -query " + input_seq + " -out " + blastout + " -db " + blastdb_name + " -evalue 0.001 -outfmt 6"
elif seq_type == "prot":
	blast_cmd = "blastp -query " + input_seq + " -out " + blastout + " -db " + blastdb_name + " -evalue 0.001 -outfmt 6"
print(blast_cmd)
blast_cmd2 = shlex.split(blast_cmd)
subprocess.call(blast_cmd2, stdout=open("test_out", "w"), stderr=open("test_err", "w"))

# Parse blastn for reps 

seq2hits = defaultdict(list)

seq2rep = {}

all_seq = []
redund_seq = []

blastout_open = open(blastout)

for i in blastout_open:
	line = i.rstrip()
	tabs = line.split("\t")
	seq = tabs[0]
	hit = tabs[1]
	aln_len = int(tabs[3])
	ident = float(tabs[2])
	seqlen = seq2len[seq]
	hitlen = seq2len[hit]
	all_seq.append(seq)
	all_seq.append(hit)
	if seq == hit:
		pass
	else:
		if ident >= id_min:
			prop_seq = aln_len / seqlen
			prop_hit = aln_len / hitlen
			if prop_seq >= cov_min or prop_hit >= cov_min:
				infolist = [seq,hit,str(prop_hit),str(prop_seq)]
				info = "\t".join(infolist)
				max_seq = str()
				if seqlen > hitlen:
					max_seq = seq
				elif seqlen < hitlen:
					max_seq = hit
				elif seqlen == hitlen:
					max_seq = seq
				seq2hits[max_seq].append(hit)
				redund_seq.append(hit)
				redund_seq.append(seq)

seq_set = set(all_seq)

for seq in seq_set:
	if seq not in redund_seq:
		seq2rep[seq] = seq

for key, values in seq2hits.items():
	max_seq = key
	hits = values
	numhits = len(hits)
	max_len = seq2len[max_seq]
	hit_lens = []
	for hit in hits:
		seq2rep[hit] = max_seq
	seq2rep[max_seq] = max_seq

rep_tbl_name = out_prefix + "_rep_tbl.tsv"
rep_tbl_open = open(rep_tbl_name,'w')

rep_seq = []

#print("sequence" + "\t" + "representative" + "\t" + "seq_len" + "\t" + "rep_len")
rep_tbl_open.write("sequence" + "\t" + "representative" + "\t" + "seq_len" + "\t" + "rep_len" + "\n")
for key, values in seq2rep.items():
	seqlen = seq2len[key]
	replen = seq2len[values]
	#print(key + "\t" + values + "\t" + str(seqlen) + "\t" + str(replen))
	rep_seq.append(values)
	rep_tbl_open.write(key + "\t" + values + "\t" + str(seqlen) + "\t" + str(replen) + "\n")

seq_open_2 = open(input_seq,'r')

rep_set = set(rep_seq)

outseq_name = out_prefix + "_rep_seq.fna"
outseq_open = open(outseq_name,'w')
keepseq = []

for record in SeqIO.parse(seq_open_2,'fasta'):
	if record.id in rep_set:
		keepseq.append(record)

SeqIO.write(keepseq,outseq_open,'fasta')
