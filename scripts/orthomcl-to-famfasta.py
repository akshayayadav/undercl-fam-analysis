#!/usr/bin/env python

import re
import sys
import os


def get_seqid_sequence_dict(fasta_fileName):
	fasta_file = open(fasta_fileName,"r")
	seqid_sequence_dict={}
	for line in fasta_file:
		line = line.rstrip()
		if(re.match(r'^>',line)):
			seqid = re.split("\s+",line)[0]
			seqid = seqid[1:]
		else:
			if(seqid_sequence_dict.has_key(seqid)):
				seqid_sequence_dict[seqid]=seqid_sequence_dict[seqid]+line
			else:
				seqid_sequence_dict[seqid]=line

	fasta_file.close()
	return(seqid_sequence_dict)

def read_orthomcl_file(orthomcl_fileName, seqid_sequence_dict, famfasta_outdirectory, min_fam_size, unclustered_seq_outfileName):
	orthomcl_file = open(orthomcl_fileName, "r")
	for line in orthomcl_file:
		line=line.rstrip()
		linearr = re.split(r'\s+',line)
		famid=linearr.pop(0)
		famid=famid[:len(famid)-1]
		seqid_arr = linearr
		if(len(seqid_arr)>min_fam_size):
			print_famfasta_file(famid, seqid_sequence_dict, seqid_arr, famfasta_outdirectory)
		else:
			print_unclustered_fasta(unclustered_seq_outfileName, seqid_sequence_dict, seqid_arr)

	orthomcl_file.close()

def print_famfasta_file(famid, seqid_sequence_dict, seqid_arr, famfasta_outdirectory):
	famfasta_outfile = open(famfasta_outdirectory+"/"+famid, "w")
	for seqid in seqid_arr:
		famfasta_outfile.write(">"+seqid+"\n"+seqid_sequence_dict[seqid]+"\n")
	
	famfasta_outfile.close()

def print_unclustered_fasta(unclustered_seq_outfileName, seqid_sequence_dict, seqid_arr):
	unclustered_seq_outfile = open(unclustered_seq_outfileName,"a")
	for seqid in seqid_arr:
		unclustered_seq_outfile.write(">"+seqid+"\n"+seqid_sequence_dict[seqid]+"\n")
	unclustered_seq_outfile.close()

################################################################
min_fam_size=0
fasta_fileName = "/home/aayadav/research/orthofinder_legumes/01_proteomes_legumes/orthofinder_legumes.fa"
orthomcl_fileName = "/home/aayadav/research/orthofinder_legumes/Results_Feb11/Orthogroups.txt"
famfasta_outdirectory = "/home/aayadav/research/orthofinder_legumes/orthofinder_legumes_family_fasta/"
unclustered_seq_outfileName = "/home/aayadav/research/orthofinder_legumes/01_proteomes_legumes/orthfinder_legumes_unclustered.fa"
seqid_sequence_dict = get_seqid_sequence_dict(fasta_fileName)
read_orthomcl_file(orthomcl_fileName, seqid_sequence_dict, famfasta_outdirectory, min_fam_size, unclustered_seq_outfileName)
