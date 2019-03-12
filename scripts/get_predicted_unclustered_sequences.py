#!/usr/bin/env python

import re
import sys
import os

def get_seqid_dict_from_unclustered_fasta(unclustered_fastafileName):
	unclustered_fasta_seqid_dict={}
	unclustered_fastafile=open(unclustered_fastafileName, "r")
	for line in unclustered_fastafile:
		line = line.rstrip()
		if(re.match(r'^>',line)):
			linearr=re.split(r'\s+',line)
			seqid = linearr[0]
			seqid = seqid[1:]
			unclustered_fasta_seqid_dict[seqid]=1

	unclustered_fastafile.close()
	return(unclustered_fasta_seqid_dict)


def count_assigned_unclustered_sequences(under_cl_results_dirName, unclustered_fasta_seqid_dict, outfile_prefix):
	predicted_missing_seqs_outfile = open(under_cl_results_dirName+"/"+outfile_prefix+".predicted_missing_seqs","w")
	predicted_missing_seqs_count_outfile = open(under_cl_results_dirName+"/"+outfile_prefix+".predicted_missing_seqs_count","w")
	
	predicted_missing_seqs_count_outfile.write("famid"+" "+"predicted_missing_seqs_count"+"\n")

	for results_dir in os.listdir(under_cl_results_dirName):
		if not (os.path.isdir(under_cl_results_dirName+"/"+results_dir)):
			continue
		missing_sequences_results_fileName=under_cl_results_dirName+"/"+results_dir+"/"+results_dir+".missing_sequences"
		if not (os.path.exists(missing_sequences_results_fileName)):
			continue
		#if(os.stat(missing_sequences_results_fileName).st_size == 0):
		#	continue
	
		assigned_unclustered_sequences_count=0
		missing_sequences_results_file = open(missing_sequences_results_fileName, "r")
		for line in missing_sequences_results_file:
			line = line.rstrip()
			if(unclustered_fasta_seqid_dict.has_key(line)):
				assigned_unclustered_sequences_count+=1
				#print '{0} {1}'.format(results_dir, line)
				predicted_missing_seqs_outfile.write(results_dir+" "+line+"\n")
			
		missing_sequences_results_file.close()

		#print '{0} {1}'.format(results_dir, assigned_unclustered_sequences_count)
		predicted_missing_seqs_count_outfile.write(results_dir+" "+str(assigned_unclustered_sequences_count)+"\n")

	predicted_missing_seqs_outfile.close()
	predicted_missing_seqs_count_outfile.close()
	

################################################################################################################
unclustered_fastafileName = "/data/orthofinder_legumes/01_proteomes_legumes/orthofinder_legumes_unclustered_size_1-8.fa"
under_cl_results_dirName = "/data/family_quality_orthofinder_legumes/orthofinder_legumes/"
outfile_prefix = "orthofinder_legumes"

unclustered_fasta_seqid_dict = get_seqid_dict_from_unclustered_fasta(unclustered_fastafileName)
count_assigned_unclustered_sequences(under_cl_results_dirName, unclustered_fasta_seqid_dict, outfile_prefix)
