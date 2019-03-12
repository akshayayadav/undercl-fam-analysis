#!/usr/bin/env python

from __future__ import division
import re
import sys
import os

#small_fam - families in the unclustered.fa file
#large_fam - families processed by under-clustering correction and detection method

##########################################################################################################################################

def fam_fasta_dicts(fam_fasta_dirName):
	seqid_famid_dict={}
	famid_famsize_dict={}
	for fam_fasta_fileName in os.listdir(fam_fasta_dirName):
		fam_fasta_file = open(fam_fasta_dirName+"/"+fam_fasta_fileName,"r")
		famid_famsize_dict[fam_fasta_fileName]=0
		for line in fam_fasta_file:
			line = line.rstrip()
			if not (re.match(r'^>',line)):
				continue
			seqid_famid_dict[line[1:]]=fam_fasta_fileName
			famid_famsize_dict[fam_fasta_fileName]+=1
		fam_fasta_file.close()

	return([seqid_famid_dict, famid_famsize_dict])		

def get_large_fam_small_fam_seqcount_dict(reassigned_uncl_seqs_results_fileName, seqid_famid_dict):
	large_fam_small_fam_seqcount_dict={}
	reassigned_uncl_seqs_results_file = open(reassigned_uncl_seqs_results_fileName,"r")
	for line in reassigned_uncl_seqs_results_file:
		line = line.rstrip()
		linearr = re.split(r'\s+',line)
		large_famid = linearr[0]
		seqid = linearr[1]
		small_famid = seqid_famid_dict[seqid]
		if(large_fam_small_fam_seqcount_dict.has_key(large_famid)):
			if(large_fam_small_fam_seqcount_dict[large_famid].has_key(small_famid)):
				large_fam_small_fam_seqcount_dict[large_famid][small_famid]+=1
			else:
				large_fam_small_fam_seqcount_dict[large_famid][small_famid]=1
		else:
			large_fam_small_fam_seqcount_dict[large_famid]={}
			large_fam_small_fam_seqcount_dict[large_famid][small_famid]=1
		
	return(large_fam_small_fam_seqcount_dict)

def get_family_mergings(large_fam_small_fam_seqcount_dict, famid_famsize_dict, small_fam_large_fam_overlap_cutoff):
	for large_fam in large_fam_small_fam_seqcount_dict:
		for small_fam in large_fam_small_fam_seqcount_dict[large_fam]:
			small_fam_seqcount=large_fam_small_fam_seqcount_dict[large_fam][small_fam]
			small_fam_size = famid_famsize_dict[small_fam]
			small_fam_large_fam_overlap = small_fam_seqcount/small_fam_size
			if (small_fam_large_fam_overlap>small_fam_large_fam_overlap_cutoff):
				print '{0} {1} {2} {3}'.format(large_fam, small_fam, small_fam_large_fam_overlap, small_fam_size)
##########################################################################################################################################

fam_fasta_dirName="/home/aayadav/research/orthofinder_eudicot/orthofinder_family_fasta/"
#reassigned_uncl_seqs_results_fileName="orthofinder_eudicot_famsize_13_22/orthofinder_eudicot_famsize_13_22.reassigned_uncl_seqs"
reassigned_uncl_seqs_results_fileName="orthofinder_eudicot_famsize_5_12/orthofinder_eudicot_famsize_5_12.reassigned_uncl_seqs"
small_fam_large_fam_overlap_cutoff=0.5

seqid_famid_dict, famid_famsize_dict = fam_fasta_dicts(fam_fasta_dirName)
large_fam_small_fam_seqcount_dict = get_large_fam_small_fam_seqcount_dict(reassigned_uncl_seqs_results_fileName, seqid_famid_dict)
get_family_mergings(large_fam_small_fam_seqcount_dict, famid_famsize_dict, small_fam_large_fam_overlap_cutoff)

