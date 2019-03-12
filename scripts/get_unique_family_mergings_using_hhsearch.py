#!/usr/bin/env python

import re
import sys
import os
import subprocess
import glob
from shutil import copyfile

#this program works on the output of get_family_mergings.py. This is used to detect single large families for small families that are detected as merging into more that one large families using the hhsearch program

def get_small_fam_large_fam_arr_dict(family_mergings_fileName):
	small_fam_large_fam_arr_dict = {}
	large_famid_small_famid_overlap_dict = {}
	small_famid_famsize_dict = {}
	family_mergings_file = open(family_mergings_fileName)
	for line in family_mergings_file:
		line = line.rstrip()
		if (re.match(r'^\#',line)):
			continue
		linearr = re.split(r'\s+',line)
		large_famid = linearr[0]
		small_famid = linearr[1]
		l_s_overlap = linearr[2]
		small_famid_famsize = linearr[3]
		small_famid_famsize_dict[small_famid] = int(small_famid_famsize)

		if(small_fam_large_fam_arr_dict.has_key(small_famid)):
			small_fam_large_fam_arr_dict[small_famid].append(large_famid)
		else:
			small_fam_large_fam_arr_dict[small_famid]=list()
			small_fam_large_fam_arr_dict[small_famid].append(large_famid)

		fill_large_famid_small_famid_overlap_dict(large_famid_small_famid_overlap_dict, large_famid, small_famid, l_s_overlap, small_famid_famsize)
	
	family_mergings_file.close()
	return([small_fam_large_fam_arr_dict, large_famid_small_famid_overlap_dict, small_famid_famsize_dict])

def fill_large_famid_small_famid_overlap_dict(large_famid_small_famid_overlap_dict, large_famid, small_famid, l_s_overlap, small_famid_famsize):
	if(large_famid_small_famid_overlap_dict.has_key(large_famid)):
		large_famid_small_famid_overlap_dict[large_famid][small_famid] = [l_s_overlap, small_famid_famsize]
	else:
		large_famid_small_famid_overlap_dict[large_famid]={}
		large_famid_small_famid_overlap_dict[large_famid][small_famid]=list()
		large_famid_small_famid_overlap_dict[large_famid][small_famid] = [l_s_overlap, small_famid_famsize]
		

def get_unique_fam_mergings(small_fam_large_fam_arr_dict, small_famid_famsize_dict, temp_dir, fam_fasta_dirName):
	unique_small_fam_large_fam_arr_dict = {}
	for small_famid in small_fam_large_fam_arr_dict:
		if(len(small_fam_large_fam_arr_dict[small_famid])>1):
			build_msa(small_famid, small_fam_large_fam_arr_dict[small_famid], small_famid_famsize_dict, temp_dir, fam_fasta_dirName)
			build_hmm(small_famid, small_fam_large_fam_arr_dict[small_famid], temp_dir, fam_fasta_dirName)
			prepare_hmm_database(small_fam_large_fam_arr_dict[small_famid], temp_dir)
			execute_hhsearch(small_famid, temp_dir)
			best_large_famid = get_best_large_famid_from_hhr(small_famid, temp_dir)
			delete_iteration_files(temp_dir)
			unique_small_fam_large_fam_arr_dict[small_famid] = best_large_famid
			#break
		else:
			large_fam_arr = small_fam_large_fam_arr_dict[small_famid]
			unique_small_fam_large_fam_arr_dict[small_famid] = large_fam_arr[0]

	return(unique_small_fam_large_fam_arr_dict)


def build_msa(small_famid, large_famid_arr, small_famid_famsize_dict, temp_dir, fam_fasta_dirName):
	if(small_famid_famsize_dict[small_famid] < 2):
		copyfile(fam_fasta_dirName+small_famid, temp_dir+small_famid+".msa")
	else:
		build_msa_from_fasta(small_famid, temp_dir, fam_fasta_dirName)

	for large_famid in large_famid_arr:
		build_msa_from_fasta(large_famid, temp_dir, fam_fasta_dirName)


def build_hmm(small_famid, large_famid_arr, temp_dir, fam_fasta_dirName):
	build_hmm_model_from_msa(small_famid, temp_dir)
	for large_famid in large_famid_arr:
		build_hmm_model_from_msa(large_famid, temp_dir)


def build_msa_from_fasta(famid, temp_dir, fam_fasta_dirName):
	msa_outfile = open(temp_dir+famid+".msa", "w")
	run_mafft = subprocess.Popen(["mafft", "--auto", "--amino", fam_fasta_dirName+famid],stdout=msa_outfile, stderr=subprocess.PIPE)
	align = run_mafft.communicate()
	msa_outfile.close()

def build_hmm_model_from_msa(famid, temp_dir):
	run_hmmbuild=subprocess.Popen(["hhmake", "-i", temp_dir+famid+".msa", "-o", temp_dir+famid+".hhm", "-name", famid], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmbuild_results=run_hmmbuild.communicate()


def prepare_hmm_database(large_famid_arr, temp_dir):
	large_famid_hmm_database_file = open(temp_dir+"database.hhm","w")
	for large_famid in large_famid_arr:
		large_famid_hmmfileName = temp_dir+large_famid+".hhm"
		large_famid_hmmfile  = open(large_famid_hmmfileName, "r")
		for line in large_famid_hmmfile:
			line = line.rstrip()
			large_famid_hmm_database_file.write(line+"\n")
		
		large_famid_hmm_database_file.write("\n")
		large_famid_hmmfile.close()

	large_famid_hmm_database_file.close()

def execute_hhsearch(small_famid, temp_dir):
	run_hhsearch  = subprocess.Popen(["hhsearch", "-i", temp_dir+small_famid+".hhm", "-d", temp_dir+"database.hhm", "-o", temp_dir+small_famid+".hhr", "-nocons", "-nopred", "-nodssp"],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	run_hhsearch.communicate()

def delete_iteration_files(directoryName):
	for fl in glob.glob(directoryName+'*'):
    		os.remove(fl)

def get_best_large_famid_from_hhr(small_famid, temp_dir):
	hhr_file = open(temp_dir+small_famid+".hhr", "r")
	
	top_hit_flag=0
	for line in hhr_file:
		line = line.rstrip()
		if(top_hit_flag==1):
			best_large_famid = line[1:]
			break
		if(re.match(r'^No\s+1',line)):
			top_hit_flag=1
	hhr_file.close()
	return(best_large_famid)

def print_unique_family_mergings(large_famid_small_famid_overlap_dict, unique_small_fam_large_fam_arr_dict):
	print "#large_famid small_famid l-s_overlap small_famid_famsize"
	for small_famid in unique_small_fam_large_fam_arr_dict:
		large_famid = unique_small_fam_large_fam_arr_dict[small_famid]
		l_s_overlap, small_famid_famsize = large_famid_small_famid_overlap_dict[large_famid][small_famid]
		print '{0} {1} {2} {3}'.format(large_famid, small_famid, l_s_overlap, small_famid_famsize)
		
##########################################################################################################################
family_mergings_fileName = "/home/aayadav/research/family_quality_orthofinder_eudicot/orthofinder_eudicot_famsize_13_22/orthofinder_eudicot_famsize_13_22.family_mergings"
temp_dir = "/home/aayadav/research/family_quality_orthofinder_eudicot/orthofinder_eudicot_famsize_13_22/temp/"
fam_fasta_dirName = "/data/orthofinder_eudicot/orthofinder_family_fasta/"

small_fam_large_fam_arr_dict, large_famid_small_famid_overlap_dict, small_famid_famsize_dict = get_small_fam_large_fam_arr_dict(family_mergings_fileName)
unique_small_fam_large_fam_arr_dict = get_unique_fam_mergings(small_fam_large_fam_arr_dict, small_famid_famsize_dict, temp_dir, fam_fasta_dirName)
print_unique_family_mergings(large_famid_small_famid_overlap_dict, unique_small_fam_large_fam_arr_dict)

