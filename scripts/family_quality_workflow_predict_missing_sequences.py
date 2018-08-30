#!/usr/bin/python

import re
import sys
import os
import pandas as pd
import numpy as np
import subprocess
import random
import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.metrics import precision_recall_fscore_support, precision_recall_curve, confusion_matrix, cohen_kappa_score, average_precision_score, auc

###############################################################################################################################################################
def execute_phmmer_familyfasta_vs_masterfasta(phmmertlbout_outfileName,family_fasta_fileName, master_fasta_fileName):
	run_phmmer=subprocess.Popen(["phmmer","--tblout",phmmertlbout_outfileName,"--noali",family_fasta_fileName,master_fasta_fileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	run_phmmer_results= run_phmmer.communicate()
################################################################################################################################################################
def get_sequences_from_phmmer_search(family_phmmertblout_fileName, family_fasta_fileName, outseqlist_fileName):
	family_seqid_dict = get_family_seqid_dict(family_fasta_fileName)
	family_phmmertblout_file = open(family_phmmertblout_fileName, "r")
	query_subject_dict={}
	for line in family_phmmertblout_file:
		line=line.rstrip()
		if(re.match(r'^\#',line)):
			continue
		linearr = re.split(r'\s+',line)
		if(query_subject_dict.has_key(linearr[2])):
			query_subject_dict[linearr[2]].append(linearr[0])
		else:
			query_subject_dict[linearr[2]]=list()
			query_subject_dict[linearr[2]].append(linearr[0])

	family_phmmertblout_file.close()
	query_subject_dict, removed_seqdict = remove_worst_nonfamily_sequences(query_subject_dict, family_seqid_dict)
	print_sequence_list(query_subject_dict, removed_seqdict, outseqlist_fileName)

def get_family_seqid_dict(family_fasta_fileName):
	family_seqid_dict={}
	family_fasta_file = open(family_fasta_fileName,"r")
	for line in family_fasta_file:
		line = line.rstrip()
		if(re.match(r'^\>',line)):
			family_seqid_dict[line[1:]]=1

	family_fasta_file.close()
	return(family_seqid_dict)

def remove_worst_nonfamily_sequences(query_subject_dict, family_seqid_dict):
	removed_seqdict = {}
	for query in query_subject_dict:
		seqlist = query_subject_dict[query]
		self_match_flag=0
		for seq in list(reversed(query_subject_dict[query])):
			if not (family_seqid_dict.has_key(seq)):
					seqlist.remove(seq)
					removed_seqdict[seq]=1
			else:
				if(seq==query):
					self_match_flag=1
				break
		if not (self_match_flag==1):
			query_subject_dict[query] = list()
			query_subject_dict[query] = seqlist
	return([query_subject_dict, removed_seqdict])
	


'''
def remove_worst_nonfamily_sequences(query_subject_dict, family_seqid_dict):
	for query in query_subject_dict:
		seqlist = query_subject_dict[query]
		for seq in list(reversed(query_subject_dict[query])):
			if not (family_seqid_dict.has_key(seq)):
					seqlist.remove(seq)
			else:
				break
		query_subject_dict[query] = list()
		query_subject_dict[query] = seqlist
	return(query_subject_dict)



def limit_no_of_nonfamily_sequences(seqid_dict, query_subject_dict):
	no_of_nonfamily_seqs_to_keep=len(query_subject_dict)
	nonfamily_seq_counter=0
	final_seqid_dict={}
	non_family_seqid_dict={}
	for seq in seqid_dict:
		if(query_subject_dict.has_key(seq)):
			final_seqid_dict[seq]=1
		else:
			non_family_seqid_dict[seq]=1
	
	selected_non_family_seqids=random.sample(non_family_seqid_dict, no_of_nonfamily_seqs_to_keep)
	for non_family_id in selected_non_family_seqids:
		final_seqid_dict[non_family_id]=1
	return(final_seqid_dict)
'''
def print_sequence_list(query_subject_dict, removed_seqdict, outseqlist_fileName):
	seqid_dict={}
	for query in query_subject_dict:
		for seq in query_subject_dict[query]:
			seqid_dict[seq]=1
	
	if(len(seqid_dict) == len(query_subject_dict)):
		removed_seq_arr = removed_seqdict.keys()
		if(len(removed_seq_arr)>=len(query_subject_dict)):
			for i in range(0, len(query_subject_dict)):
				seqid_dict[removed_seq_arr[i]]=1
		else:
			for i in range(0, len(removed_seq_arr)):
				seqid_dict[removed_seq_arr[i]]=1
	#this else is to limit the number of non-family sequences in-case huge number of non-family seqs are found
	#else:
	#	seqid_dict=limit_no_of_nonfamily_sequences(seqid_dict, query_subject_dict)
	
	outseqlist_file = open(outseqlist_fileName, "w")
	for seq in seqid_dict:
		outseqlist_file.write(seq+"\n")

	outseqlist_file.close()
##################################################################################################################################################################
def get_fasta_from_sequence_list(seqlist_fileName, master_fasta_fileName, outfasta_fileName):
	master_fasta_file = open(master_fasta_fileName,"r")
	current_seqid=""
	seqid_sequence_dict={}
	for line in master_fasta_file:
		line=line.rstrip()
		if(re.match(r'^\>',line)):
			current_seqid = line[1:]
			seqid_sequence_dict[current_seqid]=""
		else:
			seqid_sequence_dict[current_seqid] = seqid_sequence_dict[current_seqid]+line

	
	master_fasta_file.close()
	
	outfasta_file = open(outfasta_fileName,"w")
	seqlist_file = open(seqlist_fileName, "r")
	for line in seqlist_file:
		line = line.rstrip()
		outfasta_file.write('>'+line+"\n")
		outfasta_file.write(seqid_sequence_dict[line]+"\n")
	seqlist_file.close()
	outfasta_file.close()
###################################################################################################################################################################
def get_pairs_file_from_seqlist(seqlist_fileName, family_fasta_fileName, pairs_outfileName):
	family_fasta_file = open(family_fasta_fileName,"r")
	seqid_dict={}
	for line in family_fasta_file:
		line = line.rstrip()
		if(re.match(r'^\>',line)):
			seqid_dict[line[1:]]=1
	family_fasta_file.close()
	
	seqlist_file = open(seqlist_fileName, "r")
	seqlist_arr=list()
	for line in seqlist_file:
		line = line.rstrip()
		seqlist_arr.append(line)
	
	seqlist_file.close()

	pairs_outfile = open(pairs_outfileName,"w")

	for i in range(0, len(seqlist_arr)):
		for j in range(i+1, len(seqlist_arr)):
			if(seqid_dict.has_key(seqlist_arr[i]) and seqid_dict.has_key(seqlist_arr[j])):
				#print '10000\t{0}\t{1}'.format(seqlist_arr[i], seqlist_arr[j])
				pairs_outfile.write(str(10000)+"\t"+seqlist_arr[i]+"\t"+seqlist_arr[j]+"\n")

			#else:
			#	#print '99999\t{0}\t{1}'.format(seqlist_arr[i], seqlist_arr[j])
			#	pairs_outfile.write(str(99999)+"\t"+seqlist_arr[i]+"\t"+seqlist_arr[j]+"\n")

			elif (seqid_dict.has_key(seqlist_arr[i]) or seqid_dict.has_key(seqlist_arr[j])):
				pairs_outfile.write(str(99999)+"\t"+seqlist_arr[i]+"\t"+seqlist_arr[j]+"\n")
	
	pairs_outfile.close()
	return(seqid_dict)
###################################################################################################################################################################
def get_seqid_sequence_dict_from_fasta(trainingfasta_fileName):
	seqid_sequence_dict = {}
	trainingfasta_file=open(trainingfasta_fileName,'r')
	for line in trainingfasta_file:
		line=line.rstrip()
		if(re.match(r'^\>',line)):
			current_sequenceid = line[1:]
			seqid_sequence_dict[current_sequenceid]=''
		else:
			seqid_sequence_dict[current_sequenceid] = seqid_sequence_dict[current_sequenceid]+line

	trainingfasta_file.close()

	return(seqid_sequence_dict)

def get_training_pairs_dataframe(trainingpairs_fileName):
	trainingpairs_dataframe=pd.DataFrame()
	trainingpairs_dataframe=pd.read_table(trainingpairs_fileName, sep="\s+", header=None)

	return(trainingpairs_dataframe)

def attach_class_label_column(trainingpairs_dataframe, targetfamid):
	trainingpairs_dataframe['label']=""
	trainingpairs_dataframe.loc[trainingpairs_dataframe[0]==targetfamid,'label']='pos'
	trainingpairs_dataframe.loc[trainingpairs_dataframe[0]!=targetfamid,'label']='neg'

def hmm_model_training_evaluation_stratTestTrainSplit(trainingpairs_dataframe, targetfamid, famid_name):
	pos_pairs_table = trainingpairs_dataframe.loc[trainingpairs_dataframe['label']=='pos']
	neg_pairs_table = trainingpairs_dataframe.loc[trainingpairs_dataframe['label']=='neg']

	prediction_dataframe=pd.DataFrame()
	prediction_table = pd.DataFrame()
	test_train_split_iterations = 10
	for iteration in range(0,test_train_split_iterations):
		print 'iteration-no......{0}'.format(iteration)

		train_index = np.random.rand(len(pos_pairs_table)) < 0.8
		
		train_split = pos_pairs_table[train_index]
		train_split = train_split.reset_index(drop=True)
		test_split = pd.concat([pos_pairs_table[~train_index],neg_pairs_table],ignore_index=True)
		
		train_hmm_model(train_split)
		test_hmm_model(test_split)
		
		prediction_table = get_test_pair_prediction_dataframe(test_split) 
		prediction_dataframe = pd.concat([prediction_dataframe, prediction_table], ignore_index=True)

		delete_iteration_files(outPath+'/temp/')
		#break

	prediction_performance_stats = get_prediction_performance_stats(prediction_dataframe, famid_name)
	return(prediction_performance_stats)

def train_hmm_model(train_split):
	write_sequences_to_fasta_file(train_split, "training.fasta")
	build_msa_from_fasta("training.fasta")
	build_hmm_model_from_msa("training.fasta.msa")


def write_sequences_to_fasta_file(pairs_dataframe, fileName):

	fasta_outfile = open(outPath+'/temp/'+fileName,"w")
	if(fileName=='training.fasta'):
		for index, testrow in pairs_dataframe.iterrows():
			seq1 = testrow.ix[1] 
			seq2 = testrow.ix[2]
			fasta_outfile.write(">"+seq1+"#"+seq2+"\n")
			fasta_outfile.write(seqpair_consensus_dict[seq1][seq2]+"\n")
	else:
		seq_list = np.array(list(pairs_dataframe.ix[:,1])+list(pairs_dataframe.ix[:,2]))
        	seq_list = np.unique(seq_list)
        	for seq in seq_list:
        	        fasta_outfile.write(">"+seq+"\n")
        	        fasta_outfile.write(seqid_sequence_dict[seq]+"\n")

	fasta_outfile.close()

def build_msa_from_fasta(fastafileName):
	run_muscle=subprocess.Popen(["/home/aayadav/Downloads/muscle3.8.31_i86linux64","-in",outPath+"/temp/"+fastafileName,"-out",outPath+"/temp/"+fastafileName+".msa"],stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	align = run_muscle.communicate()

def build_hmm_model_from_msa(msafileName):
	run_hmmbuild=subprocess.Popen(["hmmbuild","--amino", outPath+"/temp/"+msafileName+".hmm", outPath+"/temp/"+msafileName], stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmbuild_results=run_hmmbuild.communicate()

def test_hmm_model(test_split):
	write_sequences_to_fasta_file(test_split, "testing.fasta")
	execute_hmm_search("training.fasta.msa.hmm", "testing.fasta")
	
def execute_hmm_search(hmm_model_fileName, testingfasta_fileName):
	run_hmmsearch=subprocess.Popen(["hmmsearch","--tblout",outPath+"/temp/"+testingfasta_fileName+".hmmsearch","--noali", outPath+"/temp/"+hmm_model_fileName, outPath+"/temp/"+testingfasta_fileName],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmsearch_results= run_hmmsearch.communicate()

def get_test_pair_prediction_dataframe(test_split):
	test_split.columns = ['clusterid', 'seq1', 'seq2', 'label']
	seqid_evalue_domcount_dict={}

	seq1_values=list()
	seq2_values=list()

	seq1_domcount=list()
	seq2_domcount=list()

	seq1_score=list()
	seq2_score=list()

	seqid_evalue_domcount_dict = get_seqid_evalue_domcount_dict("testing.fasta.hmmsearch")
	
	for index, testrow in test_split.iterrows():
		if(seqid_evalue_domcount_dict.has_key(testrow['seq1'])):
			seq1_values.append(seqid_evalue_domcount_dict[testrow['seq1']]['evalue'])
			seq1_domcount.append(seqid_evalue_domcount_dict[testrow['seq1']]['domcount'])
			seq1_score.append(seqid_evalue_domcount_dict[testrow['seq1']]['score'])
		else:
			seq1_values.append(100)
			seq1_domcount.append(0)
			seq1_score.append(0)

		if(seqid_evalue_domcount_dict.has_key(testrow['seq2'])):
			seq2_values.append(seqid_evalue_domcount_dict[testrow['seq2']]['evalue'])
			seq2_domcount.append(seqid_evalue_domcount_dict[testrow['seq2']]['domcount'])
			seq2_score.append(seqid_evalue_domcount_dict[testrow['seq2']]['score'])
		else:
			seq2_values.append(100)
			seq2_domcount.append(0)
			seq2_score.append(0)

	seq1_values = pd.Series(seq1_values)
	seq2_values = pd.Series(seq2_values)

	seq1_domcount = pd.Series(seq1_domcount)
	seq2_domcount = pd.Series(seq2_domcount)

	seq1_score = pd.Series(seq1_score)
	seq2_score = pd.Series(seq2_score)


	test_split = test_split.assign(seq1_evalue = seq1_values.values)
	test_split = test_split.assign(seq2_evalue = seq2_values.values)

	test_split = test_split.assign(seq1_domcount = seq1_domcount.values)
	test_split = test_split.assign(seq2_domcount = seq2_domcount.values)

	test_split = test_split.assign(seq1_score = seq1_score.values)
	test_split = test_split.assign(seq2_score = seq2_score.values)

	return(test_split)

def get_seqid_evalue_domcount_dict(hmmsearch_fileName):
	seqid_evalue_domcount_dict={}
	hmmsearch_file = open(outPath+"/temp/"+hmmsearch_fileName,"r")
	for line in hmmsearch_file:
		line=line.rstrip()
		if not (re.match(r'^\#',line)):
			linearr=re.split(r'\s+',line)
			seqid_evalue_domcount_dict[linearr[0]]={}
			seqid_evalue_domcount_dict[linearr[0]]['evalue']=float(linearr[4])
			seqid_evalue_domcount_dict[linearr[0]]['domcount']=float(linearr[15])
			seqid_evalue_domcount_dict[linearr[0]]['score']=float(linearr[5])
	
	hmmsearch_file.close()

	return(seqid_evalue_domcount_dict)

def delete_iteration_files(directoryName):
	for fl in glob.glob(directoryName+'*'):
    		os.remove(fl)


def get_prediction_performance_stats(prediction_dataframe, famid_name):
	precision_vals, recall_vals, best_cutoff, lowest_cutoff = get_precision_recall_vals_for_scores(prediction_dataframe, famid_name)
	pos_pr_auc = get_pr_auc(precision_vals, recall_vals)

	plot_pr_curve(precision_vals, recall_vals, famid_name)

	prediction_performance_stats = get_predicition_stats_for_specified_cutoff_for_scores(prediction_dataframe, best_cutoff)

	prediction_performance_stats.loc[:,'pos_pr_auc'] = list(np.repeat(pos_pr_auc,prediction_performance_stats.shape[0]))
	prediction_performance_stats.loc[:,'best_score_cutoff'] = list(np.repeat(best_cutoff, prediction_performance_stats.shape[0]))
	prediction_performance_stats.loc[:,'lowest_score_cutoff'] = list(np.repeat(lowest_cutoff, prediction_performance_stats.shape[0]))
	return(prediction_performance_stats)

def get_score_cutoffs_for_pr_curve(pos_prediction_dataframe, famid_name):
	largest_cutoff=np.sort(np.array(list(pos_prediction_dataframe['seq1_score'])+list(pos_prediction_dataframe['seq2_score'])))
        #largest_cutoff=largest_cutoff[-1]
        score_list=list()
        for index, testrow in pos_prediction_dataframe.iterrows():
                if(testrow['seq1_score']<testrow['seq2_score']):
                        score_list.append(testrow['seq1_score'])
                else:
                        score_list.append(testrow['seq2_score'])

        #score_list.append(largest_cutoff)
        score_list = np.array(score_list)
        score_list = np.unique(score_list)
        score_list = np.sort(score_list)

        #score_list = score_list-1

	score_cutoff_outfile = open(outPath+"/"+famid_name+"/"+famid_name+".score_cutoffs","w")
	score_cutoff_outfile.write(np.array_str(score_list))
	score_cutoff_outfile.close()

        return(score_list)

	

def get_precision_recall_vals_for_scores(prediction_dataframe, famid_name):
        precision_vals=list()
        recall_vals=list()
        fscore_vals=list()
        prediction_dataframe['temp_pred_label']=""
        cutoffs_arr = get_score_cutoffs_for_pr_curve(prediction_dataframe.loc[prediction_dataframe['label']=='pos'], famid_name)
        best_fscore=0
        best_cutoff=1
        for threshold in cutoffs_arr:
                prediction_dataframe['temp_pred_label']=""
                prediction_dataframe.loc[(prediction_dataframe['seq1_score']>=threshold) & (prediction_dataframe['seq2_score']>=threshold),'temp_pred_label']='pos'
                prediction_dataframe.loc[(prediction_dataframe['seq1_score']<threshold) | (prediction_dataframe['seq2_score']<threshold),'temp_pred_label']='neg'
                precision, recall, fscore, support = precision_recall_fscore_support(prediction_dataframe['label'],prediction_dataframe['temp_pred_label'], labels=['pos','neg'])
                if(fscore[0]>best_fscore):
               		best_cutoff=threshold
                	best_fscore=fscore[0]

                precision_vals.append(precision[0])
                recall_vals.append(recall[0])
                fscore_vals.append(fscore[0])



        precision_vals.append(1)
        recall_vals.append(0)

        precision_vals = np.array(precision_vals)
        recall_vals = np.array(recall_vals)
	
	lowest_cutoff = cutoffs_arr[0]
		
        #print precision_vals
        #print recall_vals
        #print fscore_vals

        return([precision_vals, recall_vals, best_cutoff, lowest_cutoff])


def get_pr_auc(precision_vals, recall_vals):
	pos_pr_auc = auc(recall_vals, precision_vals)
	return(pos_pr_auc)
	
#this function plots pr-curve for the pos class for the combined test dataset using the true labels and then predictions for the pos class.
def plot_pr_curve(precision_vals, recall_vals, famid_name):
	plt.clf()
	plt.plot(recall_vals, precision_vals, lw=1, color='navy',label='Precision-Recall curve')
	plt.xlabel('Recall')
	plt.ylabel('Precision')
	plt.ylim([0.0, 1.05])
	plt.xlim([0.0, 1.05])
	plt.savefig(outPath+"/"+famid_name+"/"+famid_name+"pr_curve",format='png')
	#plt.show()

def get_predicition_stats_for_specified_cutoff_for_scores(prediction_dataframe,best_cutoff):
        prediction_dataframe['temp_pred_label']=""
        prediction_dataframe.loc[(prediction_dataframe['seq1_score']>=best_cutoff) & (prediction_dataframe['seq2_score']>=best_cutoff),'temp_pred_label']='pos'
        prediction_dataframe.loc[(prediction_dataframe['seq1_score']<best_cutoff) | (prediction_dataframe['seq2_score']<best_cutoff),'temp_pred_label']='neg'
        precision, recall, fscore, support = precision_recall_fscore_support(prediction_dataframe['label'],prediction_dataframe['temp_pred_label'], labels=['pos','neg'])

        prediction_performance_stats = list(precision) + list(recall) + list(fscore) + list(support)


        prediction_performance_stats = pd.DataFrame([prediction_performance_stats])
        prediction_performance_stats.columns = ["pos-prec","neg-prec","pos-recl","neg-recl","pos-fsco","neg-fsco","pos-supp","neg-supp"]

        return(prediction_performance_stats)



def print_performance_stat_dataframe_to_file(prediction_performance_stats, famid_name):
	prediction_performance_stats.to_csv(outPath+"/"+famid_name+"/"+famid_name+".performance_stats",index=None,sep=' ',mode='w')


def emit_consensus_sequence(hmm_model_file, seq1, seq2):
	run_hmmemit=subprocess.Popen(["hmmemit", "-c", outPath+"temp/"+hmm_model_file ],stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	hmmemit_results= run_hmmemit.communicate()
	consensus_seq = hmmemit_results[0]
	consensus_seq_arr=re.split(r'\n',consensus_seq)
	consensus_seq_arr[0] = ">"+seq1+"#"+seq2
	consensus_seq_arr = [ consensus_seq_arr[0], ''.join(consensus_seq_arr[1:])]
	return(consensus_seq_arr)
	

def write_pair_consensus_fasta(trainingpairs_fileName, pospair_outfasta_fileName):
	trainingpairs_dataframe=pd.DataFrame()
	trainingpairs_dataframe=pd.read_table(trainingpairs_fileName, sep="\s+", header=None)
	pospair_outfasta_file=open(pospair_outfasta_fileName,"w")
	for index, testrow in trainingpairs_dataframe.iterrows():
		if(testrow.ix[0]!=99999):
			write_sequences_to_fasta_file(trainingpairs_dataframe.iloc[[index]], "pair.fasta")
			build_msa_from_fasta("pair.fasta")
			build_hmm_model_from_msa("pair.fasta.msa")
			consensus_seq_arr = emit_consensus_sequence("pair.fasta.msa.hmm", testrow.ix[1], testrow.ix[2])
			pospair_outfasta_file.write(consensus_seq_arr[0]+"\n"+consensus_seq_arr[1]+"\n")
			delete_iteration_files(outPath+'/temp/')
			#break
	pospair_outfasta_file.close()

def get_sequence_pair_consensus_dict(pospair_fasta_fileName):
	pospair_fasta_file = open(pospair_fasta_fileName, "r")
	seqpair_consensus_dict = {}
	for line in pospair_fasta_file:
		line=line.rstrip()
		if(re.match(r'^>',line)):
			seqheader = line[1:]
			seq1, seq2 = re.split(r'\#',seqheader)
		else:
			if(seqpair_consensus_dict.has_key(seq1)):
				seqpair_consensus_dict[seq1][seq2]=line
			else:
				seqpair_consensus_dict[seq1]={}
				seqpair_consensus_dict[seq1][seq2]=line
	
			if(seqpair_consensus_dict.has_key(seq2)):
				seqpair_consensus_dict[seq2][seq1]=line
			else:
				seqpair_consensus_dict[seq2]={}
				seqpair_consensus_dict[seq2][seq1]=line
	pospair_fasta_file.close()
	return(seqpair_consensus_dict)
####################################################################################################################################################################
def predict_missing_sequences(trainingpairs_dataframe, famid_name, prediction_performance_stats):
#def predict_missing_sequences(famid_name):
	
	#global seqid_sequence_dict
	#global seqpair_consensus_dict
	
	#trainingpairs_dataframe = get_training_pairs_dataframe(outPath+"/"+famid_name+"/"+famid_name+".pairs")
	#seqid_sequence_dict = get_seqid_sequence_dict_from_fasta(outPath+"/"+famid_name+"/"+famid_name+"_closest_sequences.fa")
	#seqpair_consensus_dict = get_sequence_pair_consensus_dict(outPath+"/"+famid_name+"/"+famid_name+".pospairs-consensus.fasta")
	
	targetfamid = 10000
	attach_class_label_column(trainingpairs_dataframe, targetfamid)
	
	pos_pairs_table = trainingpairs_dataframe.loc[trainingpairs_dataframe['label']=='pos']
	neg_pairs_table = trainingpairs_dataframe.loc[trainingpairs_dataframe['label']=='neg']
	
	train_hmm_model(pos_pairs_table)
	test_hmm_model(neg_pairs_table)

	neg_prediction_table = get_test_pair_prediction_dataframe(neg_pairs_table)
	
	lowest_score_cutoff = read_lowest_score_from_family_model_stats_file(outPath+"/"+famid_name+"/"+famid_name+".performance_stats")

	get_missing_sequences(neg_prediction_table, lowest_score_cutoff, famid_name)

	delete_iteration_files(outPath+'/temp/')
			

def get_missing_sequences(neg_prediction_table, lowest_score_cutoff, famid_name):
	accepted_sequences_arr=list()

	for index, testrow in neg_prediction_table.iterrows():
		if((testrow['seq1_score']>lowest_score_cutoff) and (testrow['seq2_score']>lowest_score_cutoff)):
			accepted_sequences_arr.append(testrow['seq1'])
			accepted_sequences_arr.append(testrow['seq2'])
	
	print_missing_sequences(accepted_sequences_arr, famid_name)
	

def print_missing_sequences(accepted_sequences_arr, famid_name):
	accepted_sequences_arr=list(set(accepted_sequences_arr))
	
	missing_sequences_outfile=open(outPath+"/"+famid_name+"/"+famid_name+".missing_sequences","w")

	for accepted_seq in accepted_sequences_arr:
		if not (seq_dict_from_fam_fasta.has_key(accepted_seq)):
			missing_sequences_outfile.write(accepted_seq+"\n")

	missing_sequences_outfile.close()

def read_lowest_score_from_family_model_stats_file(family_model_stats_fileName):
	family_model_stats_dataframe=pd.DataFrame()
	family_model_stats_dataframe=pd.read_table(family_model_stats_fileName, sep="\s+")
	lowest_score_cutoff = family_model_stats_dataframe.loc[[0],['lowest_score_cutoff']]
	lowest_score_cutoff = float (lowest_score_cutoff.values[0])
	return(lowest_score_cutoff)

####################################################################################################################################################################
def family_model_training_and_evaluation(famid_name):
	global seqid_sequence_dict
	global seqpair_consensus_dict

	trainingpairs_dataframe = get_training_pairs_dataframe(outPath+"/"+famid_name+"/"+famid_name+".pairs")
	seqid_sequence_dict = get_seqid_sequence_dict_from_fasta(outPath+"/"+famid_name+"/"+famid_name+"_closest_sequences.fa")
	
	write_pair_consensus_fasta(outPath+"/"+famid_name+"/"+famid_name+".pairs", outPath+"/"+famid_name+"/"+famid_name+".pospairs-consensus.fasta")
	seqpair_consensus_dict = get_sequence_pair_consensus_dict(outPath+"/"+famid_name+"/"+famid_name+".pospairs-consensus.fasta")
	
	
	
	targetfamid = 10000
	print 'processing......{0}'.format(famid_name)
	attach_class_label_column(trainingpairs_dataframe, targetfamid)
	prediction_performance_stats = hmm_model_training_evaluation_stratTestTrainSplit(trainingpairs_dataframe, str(targetfamid), famid_name)
	print_performance_stat_dataframe_to_file(prediction_performance_stats, famid_name)
	
	predict_missing_sequences(trainingpairs_dataframe, famid_name, prediction_performance_stats)
	



###################################################################################################################################################################
###################################################################################################################################################################
global outPath
global seq_dict_from_fam_fasta
outPath = "/home/aayadav/research/family_quality_ygob/delete_families/"
famid_name = str(sys.argv[1])
###################################################################################################################################################################

#execute_phmmer_familyfasta_vs_masterfasta("L.17R79/L.17R79.phmmertlbout", "/data/legume_genefams3/32_family_fasta/L.17R79", "/data/legume_genefams3/legume_genefams3.fa")
execute_phmmer_familyfasta_vs_masterfasta(famid_name+"/"+famid_name+".phmmertlbout", "/data/ygob/family_fasta_insert20/"+famid_name, "/data/ygob/ygob.fasta")
###################################################################################################################################################################

#get_sequences_from_phmmer_search("L.17R79/L.17R79.phmmertlbout","/data/legume_genefams3/32_family_fasta/L.17R79")
get_sequences_from_phmmer_search(famid_name+"/"+famid_name+".phmmertlbout","/data/ygob/family_fasta_insert20/"+famid_name, famid_name+"/"+famid_name+".seqlist")
###################################################################################################################################################################

#get_fasta_from_sequence_list("family_phmmerout/L.17R79.seqlist", "/data/legume_genefams3/legume_genefams3.fa")
get_fasta_from_sequence_list(famid_name+"/"+famid_name+".seqlist", "/data/ygob/ygob.fasta", famid_name+"/"+famid_name+"_closest_sequences.fa")
###################################################################################################################################################################

seq_dict_from_fam_fasta = get_pairs_file_from_seqlist(famid_name+"/"+famid_name+".seqlist", "/data/ygob/family_fasta_insert20/"+famid_name, famid_name+"/"+famid_name+".pairs")
###################################################################################################################################################################

family_model_training_and_evaluation(famid_name)
###################################################################################################################################################################
#seq_dict_from_fam_fasta = get_family_seqid_dict("/data/ygob/family_fasta_delete/"+famid_name)
#predict_missing_sequences(famid_name)
###################################################################################################################################################################