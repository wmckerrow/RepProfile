"""

Copyright Wilson McKerrow, 2017

This file is part of RepProfile.

RepProfile is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

RepProfile is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with RepProfile.  If not, see <http://www.gnu.org/licenses/>.


"""

"""

This is the main RepProfile code. It takes alignments preprocessed from a bamfile by
parse_PEbam_for_RepProfile.py and runs EM to estimate hyper-editing. Results are saved 
in several pickle files. Text output can be generated from these by utilities/report.py.

"""

import numpy
from Bio import SeqIO
import argparse
import pickle
import math
import copy
import scipy
import scipy.stats
import scipy.special
import datetime
import sys
from itertools import chain
from sets import Set
from multiprocessing import Pool

# ERROR AND INDEL PARAMETERS
P_ERROR = 0.001 # Probability of a given sequencing error (e.g. A->C)
P_OPEN = 0.001 # Probability of opening a new indel
P_EXTEND = 0.1 # Probability of extending an indel

# Values need to define P(H,T,G) = P(H)P(T|H)P(G|T)
class dirmixprior(object):
	def __init__(self, pseudocounts_plus,pseudocounts_minus,strands_same,hyper_dict,pos_probs):
		self.pseudocounts_plus = pseudocounts_plus # Defines P(G_i|T_i) = Dirichlet(alpha_i) where alpha_i = pseudocounts_i+1. pseudocounts_plus[seq][type][i,j] is the pseudocounts for position i in repeat seq for nucleotide j given T_i = type. pseudocounts_plus['FB4_DM_4']['edit_f'][2000,0] = 1 would mean that at position 2000 of repeat FB4_DM_4 in the position state edit_f (forward strand edit) there is 1 adenosine pseudocount.
		self.pseudocounts_minus = pseudocounts_minus # Same as pseudocounts_plus but for reverse strand
		self.strands_same = strands_same # Whether the profile is allowed to differ between forward and reverse strand. Currently it is not recommended to allow this. strands_same[type] = TRUE means that when T_i=type the forward and reverse strands must have the same profile.
		self.hyper_dict = hyper_dict # Defines P(H): P(H_i = hyper) = hyper_dict[hyper]. hyper_dict['hyper_f'] = 0.01 means that each repeat has a 0.01 chance of being in the 'hyper_f' state.
		self.pos_probs = pos_probs # Defines P(T_i|H_{k(i)}). P(T_i=pos_type|H_{k(i)}=hyper_type,ref(i)=base) = pos_probs[pos_type][hyper_type][base]. pos_probs['edit_f']['hyper_f'][0] = 0.5, means that each A position in a hyper_f (forward strand hyper editing) repeat has a 0.5 chance of being edit_f (forward strand edited).
		self.ln_dir_f = dict() # Dirichlet distribution normalizing constants
		self.ln_dir_r = dict() # Dirichlet distribution normalizing constants
		seq = pseudocounts_plus.keys()[1]
		for type in pseudocounts_plus[seq]:
			self.ln_dir_f[type] = sum(scipy.special.gammaln(pseudocounts_plus[seq][type][1,]+1.0)) - scipy.special.gammaln(sum(pseudocounts_plus[seq][type][1,]+1.0))
			self.ln_dir_r[type] = sum(scipy.special.gammaln(pseudocounts_minus[seq][type][1,]+1.0)) - scipy.special.gammaln(sum(pseudocounts_minus[seq][type][1,]+1.0))
	def get_for_seq(self,seq):
		return dirmixprior_seq(self.pseudocounts_plus[seq],self.pseudocounts_minus[seq],self.strands_same,self.hyper_dict,self.pos_probs,self.ln_dir_f,self.ln_dir_r)
		
# dirmixprior values for a single repeat.
class dirmixprior_seq(object):
	def __init__(self, pseudocounts_plus,pseudocounts_minus,strands_same,hyper_dict,pos_probs,ln_dir_f,ln_dir_r):
		self.pseudocounts_plus = pseudocounts_plus
		self.pseudocounts_minus = pseudocounts_minus
		self.strands_same = strands_same
		self.hyper_dict = hyper_dict
		self.pos_probs = pos_probs
		self.ln_dir_f = ln_dir_f
		self.ln_dir_r = ln_dir_r

# Sequence (seq) and candidate alignments (aln) for a single read end.
class readclass(object):
	def __init__(self, read,rname_list):
		self.seq = seq2array(read.seq)
		self.aln = [alnclass(rname_list[read.rname],read.pos,read.cigarstring,not read.is_reverse)]
	def addaln(self,XA):
		self.aln += [alnclass(XA[0],int(XA[1][1:])-1,XA[2],XA[1][0]=='+')]
		
# A candidate alignment
class alnclass(object):
	def __init__(self,chrom,ref_start,cigarstring,forward):
		self.chrom = chrom # repeat that the alignment is to
		self.forward = forward # TRUE if the alignment is in the forward direction
		self.read_ranges,self.ref_ranges,self.insert_starts,self.insert_extends,self.del_starts,self.del_extends = cigarstring_2_pairs(cigarstring,ref_start) # See main/parse_PEbam_for_RepProfile.py
	def aln_read(self): # aligned positions int the read
		result = []
		for x in self.read_ranges:
			result += range(x[0],x[1])
		return result
		#return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.read_ranges])))
	def aln_ref(self): # aligned positions in the reference. aln_ref()[i] is aligned to aln_read()[i]
		result = []
		for x in self.ref_ranges:
			result += range(x[0],x[1])
		return result
		#return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.ref_ranges])))

# Sequence (seq1,seq2) and alignments (aln1,aln2) for a paired end read		
class PEreadclass(object):
	def __init__(self, read1,read2,aln1,aln2):
		self.seq1 = read1.seq # Sequence of read end 1
		self.seq2 = read2.seq # Sequence of read end 2
		self.alns1 = [aln1] # Candidate alignments for read end 1
		self.alns2 = [aln1] # Candidate alignments for read end 2
		self.aln = [(0,0)] # Alignment pairs that are concordant
	def addaln1(self,aln1):
		self.alns1.append( aln1 )
	def addaln2(self,aln2):
		self.alns2.append( aln2 )
	def addaln(self,i,j):
		self.aln.append( (i,j) )

# Convert sequences from string to numpy array.
def seq2array(seq):
	nuc2num = {'A':0,'C':1,'G':2,'T':3, 'N':4}
	return numpy.array([nuc2num[c] for c in seq],dtype=numpy.int8)

# Convert sequences from numpy array to string.
def array2seq(seq_array):
	seq = ''
	for letter in seq_array:
		seq += 'ACGTN'[letter]
	return seq

# Store sequence in a fasta file a dictionary of numpy arrays
def fasta_as_array(fasta_file):
	fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file,'r'), "fasta"))
	for seq in fasta_dict:
		fasta_dict[seq] = seq2array(fasta_dict[seq].upper())
	return fasta_dict

# Get reverse complement of a numpy array sequence	
def rev_comp(seq_array):
	revcompseq = 3 - seq_array[::-1]
	revcompseq[revcompseq<0] = 4
	return revcompseq

# Build an initial profile guess for EM based on the reference genome. Each row is the provided probs vector translated mod 4 so that the first element of the probs vector is in the reference base position.	
def create_initial_guess(ref,probs):
	genome_profile = dict()
	for seq_name in ref:
		genome_profile[seq_name] = numpy.zeros((len(ref[seq_name]),4))
		for i in range(len(ref[seq_name])):
			if ref[seq_name][i] != 4:
				genome_profile[seq_name][i,] = probs[(numpy.arange(4)-ref[seq_name][i])%4]
			else:
				genome_profile[seq_name][i,] = numpy.array([0.25,0.25,0.25,0.25])
	return genome_profile
	
# Similar to the above function, but with a pseudocount vector instead of nucleotide probabilities.
def create_pseudocounts(refseq,pseudoarray):
	pseudo_counts = numpy.zeros((len(refseq),4))
	for i in range(len(refseq)):
		if refseq[i] != 4:
			pseudo_counts[i,] = pseudoarray[(numpy.arange(4)-refseq[i])%4]
		else:
			pseudo_counts[i,] = numpy.zeros(4)
	return pseudo_counts

# Perform the E step of EM, calculating E_A[U] and E_A[V].
def EStep(input):
	
	#print 'starting E step part'
	
	alignment_pickle,reference,genome_profile_f,genome_profile_r,f_prob,r_prob,flanking = input
	alignments = pickle.load(open(alignment_pickle,'rb')).values()
	ll = 0
	expU_f = dict()
	expU_r = dict()
	for seq in reference:
		expU_f[seq] = numpy.zeros((len(reference[seq])-2*flanking,4))
		expU_r[seq] = numpy.zeros((len(reference[seq])-2*flanking,4))
	alned_reads = 0
	
	for read in alignments:
		if len(read.aln) == 0:
			continue
		unnorm_p1 = numpy.zeros(len(read.alns1))
		for i in range(len(read.alns1)):
			aln1 = read.alns1[i]
			if max(aln1.aln_ref()) >= len(reference[aln1.chrom]):
				continue
			# unnorm_p1[i] = P(R_i,A_i|G,X) \propto \prod(G_{A_i}(R_i)) * P(A_i|X) * indel_probs
			if aln1.forward:
				unnorm_p1[i] = numpy.prod(genome_profile_r[aln1.chrom][aln1.aln_ref(),read.seq1[aln1.aln_read()]])*r_prob[aln1.chrom]*P_OPEN**(aln1.insert_starts+aln1.del_starts)*P_EXTEND**(aln1.insert_extends+aln1.del_extends)*0.25**(aln1.insert_starts+aln1.insert_extends)
			else:
				unnorm_p1[i] = numpy.prod(genome_profile_f[aln1.chrom][aln1.aln_ref(),rev_comp(read.seq1)[aln1.aln_read()]])*f_prob[aln1.chrom]*P_OPEN**(aln1.insert_starts+aln1.del_starts)*P_EXTEND**(aln1.insert_extends+aln1.del_extends)*0.25**(aln1.insert_starts+aln1.insert_extends)
		
		unnorm_p2 = numpy.zeros(len(read.alns2))
		for j in range(len(read.alns2)):
			aln2 = read.alns2[j]
			if max(aln2.aln_ref()) >= len(reference[aln2.chrom]):
				continue
			# unnorm_p2[i] = P(R_i,A_i|G,X) \propto \prod(G_{A_i}(R_i)) * P(A_i|X) * indel_probs
			if aln2.forward:
				unnorm_p2[j] = numpy.prod(genome_profile_f[aln2.chrom][aln2.aln_ref(),read.seq2[aln2.aln_read()]])*f_prob[aln2.chrom]*P_OPEN**(aln2.insert_starts+aln2.del_starts)*P_EXTEND**(aln2.insert_extends+aln2.del_extends)*0.25**(aln2.insert_starts+aln2.insert_extends)
			else:
				unnorm_p2[j] = numpy.prod(genome_profile_r[aln2.chrom][aln2.aln_ref(),rev_comp(read.seq2)[aln2.aln_read()]])*r_prob[aln2.chrom]*P_OPEN**(aln2.insert_starts+aln2.del_starts)*P_EXTEND**(aln2.insert_extends+aln2.del_extends)*0.25**(aln2.insert_starts+aln2.insert_extends)
			
		unnorm_p = numpy.zeros(len(read.aln))
# 			print read.aln
# 			print len(read.alns1)
# 			print len(read.alns2)
# 			print
		for i in range(len(read.aln)):
			unnorm_p[i] = unnorm_p1[read.aln[i][0]]*unnorm_p2[read.aln[i][1]] # Prob of paired alignment is product of single end probs
		if sum(unnorm_p)>0:
			# marginal_p = P(A_i|R_i,G,X) = P(R_i,A_i|G,X)/sum_{A_i}P(R_i,A_i|G,X), a k-long vector, where k is the number of candidate alignments
			marginal_p = unnorm_p/sum(unnorm_p)
			ll += math.log(sum(unnorm_p))
			alned_reads += 1
		else:
			continue
		
		# expU = E_A[U] = \sum_i U(A_i,R_i) * P(A_i|R_i,G,X)
		for i in range(len(read.aln)):
			aln1 = read.alns1[read.aln[i][0]]
			aln2 = read.alns2[read.aln[i][1]]
			if max(aln1.aln_ref()) >= len(reference[aln1.chrom]):
				continue
			if max(aln2.aln_ref()) >= len(reference[aln2.chrom]):
				continue
			aln1_ref = numpy.array(aln1.aln_ref())
			aln1_read = numpy.array(aln1.aln_read())
			aln1_2write = numpy.logical_and(aln1_ref > flanking,aln1_ref<len(reference[aln1.chrom])-flanking)
			aln1_ref -= flanking
			aln2_ref = numpy.array(aln2.aln_ref())
			aln2_read = numpy.array(aln2.aln_read())
			aln2_2write = numpy.logical_and(aln2_ref > flanking,aln2_ref<len(reference[aln1.chrom])-flanking)
			aln2_ref -= flanking
			if aln1.forward:
				expU_r[aln1.chrom][aln1_ref[aln1_2write],read.seq1[aln1_read[aln1_2write]]] += marginal_p[i]
			else:
				expU_f[aln1.chrom][aln1_ref[aln1_2write],rev_comp(read.seq1)[aln1_read[aln1_2write]]] += marginal_p[i]
			if aln2.forward:
				expU_f[aln2.chrom][aln2_ref[aln2_2write],read.seq2[aln2_read[aln2_2write]]] += marginal_p[i]
			else:
				expU_r[aln2.chrom][aln2_ref[aln2_2write],rev_comp(read.seq2)[aln2_read[aln2_2write]]] += marginal_p[i]
				
	#print 'Finished E step part'
	
	return ll,alned_reads,expU_f,expU_r

# Perform the maximization (M) step for a single repeat
def MStep(input):
	
	#print 'starting M step part'
	
	numpy.seterr(invalid='ignore')
	
	expU_f_seq,expU_r_seq,reference_seq,genome_profile_initial_seq,seq,flanking,prior = input
		
	ll = 0
	pos_type_seq = numpy.zeros(len(reference_seq),dtype=numpy.int8)
	genome_profile_f_seq = copy.copy(genome_profile_initial_seq)
	genome_profile_r_seq = copy.copy(genome_profile_initial_seq)
	
	if len(reference_seq) < 2*flanking:
		return seq, genome_profile_f_seq, genome_profile_r_seq, None, None, ll
	
	genome_profile_type_f = dict()
	genome_profile_type_r = dict()
	ll_vector = dict()
	for type in prior.pseudocounts_plus:
		# expU_plus_pseudo = E_A[U] + \alpha_T - 1, an nx4 matrix, where n is the length of the repeat. The min is used to avoid placing 0 probability on any nucleotide and to deal with negative pseudocounts.
		expU_plus_pseudo_f = expU_f_seq+prior.pseudocounts_plus[type]
		expU_plus_pseudo_f_min0 = numpy.maximum(expU_plus_pseudo_f,numpy.zeros(expU_plus_pseudo_f.shape)+P_ERROR)
		expU_plus_pseudo_r = expU_r_seq+prior.pseudocounts_minus[type]
		expU_plus_pseudo_r_min0 = numpy.maximum(expU_plus_pseudo_r,numpy.zeros(expU_plus_pseudo_r.shape)+P_ERROR)
		# genome_profile_type[T] = (E_A[U] + \alpha_T - 1) / \sum_y(E_A[U] + \alpha_T - 1) = Ghat, a nx4 matrix
		if prior.strands_same[type]:
			genome_profile_type_f[type] = expU_plus_pseudo_f_min0+expU_plus_pseudo_r_min0
			genome_profile_type_f[type] = genome_profile_type_f[type]/numpy.transpose(numpy.tile(numpy.sum(genome_profile_type_f[type],1),(4,1)))
			genome_profile_type_f[type] = numpy.maximum(genome_profile_type_f[type],numpy.zeros(genome_profile_type_f[type].shape)+P_ERROR)
			genome_profile_type_r[type] = genome_profile_type_f[type]
		else:
			genome_profile_type_f[type] = expU_plus_pseudo_f_min0/numpy.transpose(numpy.tile(numpy.sum(expU_plus_pseudo_f_min0,1),(4,1)))
			genome_profile_type_f[type] = numpy.maximum(genome_profile_type_f[type],numpy.zeros(genome_profile_type_f[type].shape)+P_ERROR)
			genome_profile_type_r[type] = expU_plus_pseudo_r_min0/numpy.transpose(numpy.tile(numpy.sum(expU_plus_pseudo_r_min0,1),(4,1)))
			genome_profile_type_r[type] = numpy.maximum(genome_profile_type_r[type],numpy.zeros(genome_profile_type_r[type].shape)+P_ERROR)
		# ll_vector[T] = log P(U=E_A[U], G=Ghat| T) = sum_{nucleotides} log Ghat * (E_A[U] + \alpha_T - 1), an n-long vector.
		ll_vector[type] = numpy.sum( numpy.log(genome_profile_type_f[type])*expU_plus_pseudo_f + numpy.log(genome_profile_type_r[type])*expU_plus_pseudo_r , 1 )
			
	types_for_hyper = dict()
	ll_max_for_hyper = numpy.zeros(len(prior.hyper_dict.keys()))
	for i in range(len(prior.hyper_dict.keys())):
		hyper = prior.hyper_dict.keys()[i]
		ll_matrix = numpy.zeros((len(reference_seq)-2*flanking,len(ll_vector.keys())))
		for j in range(len(ll_vector.keys())):
			type = ll_vector.keys()[j]
			# ll_matrix = log P(H=hyper,T,G=Ghat,U=E_A[U]) = log P(log P(U=E_A[U], G=Ghat| T) + log P(H=hyper) + log P(T|H=hyper), an nxt matrix, where t is number of position types..
			ll_matrix[:,j] = ll_vector[type][flanking:-flanking] + numpy.log(prior.pos_probs[type][hyper][reference_seq[flanking:-flanking]]) - prior.ln_dir_f[type] - prior.ln_dir_r[type]
		ll_max_for_hyper[i] = numpy.sum(numpy.amax(ll_matrix,1))
		# types_for_hyper = argmax_T log P(H=hyper,T,G=Ghat,U=E_A[U]) = T(H)
		types_for_hyper[hyper] = numpy.argmax(ll_matrix,1)
	
	# rep_type_seq = argmax_H log P(H,T=T(H),G=Ghat,U=E_A[U])
	rep_type_seq = prior.hyper_dict.keys()[numpy.argmax(ll_max_for_hyper)]
	ll = max(ll_max_for_hyper)
	pos_type_seq = numpy.array(ll_vector.keys())[types_for_hyper[rep_type_seq]]
	
	for i in range(flanking,len(reference_seq)-flanking):
		genome_profile_f_seq[i,] = genome_profile_type_f[pos_type_seq[i-flanking]][i,]
		genome_profile_r_seq[i,] = genome_profile_type_r[pos_type_seq[i-flanking]][i,]
 		
	return seq, genome_profile_f_seq, genome_profile_r_seq, pos_type_seq, rep_type_seq, ll

# Perform numEM EM steps
def nEMsteps(master_pool,numEM,alignment_pickles,reference,genome_profile_f,genome_profile_r,f_prob,r_prob,genome_profile_initial,prior,flanking,threads):
	db_count = 0
	rep_type = dict()
	pos_type = dict()
		
	for count in range(numEM):
		start = datetime.datetime.now()
		
		Estart = datetime.datetime.now()
		
		inputs = []
		for i in range(len(alignment_pickles)):
			inputs.append((alignment_pickles[i],reference,genome_profile_f,genome_profile_r,f_prob,r_prob,flanking))
		
		ll = 0
		expU_f = dict()
		expU_r = dict()
		for seq in reference:
			expU_f[seq] = numpy.zeros((len(reference[seq]),4))
			expU_r[seq] = numpy.zeros((len(reference[seq]),4))
		alned_reads = 0
		
		# Divide up reads and send to separate processes for E step
		for i in range(len(inputs)/threads):
			EStep_outputs = master_pool.map(EStep, inputs[i*threads:(i+1)*threads])
#			EStep_outputs = list(map(EStep, inputs[i*threads:(i+1)*threads]))
			for EStep_output in EStep_outputs:
				this_ll,this_alned_reads,this_expU_f,this_expU_r = EStep_output
				ll += this_ll
				alned_reads += this_alned_reads
				for seq in reference:
					expU_f[seq][flanking:-flanking] += this_expU_f[seq]
					expU_r[seq][flanking:-flanking] += this_expU_r[seq]
		
		# Combine the output from each process		
		EStep_outputs = master_pool.map(EStep, inputs[len(inputs)/threads*threads:])
#		EStep_outputs = list(map(EStep, inputs[len(inputs)/threads*threads:]))
		for EStep_output in EStep_outputs:
			this_ll,this_alned_reads,this_expU_f,this_expU_r = EStep_output
			ll += this_ll
			alned_reads += this_alned_reads
			for seq in reference:
				expU_f[seq][flanking:-flanking] += this_expU_f[seq]
				expU_r[seq][flanking:-flanking] += this_expU_r[seq]
						
#  		seq = 'FB4_DM_25'
#  		for i in range(len(reference[seq])):
#  			if reference[seq][i] == 0 and expU_f[seq][i,2] > 1:
#  				print seq,i,expU_f[seq][i,]
		
#		print 'Estep',datetime.datetime.now()-Estart
		Mstart = datetime.datetime.now()
		
		#print 'sending M steps'
		
		#MStep_outputs = master_pool.map(MStep, [(expU_f[seq],expU_r[seq],reference[seq],genome_profile_initial[seq],seq,flanking,prior.get_for_seq(seq)) for seq in reference])
		
		# Send the M step for each repeat to a process
		MStep_outputs = []
		for seq in reference:
			MStep_outputs.append(MStep((expU_f[seq],expU_r[seq],reference[seq],genome_profile_initial[seq],seq,flanking,prior.get_for_seq(seq))))
		
		#print 'pooling M step outputs'
		
		# Pool the M step results
		for MStep_output in MStep_outputs:
			seq, genome_profile_f_seq, genome_profile_r_seq, pos_type_seq, rep_type_seq, this_ll = MStep_output
			ll += this_ll
			genome_profile_f[seq] = genome_profile_f_seq
			genome_profile_r[seq] = genome_profile_r_seq
			pos_type[seq] = pos_type_seq
			rep_type[seq] = rep_type_seq
		
		#print 'Finishing EM step'
		
		# f_probs/r_probs = argmax_X P(V=E_A[V] | X)
		f_sum = 0
		for seq in expU_f:
			this_f_sum = numpy.sum(expU_f[seq])/len(expU_f[seq])
			if this_f_sum > 0:
				f_sum += this_f_sum
		for seq in f_prob:
			if len(expU_f[seq]) > 0:
				f_prob[seq] = max(numpy.sum(expU_f[seq])/len(expU_f[seq])/f_sum,10**-10)
			else:
				f_prob[seq] = 0.0
			
		r_sum = 0
		for seq in expU_r:
			this_r_sum = numpy.sum(expU_r[seq])/len(expU_r[seq])
			if this_r_sum > 0:
				r_sum += numpy.sum(expU_r[seq])/len(expU_r[seq])
		for seq in r_prob:
			if len(expU_r[seq]) > 0:
				r_prob[seq] = max(numpy.sum(expU_r[seq])/len(expU_r[seq])/r_sum,10**-10)
			else:
				r_prob[seq] = 0.0
		
#		print 'Mstep',datetime.datetime.now()-Mstart
		
		print count,ll,alned_reads,datetime.datetime.now()-start
	return genome_profile_f,genome_profile_r,f_prob,r_prob,rep_type,pos_type,expU_f,expU_r,ll

# Read command line arguments
def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-a', '--alignment_pickles_list',
							type=str,
							required=False,
							default='alignments_list.txt',
							help='File listing alignment pickles. (alignments_list.txt)')
		parser.add_argument('-r', '--reference_file',
							required=True,
							type=str,
							help='Reference genome in fasta format.')
		parser.add_argument('-p', '--prior_file',
							required=True,
							type=str,
							help='Text file containing prior.')
		parser.add_argument('-j', '--jumpSteps',
							required=False,
							default=5,
							type=int,
							help='Number of EM steps for each jump. 0 means no jumping. Jumping is only implemented for the hyper editing prior. (0)')
		parser.add_argument('-y', '--highlyConfidentOnly',
							required=False,
							default='False',
							type=str,
							help='With j>0, only loop predictions through once and report highly confident predictions.')
		parser.add_argument('-n', '--numEM',
							required=False,
							default=10,
							type=int,
							help='Number of EM steps to take (before jumping.) (10)')
		parser.add_argument('-c', '--coverage',
							required=False,
							default=None,
							type=str,
							help='Coverage of genes containing foldbacks for intial guess. Forward then reverse, comma seperated. File should be repeatname, coverage - tab seperated. (None)')
		parser.add_argument('-i', '--initial_guess',
							required=False,
							default=None,
							type=str,
							help='Specify an initial guess for genome_profile_f, genome_profile_r, f_prob and r_prob - comma separated. Eg: genome_profile_f.pkl,genome_profile_r.pkl,f_prob.pkl,r_prob.pkl. (None)')
		parser.add_argument('-k', '--flanking',
							required=False,
							default=1000,
							type=int,
							help='Amount of fa flanking repeat -- no editing in flank. (1000)')
		parser.add_argument('-t', '--threads',
							required=False,
							default=8,
							type=int,
							help='Number of threads to use. (8)')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.alignment_pickles_list, args.reference_file, args.prior_file, args.jumpSteps, args.highlyConfidentOnly, args.numEM, args.coverage, args.initial_guess, args.flanking, args.threads
	
def main():
	print ' '.join(sys.argv)
	alignment_pickles_list, reference_file, prior_file, jumpSteps, highlyConfidentOnly, numEM, coverage, initial_guess, flanking, threads = GetArgs()
	master_pool = Pool(threads,maxtasksperchild=1)
	
	if highlyConfidentOnly[0].upper() == 'F':
		highlyConfidentOnly = False
	else:
		highlyConfidentOnly = True
		
	# Candidate alignments from parse_PEbam_for_RepProfile.py
	alignment_pickles = []
	for line in open(alignment_pickles_list,'r'):
		alignment_pickles.append(line.strip())
	print alignment_pickles
	
	# Load reference genome into memory
	if reference_file[-4:] == '.pkl':
		reference = pickle.load(open(reference_file,'rb'))
	else:
		reference = fasta_as_array(reference_file)
		pickle.dump(reference,open(reference_file+'.pkl','wb'))
	
	# Read prior
	hyper_dict = dict()
	pos_probs = dict()
	alpha_plus = dict()
	alpha_minus = dict()
	strands_same = dict()
		
	line_type = 'hyper_type'
	for line in open(prior_file,'r'):
		if line.strip()[0] == '#':
			continue
		if line.strip() == 'END':
			if line_type == 'hyper_type':
				line_type = 'pos_type'
			elif line_type == 'pos_type':
				line_type = 'pos_prob'
				for pos_type in alpha_plus:
					pos_probs[pos_type] = dict()
					for hyper_type in hyper_dict:
						pos_probs[pos_type][hyper_type] = numpy.zeros(5)
			else:
				break
			continue
		if line_type == 'hyper_type':
			name,prob = line.strip().split('\t')
			hyper_dict[name] = float(prob)
		elif line_type == 'pos_type':
			name,alpha_p,alpha_m,same = line.strip().split('\t')
			alpha_plus[name] = numpy.array([float(x) for x in alpha_p.split(',')])
			alpha_minus[name] = numpy.array([float(x) for x in alpha_m.split(',')])
			strands_same[name] = bool(int(same))
		elif line_type == 'pos_prob':
			pos_names,TE_names,bases,prob = line.strip().split('\t')
			for pos_type in pos_names.split(','):
				for hyper_type in TE_names.split(','):
					for base in bases.split(','):
						pos_probs[pos_type][hyper_type][int(base)] = float(prob)
		else:
			print 'line_type error:', line_type
	
	pseudocounts_plus = dict()
	for seq in reference:
		pseudocounts_plus[seq] = dict()
		for type in alpha_plus:
			pseudocounts_plus[seq][type] = create_pseudocounts(reference[seq],alpha_plus[type]-1.0)
	pseudocounts_minus = dict()
	for seq in reference:
		pseudocounts_minus[seq] = dict()
		for type in alpha_plus:
			pseudocounts_minus[seq][type] = create_pseudocounts(reference[seq],alpha_minus[type]-1.0)
	prior = dirmixprior(pseudocounts_plus,pseudocounts_minus,strands_same,hyper_dict,pos_probs)
	
	# Initial guess for G
	genome_profile_initial = create_initial_guess(reference,numpy.array([1-3*P_ERROR,P_ERROR,P_ERROR,P_ERROR]))
	pickle.dump(genome_profile_initial,open('genome_profile_initial.pkl','wb'))
	if initial_guess:
		genome_profile_f = pickle.load(open(initial_guess.split(',')[0],'rb'))
		genome_profile_r = pickle.load(open(initial_guess.split(',')[1],'rb'))
	else:
		genome_profile_f = copy.deepcopy(genome_profile_initial)
		genome_profile_r = copy.deepcopy(genome_profile_initial)
	
	# Initial guess for X
	f_prob = dict()
	for seq in reference:
		f_prob[seq] = 1.0
	r_prob = dict()
	for seq in reference:
		r_prob[seq] = 1.0
	if coverage:
		for line in open(coverage.split(',')[0]):
			seq, genecov = line.strip().split('\t')
			f_prob[seq] = float(genecov)
		for line in open(coverage.split(',')[1]):
			seq, genecov = line.strip().split('\t')
			r_prob[seq] = float(genecov)
	if initial_guess:
		f_prob = pickle.load(open(initial_guess.split(',')[2],'rb'))
		r_prob = pickle.load(open(initial_guess.split(',')[3],'rb'))
	f_sum = 0.0
	r_sum = 0.0
	for seq in reference:
		f_sum += f_prob[seq]
		r_sum += r_prob[seq]
	for seq in reference:
		f_prob[seq] = f_prob[seq]/f_sum
		r_prob[seq] = r_prob[seq]/r_sum
	
	if jumpSteps > 0:
		flipped = True
	else:
		flipped = False
	
	# Run EM a few steps and save output
	genome_profile_f,genome_profile_r,f_prob,r_prob,rep_type,pos_type,expU_f,expU_r,ll = nEMsteps(master_pool,numEM,alignment_pickles,reference,genome_profile_f,genome_profile_r,f_prob,r_prob,genome_profile_initial,prior,flanking,threads)
	pickle.dump(genome_profile_f,open('genome_profile_f.pkl','wb'))
	pickle.dump(genome_profile_r,open('genome_profile_r.pkl','wb'))
	pickle.dump(expU_f,open('expU_f.pkl','wb'))
	pickle.dump(expU_r,open('expU_r.pkl','wb'))
	pickle.dump(f_prob,open('f_prob.pkl','wb'))
	pickle.dump(r_prob,open('r_prob.pkl','wb'))
	pickle.dump(rep_type,open('rep_type.pkl','wb'))
	pickle.dump(pos_type,open('pos_type.pkl','wb'))
		
	rep_type_last_rotation = copy.copy(rep_type)
	
	highlyConfident = list()
	
	# Jumping strategy
	while flipped:
		flipped = False
		
		print 'Try flipping these off:'
		print  numpy.array(rep_type.keys())[numpy.array(rep_type.values()) == 'hyper_f'], numpy.array(rep_type.keys())[numpy.array(rep_type.values()) == 'hyper_r']
		
		for flipseq in numpy.append( numpy.array(rep_type.keys())[numpy.array(rep_type.values()) == 'hyper_f'], numpy.array(rep_type.keys())[numpy.array(rep_type.values()) == 'hyper_r'] ):
			
			print flipseq, 'off'
			
			# Remove editing from flipseq and make every position in other repeats predicted to editing have equal probability of A/G
			last_ll = ll
			last_genome_profile_f = copy.deepcopy(genome_profile_f)
			last_genome_profile_r = copy.deepcopy(genome_profile_r)
			last_rep_type = copy.deepcopy(rep_type)
			last_pos_type = copy.deepcopy(pos_type)
			last_f_prob = copy.deepcopy(f_prob)
			last_r_prob = copy.deepcopy(r_prob)
			rep_type[flipseq] = 'not'
			for i in range(flanking,len(reference[flipseq])-flanking):
				if pos_type[flipseq][i-flanking] in ['edit_f','edit_r']:
					genome_profile_f[flipseq][i,] = genome_profile_initial[flipseq][i,]
					genome_profile_r[flipseq][i,] = genome_profile_initial[flipseq][i,]
			for seq in rep_type:
				if rep_type[seq] == 'hyper_f':
					for i in range(len(reference[seq])):
						if reference[seq][i] == 0:
							genome_profile_f[seq][i,0] = (genome_profile_f[seq][i,0]+genome_profile_f[seq][i,2])/2
							genome_profile_f[seq][i,2] = genome_profile_f[seq][i,0]
				if rep_type[seq] == 'hyper_r':
					for i in range(len(reference[seq])):
						if reference[seq][i] == 3:
							genome_profile_r[seq][i,3] = (genome_profile_r[seq][i,3]+genome_profile_r[seq][i,1])/2
							genome_profile_r[seq][i,1] = genome_profile_r[seq][i,3]
				
			# Run a EM steps		
			genome_profile_f,genome_profile_r,f_prob,r_prob,rep_type,pos_type,expU_f,expU_r,ll = nEMsteps(master_pool,jumpSteps,alignment_pickles,reference,genome_profile_f,genome_profile_r,f_prob,r_prob,genome_profile_initial,prior,flanking,threads)
			
			print last_ll, ll, last_rep_type[flipseq] != rep_type[flipseq]
			
			if last_rep_type[flipseq] == rep_type[flipseq]:
				highlyConfident.append(flipseq)
			
			# Keep the better estimate
			if last_ll > ll:
				genome_profile_f = last_genome_profile_f
				genome_profile_r = last_genome_profile_r
				f_prob = last_f_prob
				r_prob = last_r_prob
				rep_type = last_rep_type
				pos_type = last_pos_type
				ll = last_ll
			elif last_rep_type[flipseq] != rep_type[flipseq]:
				print 'Flipped', flipseq, 'off'
						
		pickle.dump(genome_profile_f,open('genome_profile_f.pkl','wb'))
		pickle.dump(genome_profile_r,open('genome_profile_r.pkl','wb'))
		pickle.dump(expU_f,open('expU_f.pkl','wb'))
		pickle.dump(expU_r,open('expU_r.pkl','wb'))
		pickle.dump(f_prob,open('f_prob.pkl','wb'))
		pickle.dump(r_prob,open('r_prob.pkl','wb'))
		pickle.dump(rep_type,open('rep_type.pkl','wb'))
		pickle.dump(pos_type,open('pos_type.pkl','wb'))
		
		#print rep_type
		#print rep_type_last_rotation
		
		# Stop when no hyper-editing is flipped off
		for seq in rep_type:
			if rep_type[seq] != rep_type_last_rotation[seq] and not highlyConfidentOnly:
				flipped=True
		rep_type_last_rotation = copy.copy(rep_type)
	
	if jumpSteps > 0:
		genome_profile_f,genome_profile_r,f_prob,r_prob,rep_type,pos_type,expU_f,expU_r,ll = nEMsteps(master_pool,numEM,alignment_pickles,reference,genome_profile_f,genome_profile_r,f_prob,r_prob,genome_profile_initial,prior,flanking,threads)
		pickle.dump(genome_profile_f,open('genome_profile_f.pkl','wb'))
		pickle.dump(genome_profile_r,open('genome_profile_r.pkl','wb'))
		pickle.dump(expU_f,open('expU_f.pkl','wb'))
		pickle.dump(expU_r,open('expU_r.pkl','wb'))
		pickle.dump(f_prob,open('f_prob.pkl','wb'))
		pickle.dump(r_prob,open('r_prob.pkl','wb'))
		pickle.dump(rep_type,open('rep_type.pkl','wb'))
		pickle.dump(pos_type,open('pos_type.pkl','wb'))
		if highlyConfidentOnly:
			print highlyConfident

if __name__ == '__main__':
	main()
