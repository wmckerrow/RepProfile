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

import pickle
import numpy
from Bio import SeqIO
import argparse
import sys
import random

# HYPERPARAMETERS
P_SNP = 0.01 # Probability of SNP at a given position
P_EDIT = 0.5 # Probability that an A(T) is edited if element is hyper edited
P_ERROR = 0.001 # Probability of a given sequencing error (e.g. A->C)
P_OPEN = 0.0 # Probability of opening a new indel. Not implemented.
P_EXTEND = 0.0 # Probability of extending an indel. Not implemented.
FLANKING = 1000

def seq2array(seq):
	nuc2num = {'A':0,'C':1,'G':2,'T':3, 'N':4}
	return numpy.array([nuc2num[c] for c in seq])

def array2seq(seq_array):
	seq = ''
	for letter in seq_array:
		seq += 'ACGTN'[letter]
	return seq

def fasta_as_array(fasta_file):
	fasta_dict = SeqIO.to_dict(SeqIO.parse(open(fasta_file,'r'), "fasta"))
	for seq in fasta_dict:
		fasta_dict[seq] = seq2array(fasta_dict[seq].upper())
	return fasta_dict
	
def rev_comp(seq_array):
	revcompseq = 3 - seq_array[::-1]
	revcompseq[revcompseq<0] = 4
	return revcompseq

def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-r', '--reference_file',
							type=str,
							required=True,
							help='Reference in fasta format.')
		parser.add_argument('-p', '--hyper_editing_prob',
							type=float,
							required=True,
							help='Probability that a given editable TE is edited.')
		parser.add_argument('-e', '--editable',
							type=str,
							required=True,
							help='Editable TEs <tab> edit strand.')
		parser.add_argument('-v', '--coverage',
							required=True,
							type=str,
							help='Coverage of genes containing foldbacks for intial guess. Forward then reverse, comma seperated. File should be repeatname, coverage - tab seperated.')
		parser.add_argument('-n', '--n_reads',
							required=True,
							type=int,
							help='Number of reads to simulate.')
		parser.add_argument('-l', '--read_length',
							required=False,
							type=int,
							default = 100,
							help='Length of simulates reads.')
		parser.add_argument('-m','--insert_mean',
							required=False,
							type=int,
							default = 200,
							help='Mean insert size.')
		parser.add_argument('-s','--insert_sd',
							required=False,
							type=int,
							default = 40,
							help='Standard deviation of insert size.')
		parser.add_argument('-f','--fastq_prefix',
							required=False,
							type=str,
							default = 'reads',
							help='Prefix for fastq files.')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.reference_file, args.hyper_editing_prob, args.editable, args.coverage, args.n_reads, args.read_length, args.insert_mean, args.insert_sd, args.fastq_prefix
	
def Main():
	reference_file, hyper_editing_prob, sense_dsrna, coverage, n_reads, read_length, insert_mean, insert_sd, fastq_prefix = GetArgs()
	reference = fasta_as_array(reference_file)
	
	hyper_edited_forward = []
	hyper_edited_reverse = []
	for line in open(sense_dsrna,'r'):
		TE, strand = line.strip().split('\t')
		if strand == '+' and random.random() < hyper_editing_prob:
			hyper_edited_forward.append(TE)
		elif strand == '-' and random.random() < hyper_editing_prob:
			hyper_edited_reverse.append(TE)
			
	pickle.dump(hyper_edited_forward,open('hyper_edited_forward.pkl','wb'))
	pickle.dump(hyper_edited_reverse,open('hyper_edited_reverse.pkl','wb'))
		
	fastq1 = open(fastq_prefix+'_R1.fastq','w')
	fastq2 = open(fastq_prefix+'_R2.fastq','w')
	
	f_prob = dict()
	for seq in reference:
		f_prob[seq] = 1.0
	r_prob = dict()
	for seq in reference:
		r_prob[seq] = 1.0
	for line in open(coverage.split(',')[0]):
		seq, genecov = line.strip().split('\t')
		f_prob[seq] = float(genecov)
	for line in open(coverage.split(',')[1]):
		seq, genecov = line.strip().split('\t')
		r_prob[seq] = float(genecov)
	for seq in reference:
		if len(reference[seq]) <= read_length:
			f_prob[seq] = 0.0
			r_prob[seq] = 0.0
		
	f_sum = 0
	for seq in f_prob:
		f_sum += f_prob[seq]
	for seq in f_prob:
		f_prob[seq] = f_prob[seq]/f_sum
	r_sum = 0
	for seq in r_prob:
		r_sum += r_prob[seq]
	for seq in r_prob:
		r_prob[seq] = r_prob[seq]/r_sum
		
	# 0=ref, 1=snp, 2=f_edit, 3=r_edit
	pos_type = dict()
	for seq in reference:
		pos_type[seq] = numpy.zeros(len(reference[seq]),dtype=int)
		for i in range(len(reference[seq])):
			if seq in hyper_edited_forward and reference[seq][i] == 0 and i >= FLANKING and i <= len(reference[seq])-FLANKING:
				pos_type[seq][i] = numpy.random.choice(range(4),p=[1.0-P_SNP-P_EDIT,P_SNP,P_EDIT,0.0])
			elif seq in hyper_edited_reverse and reference[seq][i] == 3 and i >= FLANKING and i <= len(reference[seq])-FLANKING:
				pos_type[seq][i] = numpy.random.choice(range(4),p=[1.0-P_SNP-P_EDIT,P_SNP,0.0,P_EDIT])
			else:
				pos_type[seq][i] = numpy.random.choice(range(4),p=[1.0-P_SNP,P_SNP,0.0,0.0])
				
	theta_f = dict()
	theta_r = dict()
	
	for seq in reference:
		theta_f[seq] = numpy.zeros((len(reference[seq]),4))
		theta_r[seq] = numpy.zeros((len(reference[seq]),4))
		
		for i in range(len(reference[seq])):
			if reference[seq][i] == 4:
				theta_f[seq][i,] = numpy.array([0.25,0.25,0.25,0.25])
				theta_r[seq][i,] = numpy.array([0.25,0.25,0.25,0.25])
				continue
			if pos_type[seq][i] == 0:
				theta_f[seq][i,] = numpy.array([1.0-3*P_ERROR,P_ERROR,P_ERROR,P_ERROR])[range(4)-reference[seq][i]]
				theta_r[seq][i,] = numpy.array([1.0-3*P_ERROR,P_ERROR,P_ERROR,P_ERROR])[range(4)-reference[seq][i]]
			elif pos_type[seq][i] == 1:
				theta_f[seq][i,] = numpy.array([P_ERROR,P_ERROR,P_ERROR,P_ERROR])
				p = P_ERROR+numpy.random.rand()*(1.0-4.0*P_ERROR)
				theta_f[seq][i,reference[seq][i]] = p
				theta_f[seq][i,numpy.random.choice(numpy.array(range(reference[seq][i])+range(reference[seq][i]+1,4)))] = 1.0-p-2.0*P_ERROR
				theta_r[seq][i,] = theta_f[seq][i,]
			elif pos_type[seq][i] == 2:
				p = numpy.random.rand()
				theta_f[seq][i,] = numpy.array([(1.0-2*P_ERROR)*(1.0-p),P_ERROR,(1.0-2*P_ERROR)*p,P_ERROR])
				theta_r[seq][i,] = numpy.array([1.0-3*P_ERROR,P_ERROR,P_ERROR,P_ERROR])
			else:
				p = numpy.random.rand()
				theta_f[seq][i,] = numpy.array([P_ERROR,P_ERROR,P_ERROR,1.0-3*P_ERROR])
				theta_r[seq][i,] = numpy.array([P_ERROR,(1.0-2*P_ERROR)*p,P_ERROR,(1.0-2*P_ERROR)*(1.0-p)])
				
	for read_count in range(n_reads):
		is_forward = numpy.random.rand() < f_sum/(f_sum+r_sum)
		insert_size = max(int(numpy.random.normal(insert_mean,insert_sd)),0)
		if is_forward:
			seq = numpy.random.choice(f_prob.keys(),p=f_prob.values())
			pos = numpy.random.randint(len(reference[seq])-(2*read_length + insert_size))
			read_array = numpy.zeros(read_length,dtype=int)
			for i in range(read_length):
				read_array[i] = numpy.random.choice(range(4),p=theta_f[seq][pos+i,])
		else:
			seq = numpy.random.choice(r_prob.keys(),p=r_prob.values())
			pos = numpy.random.randint(len(reference[seq])-(2*read_length + insert_size))
			read_array = numpy.zeros(read_length,dtype=int)
			for i in range(read_length):
				read_array[i] = numpy.random.choice(range(4),p=theta_r[seq][pos+i,])
			read_array = rev_comp(read_array)
		read_name = 'read_' + str(read_count) + ',' + seq+':'+str(pos)+'-'+str(pos+insert_size+2*read_length) + ',' + '-+'[int(is_forward)]
		fastq2.write( '@'+read_name+'\n')
		fastq2.write(  array2seq(read_array) +'\n')
		fastq2.write(  '+' +'\n')
		fastq2.write(  'A'*read_length +'\n')
		if is_forward:
			pos += read_length + insert_size
			read_array = numpy.zeros(read_length,dtype=int)
			for i in range(read_length):
				read_array[i] = numpy.random.choice(range(4),p=theta_f[seq][pos+i,])
			read_array = rev_comp(read_array)
		else:
			pos += read_length + insert_size
			read_array = numpy.zeros(read_length,dtype=int)
			for i in range(read_length):
				read_array[i] = numpy.random.choice(range(4),p=theta_r[seq][pos+i,])
		fastq1.write( '@'+read_name+'\n')
		fastq1.write(  array2seq(read_array) +'\n')
		fastq1.write(  '+' +'\n')
		fastq1.write(  'A'*read_length +'\n')
		
	pickle.dump(pos_type,open('pos_type_sim.pkl','wb'))
	pickle.dump(theta_f,open('theta_f_sim.pkl','wb'))
	pickle.dump(theta_r,open('theta_r_sim.pkl','wb'))
	
Main()
