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
import pysam

# ERROR AND INDEL PARAMETERS
P_ERROR = 0.001 # Probability of a given sequencing error (e.g. A->C)
P_OPEN = 0.001 # Probability of opening a new indel
P_EXTEND = 0.1 # Probability of extending an indel

class readclass(object):
	def __init__(self, read,rname_list):
		self.seq = seq2array(read.seq)
		self.aln = [alnclass(rname_list[read.rname],read.pos,read.cigarstring,not read.is_reverse)]
		self.copies = 1
	def addread(self):
		self.copies += 1
	def addaln(self,XA):
		self.aln += [alnclass(XA[0],int(XA[1][1:])-1,XA[2],XA[1][0]=='+')]
		
class alnclass(object):
	def __init__(self,chrom,ref_start,cigarstring,forward):
		self.chrom = chrom
		self.forward = forward
		self.read_ranges,self.ref_ranges,self.insert_starts,self.insert_extends,self.del_starts,self.del_extends = cigarstring_2_pairs(cigarstring,ref_start)
	def aln_read(self):
		result = []
		for x in self.read_ranges:
			result += range(x[0],x[1])
		return result
		#return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.read_ranges])))
	def aln_ref(self):
		result = []
		for x in self.ref_ranges:
			result += range(x[0],x[1])
		return result
		#return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.ref_ranges])))
		
class PEreadclass(object):
	def __init__(self, read1,read2,aln1,aln2):
		self.seq1 = read1.seq
		self.seq2 = read2.seq
		self.alns1 = [aln1]
		self.alns2 = [aln1]
		self.aln = [(1,1)]
	def addaln1(self,aln1):
		self.alns1.append( aln1 )
	def addaln2(self,aln2):
		self.alns2.append( aln2 )
	def addaln(self,i,j):
		self.aln.append( (i,j) )

def alnpairs_to_cigarstring(pairs):
	cigar = list()
	last_pair = None
	for pair in pairs:
		if not last_pair:
			if pair[1] != 0:
				cigar.append(['I',pair[1]])
			cigar.append(['M',1])
		else:
			if pair[0] - last_pair[0] > 1:
				cigar.append(['D',pair[0] - last_pair[0]-1])
			if pair[1] - last_pair[1] > 1:
				cigar.append(['I',pair[1] - last_pair[1]-1])
			if cigar[-1][0] == 'M':
				cigar[-1][1] += 1
			else:
				cigar.append(['M',1])
		last_pair = pair
	cigarstring = ''
	for pair in cigar:
		cigarstring += str(pair[1])+pair[0]
	return cigarstring

def rev_comp(seq_array):
	revcompseq = 3 - seq_array[::-1]
	revcompseq[revcompseq<0] = 4
	return revcompseq

def array2seq(seq_array):
	seq = ''
	for letter in seq_array:
		seq += 'ACGTN'[letter]
	return seq

def write_sam(alignment_pickle,genome_profile_f,genome_profile_r,f_prob,r_prob):

	alignments = pickle.load(open(alignment_pickle,'rb'))
	
	for read_id in alignments:
		read = alignments[read_id]
		if len(read.aln) == 0:
			continue
		unnorm_p1 = numpy.zeros(len(read.alns1))
		for i in range(len(read.alns1)):
			aln1 = read.alns1[i]
			if max(aln1.aln_ref()) >= len(genome_profile_f[aln1.chrom]):
				continue
			if aln1.forward:
				unnorm_p1[i] = numpy.prod(genome_profile_r[aln1.chrom][aln1.aln_ref(),read.seq1[aln1.aln_read()]])*r_prob[aln1.chrom]*P_OPEN**(aln1.insert_starts+aln1.del_starts)*P_EXTEND**(aln1.insert_extends+aln1.del_extends)*0.25**(aln1.insert_starts+aln1.insert_extends)
			else:
				unnorm_p1[i] = numpy.prod(genome_profile_f[aln1.chrom][aln1.aln_ref(),rev_comp(read.seq1)[aln1.aln_read()]])*f_prob[aln1.chrom]*P_OPEN**(aln1.insert_starts+aln1.del_starts)*P_EXTEND**(aln1.insert_extends+aln1.del_extends)*0.25**(aln1.insert_starts+aln1.insert_extends)
		
		unnorm_p2 = numpy.zeros(len(read.alns2))
		for j in range(len(read.alns2)):
			aln2 = read.alns2[j]
			if max(aln2.aln_ref()) >= len(genome_profile_f[aln2.chrom]):
				continue
			if aln2.forward:
				unnorm_p2[j] = numpy.prod(genome_profile_f[aln2.chrom][aln2.aln_ref(),read.seq2[aln2.aln_read()]])*f_prob[aln2.chrom]*P_OPEN**(aln2.insert_starts+aln2.del_starts)*P_EXTEND**(aln2.insert_extends+aln2.del_extends)*0.25**(aln2.insert_starts+aln2.insert_extends)
			else:
				unnorm_p2[j] = numpy.prod(genome_profile_r[aln2.chrom][aln2.aln_ref(),rev_comp(read.seq2)[aln2.aln_read()]])*r_prob[aln2.chrom]*P_OPEN**(aln2.insert_starts+aln2.del_starts)*P_EXTEND**(aln2.insert_extends+aln2.del_extends)*0.25**(aln2.insert_starts+aln2.insert_extends)
			
		unnorm_p = numpy.zeros(len(read.aln))
		for i in range(len(read.aln)):
			unnorm_p[i] = unnorm_p1[read.aln[i][0]]*unnorm_p2[read.aln[i][1]]
		if sum(unnorm_p)>0:
			marginal_p = unnorm_p/sum(unnorm_p)
		else:
			continue
		
		i = numpy.argmax(marginal_p)
		aln1 = read.alns1[read.aln[i][0]]
		aln2 = read.alns2[read.aln[i][1]]
		flag1 = 1+2+16*(not aln1.forward)+32*(not aln2.forward)+64
		flag2 = 1+2+16*(not aln2.forward)+32*(not aln1.forward)+128
		mapq = int(min(-10*numpy.log(1.0-marginal_p[i])/numpy.log(10),30))
		cigar1 = alnpairs_to_cigarstring(zip(aln1.aln_ref(),aln1.aln_read()))
		cigar2 = alnpairs_to_cigarstring(zip(aln2.aln_ref(),aln2.aln_read()))
		frag_length = max(aln1.aln_ref()+aln2.aln_ref()) - min(aln1.aln_ref()+aln2.aln_ref()) + 1
		if aln1.forward:
			oriented_seq1 = array2seq(read.seq1)
		else:
			oriented_seq1 = array2seq(rev_comp(read.seq1))
		if aln2.forward:
			oriented_seq2 = array2seq(read.seq2)
		else:
			oriented_seq2 = array2seq(rev_comp(read.seq2))
		print read_id +'\t'+ str(flag1) +'\t'+ aln1.chrom +'\t'+ str(min(aln1.aln_ref())+1) +'\t'+ str(mapq) +'\t'+ cigar1 +'\t'+ aln2.chrom +'\t'+ str(min(aln2.aln_ref())+1) +'\t'+ str(frag_length) +'\t'+ oriented_seq1 +'\t*'
		print read_id +'\t'+ str(flag2) +'\t'+ aln2.chrom +'\t'+ str(min(aln2.aln_ref())+1) +'\t'+ str(mapq) +'\t'+ cigar1 +'\t'+ aln1.chrom +'\t'+ str(min(aln1.aln_ref())+1) +'\t'+ str(frag_length) +'\t'+ oriented_seq2 +'\t*' 
	
	return

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
							help='File listing alignment pickles.')
		parser.add_argument('-f', '--prob_pickles',
							required=False,
							default='f_prob.pkl,r_prob.pkl',
							type=str,
							help='f_prob,r_prob.')
		parser.add_argument('-t', '--genome_profile',
							required=False,
							default='genome_profile_f.pkl,genome_profile_r.pkl',
							type=str,
							help='Profile to align to. Forward then reverse, comma seperated. Eg: genome_profile_f.pkl,genome_profile_r.pkl')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.alignment_pickles_list, args.prob_pickles, args.genome_profile
	
def main():
	alignment_pickles_list, prob_pickles, genome_profile = GetArgs()
			
	alignment_pickles = []
	for line in open(alignment_pickles_list,'r'):
		alignment_pickles.append(line.strip())
	
	genome_profile_f = pickle.load(open(genome_profile.split(',')[0],'rb'))
	genome_profile_r = pickle.load(open(genome_profile.split(',')[1],'rb'))
	
	f_prob = pickle.load(open(prob_pickles.split(',')[0],'rb'))
	r_prob = pickle.load(open(prob_pickles.split(',')[1],'rb'))
	
	
	for alignment_pickle in alignment_pickles:
		write_sam(alignment_pickle,genome_profile_f,genome_profile_r,f_prob,r_prob)

if __name__ == '__main__':
	main()
