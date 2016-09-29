import pysam
import argparse
import pickle
import numpy
from Bio import SeqIO
import re
from itertools import chain
from sets import Set

MAX_FRAG_SIZE = 1000

def aln_overlap(aln1,aln2):
	if aln1.chrom == aln2.chrom and min(aln1.aln_ref()) < max(aln2.aln_ref()) and min(aln2.aln_ref()) < max(aln1.aln_ref()):
		return True
	else:
		return False

class readclass(object):
	def __init__(self, read,rname_list,readfq):
		self.seq = seq2array(str(readfq[read.qname].seq).upper())
		self.aln = [alnclass(rname_list[read.rname],read.pos,read.cigarstring,not read.is_reverse)]
	def addaln(self,read,rname_list):
		new_aln = alnclass(rname_list[read.rname],read.pos,read.cigarstring,not read.is_reverse)
		for old_aln in self.aln:
			if aln_overlap(old_aln,new_aln):
				return False
		self.aln += [new_aln]
		return True
	def addXAaln(self,XA):
		new_aln = alnclass(XA[0],int(XA[1][1:])-1,XA[2],XA[1][0]=='+')
		for old_aln in self.aln:
			if aln_overlap(old_aln,new_aln):
				return False
		self.aln += [new_aln]
		return True
		
class alnclass(object):
	def __init__(self,chrom,ref_start,cigarstring,forward):
		self.chrom = chrom
		self.forward = forward
		self.read_ranges,self.ref_ranges,self.insert_starts,self.insert_extends,self.del_starts,self.del_extends = cigarstring_2_pairs(cigarstring,ref_start)
	def aln_read(self):
		return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.read_ranges])))
	def aln_ref(self):
		return numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in self.ref_ranges])))
		
class PEreadclass(object):
	def __init__(self, read1,read2,aln1,aln2):
		self.seq1 = read1.seq
		self.seq2 = read2.seq
		self.alns1 = [aln1]
		self.alns2 = [aln2]
		self.aln = [(0,0)]
	def addaln1(self,aln1):
		self.alns1.append( aln1 )
	def addaln2(self,aln2):
		self.alns2.append( aln2 )
	def addaln(self,i,j):
		self.aln.append( (i,j) )
		
def cigarstring_2_pairs(cigarstring,ref_start):
	aln_read = []
	aln_ref = []
	read_ranges = []
	ref_ranges = []
	insert_starts = 0
	insert_extends = 0
	del_starts = 0
	del_extends = 0
	seq_pos = 0
	ref_pos = ref_start
	
	if cigarstring[0] in '1234567890':
		for cigar in re.findall('(\d+[DMI])', cigarstring):
			if cigar[-1] == 'M':
				aln_ref += range(ref_pos,ref_pos+int(cigar[:-1]))
				ref_ranges.append((ref_pos,ref_pos+int(cigar[:-1])))
				ref_pos += int(cigar[:-1])
				aln_read += range(seq_pos,seq_pos+int(cigar[:-1]))
				read_ranges.append((seq_pos,seq_pos+int(cigar[:-1])))
				seq_pos += int(cigar[:-1])
			elif cigar[-1] == 'D':
				ref_pos += int(cigar[:-1])
				del_starts += 1
				del_extends += int(cigar[:-1])-1
			elif cigar[-1] == 'I':
				seq_pos += int(cigar[:-1])
				insert_starts += 1
				insert_extends += int(cigar[:-1]) -1
	else:
		for cigar in re.findall('([DMI]\d+)', cigarstring):
			if cigar[0] == 'M':
				aln_ref += range(ref_pos,ref_pos+int(cigar[1:]))
				ref_ranges.append((ref_pos,ref_pos+int(cigar[1:])))
				ref_pos += int(cigar[1:])
				aln_read += range(seq_pos,seq_pos+int(cigar[1:]))
				read_ranges.append((seq_pos,seq_pos+int(cigar[1:])))
				seq_pos += int(cigar[1:])
			elif cigar[0] == 'D':
				ref_pos += int(cigar[1:])
				del_starts += 1
				del_extends += int(cigar[1:])-1
			elif cigar[0] == 'I':
				seq_pos += int(cigar[1:])
				insert_starts += 1
				insert_extends += int(cigar[1:]) -1
				
	if not numpy.array_equal( numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in read_ranges]))) , aln_read ):
		print read_ranges
		print numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in read_ranges])))
		print aln_read
		exit()
	if not numpy.array_equal( numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in ref_ranges]))) , aln_ref ):
		print ref_ranges
		print numpy.array(list(chain.from_iterable([range(x[0],x[1]) for x in ref_ranges])))
		print aln_ref
		exit()
				
	return read_ranges,ref_ranges,insert_starts,insert_extends,del_starts,del_extends

def seq2array(seq):
	nuc2num = {'A':0,'C':1,'G':2,'T':3, 'N':4}
	return numpy.array([nuc2num[c] for c in seq],dtype=numpy.int8)


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
	
def check_alignments(read,ref,max_diff):
	diffs = sum(read != ref)
	edits = sum(numpy.logical_and(ref==0,read==2))
	return diffs-edits <= max_diff

def count_diffs(read,ref):
	diffs = sum(read != ref)
	edits = sum(numpy.logical_and(ref==0,read==2))
	return edits,diffs-edits

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

		parser.add_argument('-b', '--bam_file',
							type=str,
							required=True,
							help='Bam file with all alignments. Sorted by read name.')
		parser.add_argument('-p', '--prefix',
							type=str,
							required=False,
							default='alignments',
							help='Prefix for alignment files.')
		parser.add_argument('-r', '--read_fastqs',
							type=str,
							required=True,
							help='Fastq files for reads, comma seperated. Bamfiles should contain disjoint reads.')
		parser.add_argument('-g', '--reference_file',
							required=True,
							type=str,
							help='Reference in fasta format.')
		parser.add_argument('-m', '--max_diff',
							required=True,
							type=int,
							help='Maximum number of non A->G changes allowed.')
		parser.add_argument('-k', '--flanking',
							required=False,
							default=1000,
							type=int,
							help='Length of extra sequence flanking TE.')
		parser.add_argument('-n', '--reads_per_pickle',
							required=False,
							default=5000,
							type=int,
							help='Number of reads in each alignment pickle.')
		parser.add_argument('-q', '--qcutoff',
							required=False,
							default=30,
							type=int,
							help='Remove reads if any base is below the quality cutoff.')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.bam_file, args.prefix, args.read_fastqs, args.reference_file, args.max_diff, args.flanking, args.reads_per_pickle, args.qcutoff
	
def Main():
	bam_file, prefix, read_fastqs, reference_file, max_diff, flanking, reads_per_pickle, qcutoff = GetArgs()
	reference = fasta_as_array(reference_file)
	
	readfq1 = SeqIO.to_dict(SeqIO.parse(open(read_fastqs.split(',')[0],'r'), "fastq"))
	readfq2 = SeqIO.to_dict(SeqIO.parse(open(read_fastqs.split(',')[1],'r'), "fastq"))
	
	read_count = 0
	alns_created = 0
	alns_kept = 0
	
	reads = dict()
	pickle_num = 0
	
	bam = pysam.Samfile(bam_file, "rb")
	rnames = bam.references
	
	read_id = None
	
	alignments_list_file = open('alignments_list.txt','w')
	
	for read in bam:
		if read.is_unmapped:
			continue
		if 'N' in read.cigarstring or 'S' in read.cigarstring or 'H' in read.cigarstring or 'P' in read.cigarstring or '=' in read.cigarstring or 'X' in read.cigarstring:
			continue
		if 'N' in str((readfq1,readfq2)[read.is_read2][read.qname].seq).upper():
			continue
		if min(read.query_qualities) < qcutoff:
			continue
			
		if not read_id:
			read_id = read.qname
			new_read_id1 = True
			new_read_id2 = True
			
		if read_id != read.qname:
			if not (new_read_id1 or new_read_id2):
				accepted_read = None
				kept1_pos = [-1]*len(read1.aln)
				kept2_pos = [-1]*len(read2.aln)
				for i in range(len(read1.aln)):
					aln1 = read1.aln[i]
					for j in range(len(read2.aln)):
						aln2 = read2.aln[j]
						if max(max(aln1.aln_ref()),max(aln2.aln_ref())) < flanking or min(min(aln1.aln_ref()),min(aln2.aln_ref())) > len(reference[aln1.chrom])-flanking and aln1.forward != aln2.forward:
							continue
						if (aln1.chrom == aln2.chrom) and (max(max(aln1.aln_ref()),max(aln2.aln_ref()))-min(min(aln1.aln_ref()),min(aln2.aln_ref())) < MAX_FRAG_SIZE):
							if not accepted_read:
								accepted_read = PEreadclass(read1,read2,aln1,aln2)
								alns_kept += 1
								kept1_pos[i] = max(kept1_pos)+1
								kept2_pos[j] = max(kept2_pos)+1
							else:
								if kept1_pos[i] == -1:
									accepted_read.addaln1(aln1)
									kept1_pos[i] = max(kept1_pos)+1
								if kept2_pos[j] == -1:
									accepted_read.addaln2(aln2)
									kept2_pos[j] = max(kept2_pos)+1
								accepted_read.addaln(kept1_pos[i],kept2_pos[j])
								alns_kept += 1
				if accepted_read:
					reads[read_id] = accepted_read
				if len(reads) >= reads_per_pickle:
					pickle.dump(reads,open(prefix+'_'+str(pickle_num)+'.pkl','wb'))
					alignments_list_file.write(prefix+'_'+str(pickle_num)+'.pkl\n')
					pickle_num += 1
					reads = dict()
			
			read_id = read.qname
			new_read_id1 = True
			new_read_id2 = True
			
		read_count+=1
		
		candidate = readclass(read,rnames,(readfq1,readfq2)[read.is_read2])
		if read.is_read2 and not read.is_reverse:
			if check_alignments(candidate.seq[candidate.aln[0].aln_read()],reference[candidate.aln[0].chrom][candidate.aln[0].aln_ref()],max_diff):
				if new_read_id2:
					read2 = candidate
					new_read_id2 = False
					alns_created += 1
				else:
					alns_created += read2.addaln(read,rnames)
		elif read.is_read1 and read.is_reverse:
			if check_alignments(rev_comp(candidate.seq[candidate.aln[0].aln_read()]),reference[candidate.aln[0].chrom][candidate.aln[0].aln_ref()],max_diff):
				if new_read_id1:
					read1 = candidate
					new_read_id1 = False
					alns_created += 1
				else:
					alns_created += read1.addaln(read,rnames)
		elif read.is_read2 and read.is_reverse:
			if check_alignments(candidate.seq[candidate.aln[0].aln_read()],rev_comp(reference[candidate.aln[0].chrom][candidate.aln[0].aln_ref()]),max_diff):
				if new_read_id2:
					read2 = candidate
					new_read_id2 = False
					alns_created += 1
				else:
					read2.addaln(read,rnames)
					alns_created += 1
		elif read.is_read1 and not read.is_reverse:
			if check_alignments(rev_comp(candidate.seq[candidate.aln[0].aln_read()]),rev_comp(reference[candidate.aln[0].chrom][candidate.aln[0].aln_ref()]),max_diff):
				if new_read_id1:
					read1 = candidate
					new_read_id1 = False
					alns_created += 1
				else:
					read1.addaln(read,rnames)
					alns_created += 1
					
		if 'XA' in dict(read.tags):
			for aln in [x.split(',') for x in dict(read.tags)['XA'].split(';')[:-1]]:
				read_count+=1
				if not ('N' in aln[2] or 'S' in aln[2] or 'H' in aln[2] or 'P' in aln[2] or '=' in aln[2] or 'X' in aln[2]):
					candidate.addXAaln(aln)
					if max(candidate.aln[-1].aln_ref()) >= len(reference[candidate.aln[-1].chrom]):
						continue
					if read.is_read2 and candidate.aln[-1].forward:
						if check_alignments(candidate.seq[candidate.aln[-1].aln_read()],reference[candidate.aln[-1].chrom][candidate.aln[-1].aln_ref()],max_diff):
							if new_read_id2:
								read2 = candidate
								new_read_id2 = False
								candidate.aln = [candidate.aln[-1]]
								alns_created += 1
							else:
								read2.addXAaln(aln)
								alns_created += 1
					elif read.is_read1 and not candidate.aln[-1].forward:
						if check_alignments(rev_comp(candidate.seq[candidate.aln[-1].aln_read()]),reference[candidate.aln[-1].chrom][candidate.aln[-1].aln_ref()],max_diff):
							if new_read_id1:
								read1 = candidate
								new_read_id1 = False
								candidate.aln = [candidate.aln[-1]]
								alns_created += 1
							else:
								read1.addXAaln(aln)
								alns_created += 1
					elif read.is_read2 and read.is_reverse:
						if check_alignments(candidate.seq[candidate.aln[-1].aln_read()],rev_comp(reference[candidate.aln[-1].chrom][candidate.aln[-1].aln_ref()]),max_diff):
							if new_read_id2:
								read2 = candidate
								new_read_id2 = False
								candidate.aln = [candidate.aln[-1]]
								alns_created += 1
							else:
								read2.addXAaln(aln)
								alns_created += 1
					elif read.is_read1 and not read.is_reverse:
						if check_alignments(rev_comp(candidate.seq[candidate.aln[-1].aln_read()]),rev_comp(reference[candidate.aln[-1].chrom][candidate.aln[-1].aln_ref()]),max_diff):
							if new_read_id1:
								read1 = candidate
								new_read_id1 = False
								candidate.aln = [candidate.aln[-1]]
								alns_created += 1
							else:
								read1.addXAaln(aln)
								alns_created += 1
	
	if not (new_read_id1 or new_read_id2):
		accepted_read = None
		kept1_pos = [-1]*len(read1.aln)
		kept2_pos = [-1]*len(read2.aln)
		for i in range(len(read1.aln)):
			aln1 = read1.aln[i]
			for j in range(len(read2.aln)):
				aln2 = read2.aln[j]
				if max(max(aln1.aln_ref()),max(aln2.aln_ref())) < flanking or min(min(aln1.aln_ref()),min(aln2.aln_ref())) > len(reference[aln1.chrom])-flanking and aln1.forward != aln2.forward:
					continue
				if (aln1.chrom == aln2.chrom) and (max(max(aln1.aln_ref()),max(aln2.aln_ref()))-min(min(aln1.aln_ref()),min(aln2.aln_ref())) < MAX_FRAG_SIZE):
					if not accepted_read:
						accepted_read = PEreadclass(read1,read2,aln1,aln2)
						alns_kept += 1
						kept1_pos[i] = max(kept1_pos)+1
						kept2_pos[j] = max(kept2_pos)+1
					else:
						if kept1_pos[i] == -1:
							accepted_read.addaln1(aln1)
							kept1_pos[i] = max(kept1_pos)+1
						if kept2_pos[j] == -1:
							accepted_read.addaln2(aln2)
							kept2_pos[j] = max(kept2_pos)+1
						accepted_read.addaln(kept1_pos[i],kept2_pos[j])
						alns_kept += 1
		if accepted_read:
			reads[read_id] = accepted_read
	
	pickle.dump(reads,open(prefix+'_'+str(pickle_num)+'.pkl','wb'))
	alignments_list_file.write(prefix+'_'+str(pickle_num)+'.pkl\n')
	
	print 'alignments read', read_count
	print 'alignments created', alns_created
	print 'reads kept', len(reads)+pickle_num*reads_per_pickle
	print 'alignments kept', alns_kept
			
Main()
					
				
		
