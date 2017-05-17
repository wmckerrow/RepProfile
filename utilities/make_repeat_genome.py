"""

Copyright Wilson McKerrow, William Thompson 2017

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

Extract sequences tagged as a particular repeat along with flanking sequence.

"""

import sys
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
import argparse

class RepeatMasker(object):
	"""
	RepeatMasker - class for holding Repeatmasker data as read from a bed file
	"""
	def __init__(self, rep_str):
		rep_list = rep_str.split('\t')
		
		self.genoName = rep_list[0]
		self.genoStart = int(rep_list[1])
		self.genoEnd = int(rep_list[2])
		self.repName = rep_list[3]
		self.score = int(rep_list[4])
		self.strand = rep_list[5]


def ReadGenomeFile(fasta_file):
	"""
	ReadGenomeFasta - read a genome file in fasta format
	returns
	a dictionary of SeqRecored objects, key is teh fasta id
	"""

	seq_recs = dict()
	f = open(fasta_file, 'r')
	for rec in SeqIO.parse(f, 'fasta'):
		seq_recs[rec.id] = rec

	return seq_recs
	
def ReadRepeatMaskerFile(repeat_file):
	"""
	ReadRepeatMaskerFile - read the UCSC Repeatmasker file
	returns
	repeats - a dictionary of dictionaries of RepeatMasker objects,
			  the primary key is the chromosome, the secondary is the repeat name
			  the values are a list of RepeatMasker objects
	"""

	repeats = list()

	with open(repeat_file, 'r') as f:
		for line in f:
			repeats.append( RepeatMasker(line.strip()) )
	return repeats
	
# Read command line arguments
def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-r', '--repeat_name',
							type=str,
							required=False,
							default=None,
							help='Only build a genome from repeats that match this name.')
		parser.add_argument('-f', '--flanking',
							required=False,
							default=1000,
							type=int,
							help='Amount of sequence flanking repeats. (1000)')
		parser.add_argument('-m', '--repeat_masker_bed',
							required=True,
							type=str,
							help='A bed file listing repeats to build the genome from.')
		parser.add_argument('-g', '--genome_fasta',
							required=True,
							type=str,
							help='The full genome in fasta format. The repeat genome will be cut out of the full genome.')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.repeat_name, args.flanking, args.repeat_masker_bed, args.genome_fasta
	
def Main():
	repeat_name, flanking, repeat_masker_bed, genome_fasta = GetArgs()

	seq_recs = ReadGenomeFile(genome_fasta)
	
	repeats = ReadRepeatMaskerFile(repeat_masker_bed)
	
	# If a repeat overlaps a previously stored repeat, merge them together. Otherwise record the interval the repeat + flanking sequence covers.
	intervals = dict()
	for repeat in repeats:
		if repeat_name and repeat_name.upper() != repeat.repName[:len(repeat_name)].upper():
			continue
		if not repeat.genoName in intervals:
			intervals[repeat.genoName] = [(repeat.genoStart-flanking,repeat.genoEnd+flanking)]
			continue
		overlaps = []
		new_start = repeat.genoStart-flanking
		new_end = repeat.genoEnd+flanking
		for i in range(len(intervals[repeat.genoName])):
			if new_start < intervals[repeat.genoName][i][1] and new_end > intervals[repeat.genoName][i][0]:
				overlaps.append(i)
				new_start = min(new_start,intervals[repeat.genoName][i][0])
				new_end = max(new_end,intervals[repeat.genoName][i][1])
		for i in overlaps[::-1]:
			intervals[repeat.genoName].pop(i)
		intervals[repeat.genoName].append((new_start,new_end))
	
	# Sort intervals
	sorted_intervals = dict()
	for seq in intervals:
		sorted_intervals[seq] = sorted(intervals[seq], key=lambda pos: pos[0])
	
	# Print sequence in the recorded intervals
	count = 0
	for seq in sorted_intervals:
		for interval in sorted_intervals[seq]:
			count+=1
			if interval[0] < 0:
				start = 0
				added_up_flank = 'N'*(-interval[0])
			else:
				start = interval[0]
				added_up_flank = ''
			if interval[1] > len(seq_recs[seq].seq):
				end = len(seq_recs[seq].seq)
				added_down_flank = 'N'*(interval[1]-end)
			else:
				end = interval[1]
				added_down_flank = ''
			if repeat_name:
				print '>'+repeat_name+'_'+str(count) + ' ' + seq+':'+str(interval[0])+'-'+str(interval[1])
			else:
				print '>repeat_'+str(count) + ' ' + seq+':'+str(interval[0])+'-'+str(interval[1])
			print added_up_flank + seq_recs[seq][start:end].seq.upper() + added_down_flank
		
Main()
