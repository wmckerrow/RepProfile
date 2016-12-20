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

import os
import re
import sys
import resource
from Bio import SeqIO
from Bio import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import generic_dna
from Bio.Alphabet import IUPAC
import pickle
from sets import Set
import pysam
import argparse

class RepeatMasker(object):
	"""
	RepeatMasker - class for holding UCSC Repeatmasker records
	"""
	def __init__(self, rep_str):
		rep_list = rep_str.split('\t')
		
		self.genoName = rep_list[0]
		self.genoStart = int(rep_list[1])
		self.genoEnd = int(rep_list[2])
		self.repName = rep_list[3]
		self.score = int(rep_list[4])
		self.strand = rep_list[5]

	
def ReadRepeatMaskerFile(repeat_file):
	"""
	ReadRepeatMaskerFile - read the UCSC Repeatmasker file
	returns
	repeats - a dictionary of dictionaries of RepeatMasker objects,
			  the primary key is the chromosome, the secondary is the repeat name
			  the values are a list of RepeatMasker objects
	"""

	repeats = []

	with open(repeat_file, 'r') as f:
		f.readline()  # throw away header
		for line in f:
			repeats.append( RepeatMasker(line.strip()) )
	return repeats

def MakeHoppelBam(repeats,samfile,id,flanking):
	
	read_ids = set([])
	
	# iterate through the sequences, then through the repeats
	# for each Hoppel repeat, insert its sequence in the correct location
	for repeat in repeats:
		if id and not id in repeat.repName:
			continue
		for read in samfile.fetch(repeat.genoName,max(repeat.genoStart-flanking,0),repeat.genoEnd+flanking):
			if not read.qname in read_ids:
				read_ids.add(read.qname)
	
	for read_id in read_ids:
		print read_id
					
def GetArgs():
	"""
	GetArgs - read the command line
	returns - a bam file name, a genome file name, and an output file name
	bam file is required
	genome file is optional and has a default
	output file is option, if not given, output is written to stdout

	typing python phiX.count.py will show the options
	"""

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser = Parser(description='Save only reads aligned to named repeat.')
		parser.add_argument('-r', '--repeat_masker',
							type=str,
							required=True,
							help='Repeat masker file (required).')
		parser.add_argument('-b', '--bam_file',
							type=str,
							required=True,
							help='Bam file (required).')
		parser.add_argument('-i', '--id',
							required=False,
							default=None,
							type=str,
							help='Name of the repeat to collect (all).')
		parser.add_argument('-f', '--flanking',
							required=False,
							type=int,
							default=0,
							help='Include reads this far away from repeats (0).')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.repeat_masker, args.bam_file, args.id, args.flanking


def Main():
	"""
	main routine - get the input and output files and print the output
	"""

	repeat_masker_file, bam_file, id, flanking = GetArgs()
	
	repeats = ReadRepeatMaskerFile(repeat_masker_file)
	samfile =  pysam.Samfile(bam_file,'rb')
	
	MakeHoppelBam(repeats,samfile,id,flanking)
	
Main()
