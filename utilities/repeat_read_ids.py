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

Print unique read ids for reads that overlap a given repeat

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

# iterate through repeats for each repeat, fetch any overlapping reads and save ids in a set object
def MakeHoppelBam(repeats,samfile,id,flanking):
	
	read_ids = Set([])
	
	for repeat in repeats:
		if id and id.upper() != repeat.repName[:len(id)].upper():
			continue
		for read in samfile.fetch(repeat.genoName,max(repeat.genoStart-flanking,0),repeat.genoEnd+flanking):
			if not read.qname in read_ids:
				read_ids.add(read.qname)
	
	# Print read ids
	for read_id in read_ids:
		print read_id
			
# Read command line arguments		
def GetArgs():

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
