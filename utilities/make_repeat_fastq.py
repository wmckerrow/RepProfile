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

import gzip
import argparse
import sys
from sets import Set

def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)

		parser.add_argument('-l', '--lists_of_read_ids',
							type=str,
							required=True,
							help='Text files listing the read ids you want - comma separated.')
		parser.add_argument('-f', '--fastq',
							required=False,
							type=str,
							default=None,
							help='fastq file not gzipped. Specify exactly one of fastq,gzipped_fastq.')
		parser.add_argument('-g', '--gzipped_fastq',
							required=False,
							type=str,
							default=None,
							help='fastq file gzipped. Specify exactly one of fastq,gzipped_fastq.')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.lists_of_read_ids, args.fastq, args.gzipped_fastq
	
def Main():
	lists_of_read_ids, fastq, gzipped_fastq = GetArgs()
	
	if fastq and gzipped_fastq:
		print 'Specify only one type of read file.'
		exit()
	if not fastq and not gzipped_fastq:
		print 'Specify a read file with -f or -g. Use -h/--help for more information.'
		exit()
	
	read_ids = Set([])
	
	for list_of_read_ids in lists_of_read_ids.split(','):
		for line in open(list_of_read_ids,'r'):
			read_ids.add(line.strip())
		
	#print read_ids
	
	if fastq:
		with open(fastq, 'r') as f:
			line_num = 0
			print_line = False
			for line in f:
				if line.strip().split(' ')[0].split('/')[0][1:] in read_ids:
					print_line = True
			
				if print_line:
					print line.strip()
				
				line_num = (line_num+1)%4
				if line_num == 0:
					print_line = False
	
	if gzipped_fastq:
		with gzip.open(gzipped_fastq, 'rb') as f:
			line_num = 0
			print_line = False
			for line in f:
				if line.strip().split(' ')[0].split('/')[0][1:] in read_ids:
					print_line = True
			
				if print_line:
					print line.strip()
				
				line_num = (line_num+1)%4
				if line_num == 0:
					print_line = False
	
Main()
