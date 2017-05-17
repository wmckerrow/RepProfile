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

Print RepProfile results.

"""

import pickle
import sys
import argparse
import numpy

# Get command line arguments
def GetArgs():

	def ParseArgs(parser):
		class Parser(argparse.ArgumentParser):
			def error(self, message):
				sys.stderr.write('error: %s\n' % message)
				self.print_help()
				sys.exit(2)
		parser.add_argument('-r', '--report_rep_types',
							type=str,
							required=False,
							default=None,
							help='List of repeat (hyper) types to report - comma separated. (None)')
		parser.add_argument('-p', '--report_pos_types',
							type=str,
							required=False,
							default=None,
							help='List of position types to report - comma separated. (None)')
		parser.add_argument('-n', '--report_base_prob_at_pos_types',
							type=str,
							required=False,
							default=None,
							help='Report fraction of these bases at corresponding position types. Separate sets of bases to report for each reported position type by comma. For example -p edit_f,edit_r -n G,C would report that fraction G at forward edited sites and the fraction C at reversed edited sites. (None)')
		parser.add_argument('-R', '--rep_type_pickle',
							type=str,
							required=False,
							default='rep_type.pkl',
							help='Repeat (hyper) types saved as pickle by RepProfile. (rep_type.pkl)')
		parser.add_argument('-P', '--pos_type_pickle',
							type=str,
							required=False,
							default='pos_type.pkl',
							help='Position types saved as pickle by RepProfile. (pos_type.pkl)')
		parser.add_argument('-G', '--genome_profile_pickles',
							type=str,
							required=False,
							default='genome_profile_f.pkl,genome_profile_r.pkl',
							help='Genome profiles pickle files. Forward and reverse comma separated. (genome_profile_f.pkl,genome_profile_r.pkl)')
		parser.add_argument('-g', '--genomic_positions',
							type=str,
							required=False,
							default=None,
							help='If the repeat fasta built by make_repeat_genome.py is provided, locations will be reported in genomic coordinates (1-based). (None)')
		parser.add_argument('-f', '--flanking',
							type=int,
							required=False,
							default=1000,
							help='Amount of flanking sequence used. (1000)')
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.pos_type_pickle, args.rep_type_pickle, args.genome_profile_pickles, args.report_rep_types, args.report_pos_types, args.report_base_prob_at_pos_types, args.flanking, args.genomic_positions
	
def main():
	pos_type_pickle, rep_type_pickle, genome_profile_pickles, report_rep_types, report_pos_types,report_base_prob_at_pos_types, flanking, genomic_positions = GetArgs()
	
	nuc2num = {'A':0, 'C':1, 'G':2, 'T':3,'a':0, 'c':1, 'g':2, 't':3}
		
	genome_profile_f = pickle.load(open(genome_profile_pickles.split(',')[0],'rb'))
	genome_profile_r = pickle.load(open(genome_profile_pickles.split(',')[1],'rb'))
	
	# Convert from repeat genome coordinates to reference coordinates
	if genomic_positions:
		genomic_positions_dict = dict()
		for line in open(genomic_positions,'r'):
			if line[0] == '>':
				name = line.strip().split(' ')[0][1:]
				region = line.strip().split(' ')[1].split(':')
				genoChrom = region[0]
				if region[1][0] == '-':
					genoStart = -int(region[1][1:].split('-')[0])
					genoEnd = -int(region[1][1:].split('-')[1])
				else:
					genoStart = int(region[1].split('-')[0])
					genoEnd = int(region[1].split('-')[1])
				genomic_positions_dict[name] = (genoChrom,genoStart,genoEnd)
	
	# Print repeats that are in specified repeat types
	if report_rep_types:
		rep_types = pickle.load(open(rep_type_pickle,'rb'))
		for seq in rep_types:
			if rep_types[seq] in report_rep_types.split(','):
				if genomic_positions:
					print seq, genomic_positions_dict[seq][0]+':'+str(genomic_positions_dict[seq][1]+1)+'-'+str(genomic_positions_dict[seq][2]+1), rep_types[seq]
				else:
					print seq, rep_types[seq]
	
	# Print positions that are in specified position types, along with specified profile values at those positions.
	if report_pos_types:
		pos_types = pickle.load(open(pos_type_pickle,'rb'))
		for seq in pos_types:
			if pos_types[seq] is None:
				continue
			for i in range(len(pos_types[seq])):
				for j in range(len(report_pos_types.split(','))):
					if pos_types[seq][i] == report_pos_types.split(',')[j]:
						if report_base_prob_at_pos_types:
							if genomic_positions:
								print genomic_positions_dict[seq][0], genomic_positions_dict[seq][1]+1+i+flanking, pos_types[seq][i], ' '.join([str(genome_profile_f[seq][i+flanking][nuc2num[x]]) for x in report_base_prob_at_pos_types.split(',')[j]])
							else:
								print seq, i+flanking, pos_types[seq][i], ' '.join([str(genome_profile_f[seq][i+flanking][nuc2num[x]]) for x in report_base_prob_at_pos_types.split(',')[j]])
						else:
							if genomic_positions:
								print genomic_positions_dict[seq][0], genomic_positions_dict[seq][1]+1+i+flanking, pos_types[seq][i]
							else:
								print seq, i+flanking, pos_types[seq][i]


if __name__ == '__main__':
	main()
