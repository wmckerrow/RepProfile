import pickle
import sys
import argparse
import numpy

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
		return parser.parse_args()

	parser = argparse.ArgumentParser()
	args = ParseArgs(parser)

	return args.pos_type_pickle, args.rep_type_pickle, args.report_rep_types, args.report_pos_types
	
def main():
	pos_type_pickle, rep_type_pickle, report_rep_types, report_pos_types = GetArgs()
	
	print pos_type_pickle, rep_type_pickle, report_rep_types, report_pos_types
	
	if report_rep_types:
		rep_types = pickle.load(open(rep_type_pickle,'rb'))
		for seq in rep_types:
			if rep_types[seq] in report_rep_types.split(','):
				print seq, rep_types[seq]
				
	if report_pos_types:
		pos_types = pickle.load(open(pos_type_pickle,'rb'))
		for seq in pos_types:
			if pos_types[seq] == None:
				continue
			for i in range(len(pos_types[seq])):
				if pos_types[seq][i] in report_pos_types.split(','):
					print seq, i, pos_types[seq][i]


if __name__ == '__main__':
	main()