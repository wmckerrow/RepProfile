import sys

file = open(sys.argv[1],'r')
mask_from = sys.argv[2]
mask_to = sys.argv[3]

for line in file:
	if line[0] == '>':
		print line.strip()
	else:
		seq = line.strip()
		newseq = ''
		for x in seq:
			if x == mask_from:
				newseq += mask_to
			else:
				newseq += x
		print newseq