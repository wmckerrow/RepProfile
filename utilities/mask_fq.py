import sys

file = open(sys.argv[1],'r')
mask_from = sys.argv[2]
mask_to = sys.argv[3]

while True:
	line = file.readline().strip()
	if not line:
		break
	print line
	seq = file.readline().strip()
	newseq = ''
	for x in seq:
		if x == mask_from:
			newseq += mask_to
		else:
			newseq += x
	print newseq
	print file.readline().strip()
	print file.readline().strip()