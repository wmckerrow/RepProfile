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

Find and replace one character with another in sequence portion(s) of a fasta file.

"""

import sys

file = open(sys.argv[1],'r')
mask_from = sys.argv[2].upper()
mask_to = sys.argv[3].upper()

for line in file:
	if line[0] == '>':
		print line.strip()
	else:
		seq = line.strip()
		newseq = ''
		for x in seq:
			x = x.upper()
			if x == mask_from:
				newseq += mask_to
			else:
				newseq += x
		print newseq
