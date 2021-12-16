#Samantha Shrum 2021.8.18
#Removes any lines from a .bed12 file that do not have chrX, chrY, or chr# in the first column.

import csv, re

filein = "output/merged.bed12"
fileout = "output/merged_trimmed.bed12"

#open file in
with open(filein, newline = '') as f_in:
	reader = csv.reader(f_in, delimiter = '	')
	#open file out
	with open(fileout, "w") as f_out:
		writer = csv.writer(f_out, delimiter = '\t')
		#go through file in line by line
		for row in reader:
			#regex check for appropriate characters in first column
			if re.search("(chr[0-9xyXY][0-9]?)$", row[0]):
				#if good, write to file out
				writer.writerow(row)
