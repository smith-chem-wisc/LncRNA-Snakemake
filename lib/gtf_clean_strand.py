#Samantha Shrum 2021.8.20
#Removes any line from a gtf which has a strand value of '.'

import csv, re

filein = "output/merged.gtf"
fileout = "output/merged_strand_clean.gtf"

#open file in
with open(filein, newline = '') as f_in:
	reader = csv.reader(f_in, delimiter = '\t')
	#open file out
	with open(fileout, "w") as f_out:
		writer = csv.writer(f_out, delimiter = '\t', quotechar = '', quoting = csv.QUOTE_NONE, escapechar = '\\')
		#go through file in line by line
		for row in reader:
			if re.match("#", row[0]):
				writer.writerow(row)
				continue
			#regex check for appropriate characters in first column
			if re.search("[+-]$", row[6]):
				#if good, write to file out
				writer.writerow(row)