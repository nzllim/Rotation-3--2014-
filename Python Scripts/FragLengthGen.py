# Python Script
# Generation of MetaFile containing gene lengths according to rd.bed
# Written by Nathaniel (March 2014)

import csv

GeneBed = 'Data\\rd.bed'
Mtable = []

# Opens DataFile
with open(GeneBed) as geneflag:
	for Everyline in geneflag:
		Mtable.append(Everyline.split(sep = '\t'))

Ctable = [['Gene', 'Length']]
tmptable = []
genelength = 0

# Parsing Gene Length for feeding into R
for EverySample in Mtable:
	tmptable = []
	tmptable.append(EverySample[3])
	genelength = 0
	
	if EverySample[5] == '+':
		genelength = int(EverySample[2]) - int(EverySample[1])
	elif EverySample[5] == '-':
		genelength = int(EverySample[2])- int(EverySample[1])
	
	tmptable.append(genelength)
	Ctable.append(tmptable)

with open('FragmentLengths.txt', mode = 'w', newline = '\n') as outlog:
	writer = csv.writer(outlog, delimiter = '\t', lineterminator = '\n')
	writer.writerows(Ctable)