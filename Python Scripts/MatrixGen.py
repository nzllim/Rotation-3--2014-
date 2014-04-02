# Python Script
# Generation of MetaFile Design Matrix for feeding into EdgeR or limma pipelines (R).
# Output Matrix File has slightly formatting errors in the last line. Manual formatting necessary.
# Written by Nathaniel (March 2014)

import csv

FragmentFile = 'Data\\fragment_counts.txt'
Mtable = []

# Opens DataFile
with open(FragmentFile) as fragmentflag:
	for Everyline in fragmentflag:
		Mtable.append(Everyline.split(sep = '\t'))

HeaderList = Mtable[0]

Ctable = [['Sample', 'Type', 'Set', 'Timepoint']]
tmptable = []
breakuptable = []

# Recreating the HeaderList into a MetaFile Matrix for feeding into R
for EverySample in HeaderList:
	tmptable = []
	breakuptable = []
	tmptable.append(EverySample)
	breakuptable = str(EverySample).split(sep = '_')
	
	tmptable.extend(breakuptable)
	Ctable.append(tmptable)

with open('FragmentMatrix.txt', mode = 'w', newline = '\n') as outlog:
	writer = csv.writer(outlog, delimiter = '\t', lineterminator = '\n')
	writer.writerows(Ctable)