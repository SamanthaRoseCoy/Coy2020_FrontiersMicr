

#imports
import re

#inputs
genome_input = input('What is the genome file?  ')
motif_input = input('What is the motif you are searching for?  ')
methylated_input = input('What position is methylated? (the first position is 0 not 1!)  ')
out_input = input('What is your output file?  ')

#open genome fasta file, strip header, make genome a single string
f = open(genome_input,'r')
genome_infile = f.readlines() 
f.close() 

current_GenBank_ID = ''
dict_fasta ={}

for line in genome_infile:
	#Clear out leading/trailing whitespace (e.g. '\n')
	line = line.strip() 
	#Check to see if the line is a fasta header
	if line.startswith(">"):
		#Get the GenBank_ID 
		current_GenBank_ID = line
		#Add the GenBank_ID to the dictionary 
		dict_fasta[current_GenBank_ID]=""
	else:
		#Add the seqeunce as values to the GenBank ID
		dict_fasta[current_GenBank_ID] += line

for x in dict_fasta:
	genome = dict_fasta[x]

#using regular expressions find the position of each motif
out = []

p = re.compile(motif_input)
for m in p.finditer(genome):
	methylated_motif = m.start() + int(methylated_input)
	out.append(methylated_motif)


o = open(out_input, 'w')

for x in out:
	o.write(str(x))
	o.write('\n')

quit()