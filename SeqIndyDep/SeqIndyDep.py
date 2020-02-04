methylation_motif_and_position_input = input('What is the file name with methylation motif and position?    ')
genome_input = input('What is the name of the genome file?  ')
gene_coding_input = input('What is the name of the gene annotation table from NCBI?'  )
out_input = input('What is the name of the output file?  ')

#open file with methylation motif & postion 
f = open(methylation_motif_and_position_input,'r')
infile = f.readlines()
f.close()

#make each row a list
methylation_positions = []

for x in infile: 
	methylation_positions.append(int(x.strip()))

#put each methylation position into a single position

methylation_positions.sort()


#make contigs based on methylation positions 
start_pos = []
end_pos = []

for x in methylation_positions: 
	start_pos.append(x) 
	end_pos.append(x) 

#pull in fasta genome file 
g = open(genome_input,'r')
genome_in = g.readlines() 
g.close() 

dict_genome = dict() 
key = ""

for x in genome_in:
	x = x.strip('\n')
	if x.startswith('>'):
		key = x
		dict_genome[key] = ""
	else:
		dict_genome[key] += x
#get length of the genome 
for x in dict_genome:
	y = dict_genome[x]
	len_genome = len(y)

#Add in the start and end of the genome
start_pos.insert(0,0)
end_pos.append(len_genome) 

coordinates_contigs = zip(start_pos,end_pos)

#add number to contigs pairing to keep track 

count = 0
contigs = []
for x in coordinates_contigs:
	listx = list(x) 	
	listx.insert(0, count) 
	count = count + 1
	contigs.append(listx)

#Determine which contigs have a start or end that fall
#within a coding sequence 

#import coding sequence 
h = open(gene_coding_input,'r') 
infile_gene = h.readlines()
h.close() 

#make rows as lists 
per_row_gene = []

for x in infile_gene:
	per_row_gene.append(x.strip().split('\t'))

#delete headers
del per_row_gene[0]

#get only the information needed 
gene_calls = []

for x in per_row_gene:
	z = []	
	alias = x[6]
	start = x[12]
	end = x[13]
	z.append(alias)
	z.append(start)
	z.append(end)
	gene_calls.append(z)

out = []

#Determine whether the CDS start, end, surround, or are surrounded by the contigs generated
for x in contigs:
	z = []	
	start_contig = int(x[1])
	end_contig = int(x[2])
	z.append(x[0])
	z.append(start_contig)
	z.append(end_contig)
	len_contig = end_contig - start_contig
	z.append(len_contig)
	ee = []	
	aa = ['Starts in Contig, Ends Outside']
	bb = ['Starts Outside, Ends in Contig']
	cc = ['Gene Surrounds Contig']
	dd = ['Contig Surrounds Gene']
	for y in gene_calls:
		start_gene = int(y[1])
		end_gene = int(y[2])

		#Does Gene Start in/at Contig, but end elsewhere?
		if start_gene >= start_contig and start_gene <= end_contig and end_gene >= end_contig:
			aa.append(y[0])
		#Does Gene End in/at contig Contig but start elsewhere?
		if start_gene <= start_contig and end_gene >= start_contig and end_gene <= end_contig:
			bb.append(y[0])
		#Does Gene Surround Contig?
		if start_gene < start_contig and end_gene > end_contig:
			cc.append(y[0])
		#Does Contig Surround Gene?
		if start_gene >= start_contig and end_gene <= end_contig:
			dd.append(y[0])
	ee.append(aa)
	ee.append(bb)
	ee.append(cc)
	ee.append(dd) 
	z.append(ee)
	out.append(z) 

#Write to file 
final = []

headers = ['Contig','Start Pos','End Pos','Length','Starts in Contig, Ends Outside','Starts in, Ends Out Count','Starts Outside, Ends in Contig','Starts Out, Ends in Count','Gene Surrounds Contig','Gene Surrounds Count','Contig Surrounds Gene','Contig Surrounds Count','Number of CDS Impacted']

for x in out:
	final_per_contig = []
	final_per_contig.append(x[0])
	final_per_contig.append(x[1])
	final_per_contig.append(x[2])
	final_per_contig.append(x[3])
	CDS_methylated_by_motif = []
	for y in x[4]:
		del y[0] 
		final_per_contig.append(y)
		final_per_contig.append(len(y))
		CDS_methylated_by_motif.append(len(y))
	end = 0	
	for x in CDS_methylated_by_motif:
		end = end + x
	final_per_contig.append(str(end))
	final.append(final_per_contig)

final.insert(0,headers)

import csv

with open(out_input,'w') as i:
	writer = csv.writer(i,delimiter='\t')
	writer.writerows(final)

quit() 





