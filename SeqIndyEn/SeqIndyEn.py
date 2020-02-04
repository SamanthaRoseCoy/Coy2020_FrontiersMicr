def chunkstring(string, length, motif_size):
	#get the coordinates of each fragment with the window sliding 
	#based on the length of the motif, and the length of the fragment
	#is determined by a user input
	start_pos = 0
	coordinates_of_all_windows = []
	while start_pos < len(string):
		coordinate_of_one_window = [start_pos,start_pos+length]
		coordinates_of_all_windows.append(coordinate_of_one_window)
		start_pos = start_pos + motif_size
	#pull out each fragment of the genome 
	genome_fragments = []
	
	for coordinates in coordinates_of_all_windows:
		start = coordinates[0]
		end = coordinates[1]
		if end > len(string):
			genome_fragments.append(string[start:len(string)])
		else: 
			genome_fragments.append(string[start:end])

	return (genome_fragments)

#This function gets the fragment coordinates, and compiles
#a list of list with an assigned fragment, contig number, and coordinates
def coordinate_finder(genome_fragments,genome):
	fragment_coordinates_number = []
	#count assigns a number to the genome fragment 
	count = 0
	for fragment in genome_fragments: 
		count = count + 1
		#get the starting position of the string using .find, and get the 
		#rest by the addition of the length of the fragment
		fragment_coordinate = []
		fragment_coordinate.append(genome.find(fragment))
		fragment_coordinate.append(genome.find(fragment) + len(fragment))
		#assemble list of count, coordiantes, and fragment 
		by_fragment = [fragment,genome.find(fragment),fragment_coordinate]
		fragment_coordinates_number.append(by_fragment)

	return(fragment_coordinates_number)

#this function finds the motifs in the order in the motif list 
#and returns the positions of those motifs in the genome 
def motif_positions_of_fragment(motifs_to_search,fragment_list):

	out_list_fragments = []
		
	for fragment in fragment_list:
		motif_coordinates = []
		motif_count = []
		out_list_fragment = [fragment[0],fragment[1],fragment[2]]
		for motif in motifs_to_search:
				
				by_motif_position = []
				by_genome_position = []
				#find the start and end of each position in the fragment
				p = re.compile(motif)
				for m in p.finditer(fragment[0]):
					motif_coordinate = [m.start(),m.end()]
					by_motif_position.append(motif_coordinate)
				#find the start and end of each position in the genome 
				for motif_position in by_motif_position:
					#add to the position the start of the fragment in the genome
					motif_start_genome = motif_position[0] + fragment[2][0]
					motif_end_genome = motif_position[1] + fragment[2][0]
					motif_coordinate_genome = [motif_start_genome,motif_end_genome]
					by_genome_position.append(motif_coordinate_genome)
				#write the position of each to a file 
				if len(by_genome_position) != 0:
					motif_coordinates.append(by_genome_position)
				#record the number of that one type of motif on the 
				out_list_fragment.append(len(by_genome_position))
				out_list_fragment.append(by_genome_position)
		
		#pool together all motifs for total per fragment 
		motif_coordinates_fin = []

		for motif_coordinate_groups in motif_coordinates:
			for coordinate in motif_coordinate_groups:
				motif_coordinates_fin.append(coordinate)

		out_list_fragment.insert(3,len(motif_coordinates_fin))
		out_list_fragment.insert(4,motif_coordinates_fin)

		out_list_fragments.append(out_list_fragment)
	
	#insert headers 
	headers = ['genome fragment','fragment starting position','fragment position','total count','all coordinates']

	for motif in motifs_to_search:
		motif_count = motif + ' count'
		motif_positions = motif + ' coordinates'
		headers.append(motif_count)
		headers.append(motif_positions)


	out_list_fragments.insert(0, headers)


	return out_list_fragments

#this function removes the redundancies 
def removal_of_redundancies(list_of_list_with_redundancies,genome, motifs_to_search):
	#get all of the sets of coordinates 
	del list_of_list_with_redundancies[0]
	all_coordinates = []
	for line in list_of_list_with_redundancies:
		if len(line[4]) != 0:
			if line[4] not in all_coordinates:
				all_coordinates.append(line[4])
	
	#call new out list
	out_list = []

	for coordinates in all_coordinates:
		#get lines with those exact coordinates
		lines_with_those_coordinates = []
		for line in list_of_list_with_redundancies:
			if coordinates == line[4]:
				lines_with_those_coordinates.append(line)
		#get the starting position and end position of each line 
		starting_pos = []
		ending_pos = []
		for line in lines_with_those_coordinates:
			starting_pos.append(line[2][0])
			ending_pos.append(line[2][1])
		starting_pos.sort()
		ending_pos.sort()


		data_to_add = lines_with_those_coordinates[0][3:]
		#generate new lines 
		new_line = [genome[starting_pos[0]:ending_pos[-1]],starting_pos[0],ending_pos[-1]]
		
		for data in data_to_add:
			new_line.append(data)

		out_list.append(new_line)
		


	#Add in the headers 
	headers = ['genome fragment','fragment starting position','fragment end position','total count','all coordinates']

	for motif in motifs_to_search:
		motif_count = motif + ' count'
		motif_positions = motif + ' coordinates'
		headers.append(motif_count)
		headers.append(motif_positions)


	out_list.insert(0, headers)
	return(out_list)


#import 
import re
import csv 

#input files 
input_genome = input('What is the name of the input genome fasta file?  ')
input_chunk_size = input('What size fragment do you want the genome broken into?   ')
input_motif = input('What is the input motif flie? (one motif per line)   ')
input_motif_size = input('What is the length of the motif you are using? (If motifs are different length one must be picked to ensure the same starting position of each fragment for analysis)  ')
input_out = input('What is the name of the out file? (save as a .txt)  ')

#open genome fasta file, strip header, make genome a single string
f = open(input_genome,'r')
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

#open motif file and make a list where each term is the motif
g = open(input_motif, 'r')
motif_infile = g.readlines()
g.close() 

motif_list = []

for line in motif_infile:
	motif_list.append(line.strip())

#divide the genome into multiple chunks depending on the size wanted 
genome_window = chunkstring(genome, int(input_chunk_size), int(input_motif_size))
#get coordinates for each of these genome fragments 
fragment_with_coordinates = coordinate_finder(genome_window, genome)
#find motifs in each fragment and get coordinates of each motif in genome
with_all_redundancies = motif_positions_of_fragment(motif_list,fragment_with_coordinates)

#write to an out file 
with_redundancies_out_file_name = input_out.replace('.txt','') + '_with_redundancies.txt'

with open(with_redundancies_out_file_name,'w') as o1:
	writer = csv.writer(o1,delimiter='\t')
	writer.writerows(with_all_redundancies)


#remvoe redundancies
out = removal_of_redundancies(with_all_redundancies,genome, motif_list)

#write to an out file 
with open(input_out,'w') as o2:
	writer = csv.writer(o2,delimiter='\t')
	writer.writerows(out)

quit()