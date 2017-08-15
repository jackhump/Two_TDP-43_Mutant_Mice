'''
Intron length calculator
Finds the average length of intronic sequence for each gene in a gtf
Each transcript's intronic length is found by subtracting the total exon length from the transcript length
'''

import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("infile", type = str )
parser.add_argument("outfile", type = str )
args = parser.parse_args()
print(args.infile)
print(args.outfile)

infile = args.infile
outfile = args.outfile

#infile = "/Users/Jack/Documents/Misc/sample.gtf"
#outfile = "/Users/Jack/Documents/Misc/out.txt"

ofile =  open(outfile, 'w')

def gtf_parser(string, pattern):
	'''
	Returns value from key-value pair stored in GTF column 9
	pattern refers to the key
	'''
	#print string
	#print pattern
	colnine = string[8].split("; ")
	#print "colnine is"
	#print colnine
	for i in colnine:
    		if pattern in i:
			phrase = i
	word = phrase.split(" ")[1].strip('"')
	#print "%s is %s " % ( pattern, word )
	return(word)

with open(infile) as f:
    for line in f:
      if "chr" in line: # to ignore comments at start	
        l = line.split("\t")
	if l[2] == "gene":  # capture gene variable
		if 'gene_name' in locals():
			#print gene_dic
			# work out mean intron length from the dictionary
			i_lengths = []
			for t in gene_dic.keys():
        			if gene_dic[t][2] == "protein_coding":
					i_lengths.append(gene_dic[t][1] - gene_dic[t][0])
			#print i_lengths
			mean_i_length = 0
			if len(i_lengths) > 0:
				mean_i_length = int( np.mean(i_lengths) )
			# write out gene name and intron length to file
			outline = "\t".join( [gene_name, gene_id, gene_type, str(mean_i_length),"\n" ] )
			#print "mean intron length for %s is %i" % ( gene_name, mean_i_length )
			ofile.write(outline) 
		vars = l[8].split("; ")
		gene_name = gtf_parser(string = l, pattern = "gene_name")
		gene_id = gtf_parser(l, "gene_id")
		gene_type = gtf_parser(l, "gene_type")
		#print gene_name
		t_count = 0 # t for transcript
		
	if l[2] == "transcript":
		t_status = gtf_parser(l, "transcript_type" )
		t_name = gtf_parser(l, "transcript_id" )
		#basic = gtf_parser(l,"basic" )	
		#print basic
		gene_dic = {} # initialise dictionary
		gene_dic[t_name] = [0,0,0]  # exon length, transcript length,t_status
		gene_dic[ t_name ][1] = int(l[4]) - int(l[3]) # transcript length
		gene_dic[ t_name ][2] = t_status
		
	if l[2] == "exon":
		#if t_status == "protein_coding":
			e_length = int(l[4]) - int(l[3])
			gene_dic[ t_name ][0]  += int(e_length)

# final gene in the file
i_lengths = []
#print gene_dic
for t in gene_dic.keys():
	i_lengths.append(gene_dic[t][1] - gene_dic[t][0])  # subtract total exon length from gene length
mean_i_length = int( np.mean(i_lengths) )
outline = "\t".join( [gene_name, gene_id, gene_type, str( mean_i_length), "\n" ] )
print "mean intron length for %s is %i" % ( gene_name, mean_i_length )
ofile.write(outline) 

ofile.close()
