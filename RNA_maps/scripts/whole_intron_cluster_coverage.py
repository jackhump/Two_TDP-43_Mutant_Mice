from pybedtools import BedTool
import csv
import os
import argparse
parser = argparse.ArgumentParser()

parser.add_argument("intronBed", help = "A bed file of introns that each include an exon")
parser.add_argument("clusterBed", help = "A bed file of iCLIP clusters")
parser.add_argument("outFile", help = "The output file")
args = parser.parse_args()

print( "Intron bed file is %s" % args.intronBed)
print( "Cluster bed file is %s" % args.clusterBed)
print( "output file is %s" % args.outFile)


# introns = BedTool("/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/F210I_embryonic_brain_se_intron_included.bed")
# clusters = BedTool("/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/iCLIP/All_TDP_iCLIP_merged.bed")
# outFile="/Users/Jack/Documents/Misc/intron_coverage.csv"

outFile = args.outFile
introns = BedTool(args.intronBed)
clusters = BedTool(args.clusterBed)


flank = 100

# first flank the introns by 100nt either side

def flank_introns(feature):
	feature.start = feature.start - flank
	feature.end = feature.end + flank
	return(feature)

if not os.path.exists(os.path.dirname(outFile)):
	os.mkdir(os.path.dirname(outFile))

if not outFile.endswith("csv"):
	outFile = outFile + "_" + str(flank) + "_coverage.csv"

# flank each intron by 100 to capture exons either side
flanked = introns.each(flank_introns)

# intersect the flanked introns with the iCLIP clusters
#intersect = introns.intersect(clusters, s=True,wa=True, wb = True)
intersect = flanked.intersect(clusters, s=True,wa=True, wb = True)

# store as dictionary
intervals = dict()
for feature in intersect:
		coord = str(feature.chrom) + ":" + str(feature.start) + "-" + str(feature.stop)
		length = int(feature.stop) - int(feature.start)
		strand = str(feature.strand)
		start_pos = int(feature.start)
		if coord not in intervals:
			intervals[coord] = [ strand, length, [ int(feature[7]) - start_pos, int(feature[8]) - start_pos, ]]
		if coord in intervals:
			intervals[coord].append( [ int(feature[7]) - start_pos,int(feature[8]) - start_pos, ] )

# for each entry in the dictionary I now have the feature length and a set of tuples that are the positions of the clusters
output_list = []
for feature in intervals:
	strand = intervals[feature][0]
	length = intervals[feature][1]
	vector = [0] * length
	for cluster in intervals[feature][2:]:
		if cluster[0] < 0: # deal with clusters that overlap the edges
			cluster[0] = 0
		if cluster[1] > length:
			cluster[1] = length
		cluster_length = cluster[1] - cluster[0]
		vector[ cluster[0]:cluster[1] ] = [1] * int(cluster_length)
	if strand == "-":
		vector.reverse()
	vector.insert(0, str(feature))
	#vector.append("\n")
	output_list.append(vector)

# write out the intervals to a comma separated file
print( "Writing to %s" % outFile)

with open(outFile, 'w') as f:
	for vector in output_list:
		csv.writer(f).writerow(vector)



# testing
feature="chr6:25794574-25800384"
clusters=BedTool("/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/iCLIP/All_TDP_iCLIP_merged.bed")
introns=BedTool("/Users/Jack/SAN/IoN_RNAseq/RNA_Maps/data/F210I_embryonic_brain_notsig_se_intron_all.bed")




