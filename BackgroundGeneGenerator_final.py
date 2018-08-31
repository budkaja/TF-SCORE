import sys
import math
import os
import string
from operator import itemgetter
import numpy as np
from scipy import stats

def TSS_TF_parser(TSSFileName, in_bed_file, upstream, downstream, output_file_name, fasta):

        TSSfile = open(TSSFileName, "r")
        BindingSiteDict = {}

        GeneList = []
        InputGeneCount = 0
        MatchGeneCount = 0

	DNAfile = output_file_name
	out = open(DNAfile, "w")
	out.write('Gene' + '\t' + 'total open DNA' + '\t' + 'Open Regions' + '\n')

        infile = open(in_bed_file, "r")

	BedDict = {'chr1':[], 'chr2':[], 'chr3':[], 'chr4':[], 'chr5':[], 'chr6':[], 'chr7':[], 'chr8':[], 'chr9':[], 'chr10':[], 'chr11':[], 'chr12':[], 'chr13':[], 'chr14':[], 'chr15':[], 'chr16':[], 'chr17':[], 'chr18':[], 'chr19':[], 'chr20':[], 'chr21':[], 'chr22':[], 'chrX':[], 'chrY':[]}

	if fasta == "Yes":
		totalsites = 0
	        for lines in infile:
        	        line = lines.strip()
                	if(line[0].startswith('>')):
                        	information = line[1::]
                        	chromosome, start, stop = information.replace(':', ' ').replace('-', ' ').split()
				if chromosome not in BedDict.keys():
					BedDict.update({chromosome:[]})
				BedDict[chromosome].append((start, stop))
				totalsites = totalsites + 1
		print "Total sites: ", totalsites
	else:
		totalsites = 0
		for lines in infile:
			line = lines.strip()
			BindingSite = line.split()
			chromosome = BindingSite[0]
			start = BindingSite[1]
			stop = BindingSite[2]
			if chromosome not in BedDict.keys():
				BedDict.update({chromosome:[]})
			BedDict[chromosome].append((start, stop))
			totalsites = totalsites + 1
		print "Total sites: ", totalsites
	
	for chrom in BedDict.keys():
		print chrom, len(BedDict[chrom])
	for openRegionChromosome, sitelist in BedDict.iteritems():
		sortedlist = sorted(sitelist, key = lambda x: x[1])
		BedDict[openRegionChromosome] = sortedlist

	inlist = 0
	outlist = 0

        for TSSline in TSSfile:
                Total_DNA_Count = 0
		sitelist = ''
		remove = TSSline.strip()
                TSS = remove.split()
                TSS_gene = TSS[0]
                TSS_chromosome = TSS[1]
                TSS_location = int(TSS[2])
                TSS_strand = TSS[3]
                if TSS_strand == '+':
                        if TSS_location > upstream:
                                RegulatoryRange = [TSS_chromosome, TSS_location - upstream, TSS_location + downstream]
                        else:
                                RegulatoryRange = [TSS_chromosome, 0, TSS_location + downstream]

                else:
                        if TSS_location > downstream:
                                RegulatoryRange = [TSS_chromosome, TSS_location - downstream, TSS_location + upstream]
                        else:
                                RegulatoryRange = [TSS_chromosome, 0, TSS_location + upstream]

		if TSS_chromosome not in BedDict.keys():
			continue
		else:
			for site in BedDict[TSS_chromosome]:
				a = int(RegulatoryRange[1])

                        	b = int(RegulatoryRange[2])
	                        c = int(site[0])
        	                d = int(site[1])
                	        distance = c - TSS_location
                        	motifsize = d - c
                        	if a <= c and c <= b:
                                	Total_DNA_Count = Total_DNA_Count + motifsize
					if sitelist == '':
						sitelist = str((TSS_chromosome, c, d))
					else:
						sitelist = sitelist + '//' + str((TSS_chromosome, c, d))
				if RegulatoryRange[2] > site[1]:
					break
		
						
		if Total_DNA_Count > 0:
			outputString = "%s\t%d\t%s\n" % (TSS_gene, Total_DNA_Count, sitelist)
			out.write(outputString)
			inlist = inlist + 1
		else:
			outlist = outlist + 1
			continue
	print "In: ", inlist
	print "Out: ", outlist

def ProcessCLI(args):
        print args
	in_bed_file = "/N/u/jubudka/Mason/EncodeMergedBed.bed"
        in_filename = "/N/u/jubudka/Mason/hg19_TSS_list.txt"

	upstream = 20000
	downstream = 10000
	fasta = "No"
	string = "unnamed"
	for i in xrange(len(args)):
		if args[i] == "-u":
			upstream = int(args[i+1])
		elif args[i] == "-d":
			downstream = int(args[i+1])
		elif args[i] == "-i":
			in_bed_file = str(args[i+1])
		elif args[i] == "-f":
			fasta = "Yes"	
		elif args[i] == "-o":
			string = str(args[i+1])
		elif args[i] == "-a":
			in_filename = str(args[i+1])
	
	if string == "unnamed":
		string = "BackgroundDNA_up_" + str(upstream) + "_do_" + str(downstream) + ".txt"

	output_file_name = string

	TSS_TF_parser(in_filename, in_bed_file, upstream, downstream, output_file_name, fasta)

	print "finished BindingInformation"


if __name__ == '__main__':
     ProcessCLI(sys.argv)


