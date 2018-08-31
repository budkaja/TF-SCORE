import sys
import math
import os
import string
from operator import itemgetter
import numpy as np
import random
import scipy
from scipy import stats
from scipy import special
from sets import Set

def IntersectRegion(siteA, siteB, tuplelist):
	output = 0
	for region in tuplelist:
		if siteA[0] > region[1]:
			continue
		elif siteA[0] <= region[1] and siteA[0] >= region[0]:
			if siteB[0] <= region[1] and siteB[0] >= region[0]:
				output = region[1] - region[0] + 1
				break
			else:
				break
	
	return output

def binarySearch(listOlists, siteStart, siteEnd, chromosome):
	first = 0
        last = len(listOlists) - 1
        while first < last:
                midpoint = (first + last) // 2
                a = int(listOlists[midpoint][0])
                b = int(listOlists[midpoint][1])
                if siteEnd > a:
			first = midpoint + 1
                else:
			last = midpoint	
	return last

def Overlap(filename, TSSDict, upstream, downstream):
	#Overlap the regions in the TSS Start Site Dictionary, within the provided upstream and downstream distances, with each Transcription Factor site file

	#This outputs a dictionary with gene names as the keys and
	# a list of site locations, distances from TSS, motif size, and strand for that gene as the values

	#infile is the Transcription Factor Site file
	infile = open(filename, "r")
	sitelist = []
	TFDict = {}
	sitecounter = 0
	#Prase through the TFsite file and find overlapping regions with the TSS dictionary +/- the upstream and downstream distances
	for line in infile:
		remove = line.strip()
		BindingSite = remove.split('\t')
		if BindingSite[0] not in TSSDict.keys():
			continue
		sites = TSSDict[BindingSite[0]]
		if len(sites) == 0:
			continue
		siteStart = int(BindingSite[1])
		siteEnd = int(BindingSite[2])
		x = binarySearch(sites, siteStart, siteEnd, BindingSite[0])
		while siteStart < int(sites[x][1]) and x >= 0:
			if int(sites[x][0]) <= siteStart and siteEnd <= int(sites[x][1]):
				sitecounter = sitecounter + 1
				gene = sites[x][3]
				TSS_location = int(sites[x][2])
				distance = min((siteStart - TSS_location),(siteEnd - TSS_location))
				motifsize = siteEnd - siteStart + 1
				summary = BindingSite[0] + ':' + str(siteStart) + '-' + str(siteEnd)
				#Add the sites thatoverlap with a TSS regulatory region to the output dictionary 
				#[chr #, start, stop, distance TSS, motif size, and strand]
				if gene in TFDict.keys():
					TFDict[gene].append((BindingSite[0], siteStart, siteEnd, distance, motifsize, BindingSite[4]))
				else:
					TFDict[gene] = [(BindingSite[0], siteStart, siteEnd, distance, motifsize, BindingSite[4])]

			x = x - 1

	
	return TFDict 


def Gene_Cooccurrence(TFa, TFb, BackgroundTFDict, TFDict, DNADict, BackgroundTotalDNACount, GoITotalDNACount, GoIGenes, BackgroundRegionDict):
	#print "Gene Cooccurrence for: ", TFa, TFb
	#initiate sets to hold the regions that have a TF site for A or B
        TFa_regionlist = Set([])
        TFb_regionlist = Set([])
        All_regionlist = Set([])

	#initiates a dictionary that contains each open region and the sites that map to that region
        TFaRegion_SiteDict = {}
        TFbRegion_SiteDict = {}
	
	#initiate counters for TFa and TFb sites in genes of interest set
	TFa_sites_actual = 0
	TFb_sites_actual = 0
	#get the genes with TFa and TFb sites
	TFa_genelist = TFDict[TFa].keys()
	TFb_genelist = TFDict[TFb].keys()
	#initiate list for genes with TFa or TFb sites
	geneID_list = []
	#initiate list to first hold amount of open DNA per gene and then gets converted into a probability of a random site hitting that gene
	prob_list = []
	#counter for the total amount of open DNA in the gene of interest set
	totalDNA = 0

	#beginning checking if the gene has a binding site for TFa and/or TFb and how many binding sites are there
	for gene in GoIGenes:
		gene_regions = BackgroundRegionDict[gene]
		All_regionlist.update(gene_regions)

		if gene in TFa_genelist:
			TFa_sites_actual = TFa_sites_actual + len(TFDict[TFa][gene])
			TFa_sitelist = TFDict[TFa][gene]
			for region in gene_regions:
				region_start = int(region[1])
				region_stop = int(region[2])

				for site in TFa_sitelist:
					site_start = int(site[1])
					site_stop = int(site[2])
					if site_start >= region_start and site_start < region_stop:
						TFa_regionlist.add(region)
						if region not in TFaRegion_SiteDict.keys():
							TFaRegion_SiteDict[region] = [site]
						else:
							TFaRegion_SiteDict[region].append(site)
					else:
						pass
		else:
			pass

		if gene in TFb_genelist:
			TFb_sites_actual = TFb_sites_actual + len(TFDict[TFb][gene])
			TFb_sitelist = TFDict[TFb][gene]
                        for region in gene_regions:
                                region_start = int(region[1])
                                region_stop = int(region[2])

                                for site in TFb_sitelist:
                                        site_start = int(site[1])
                                        site_stop = int(site[2])
                                        if site_start >= region_start and site_start < region_stop:
                                                TFb_regionlist.add(region)
                                                if region not in TFbRegion_SiteDict.keys():
                                                        TFbRegion_SiteDict[region] = [site]
                                                else:
                                                        TFbRegion_SiteDict[region].append(site)
                                        else:
                                                pass
		else:
			pass

		geneID_list.append(gene)
		prob_list.append(DNADict[gene])
		totalDNA = totalDNA + DNADict[gene]

        RandomIntegerOverlaps = []
        TotalRegionCount = len(All_regionlist)

	All_regions = list(All_regionlist)

	#create a region dictionary with numbers as keys for the np.random.choice of regions later. Choice can't use 2d objects
	region_value_counter = 0
	region_value_dictionary = {}
	
	region_problist = []
	for region in All_regions:
		region_start = int(region[1])
		region_stop = int(region[2])
		region_dna = region_stop - region_start
		region_problist.append(region_dna)

		region_value_counter = region_value_counter + 1
		region_value_dictionary[region_value_counter] = region
	for x in xrange(len(region_problist)):
		region_problist[x] = region_problist[x] / float(totalDNA)

        TFaRegionCount = len(TFa_regionlist)
        TFbRegionCount = len(TFb_regionlist)
        OverlapRegionCount = len(TFa_regionlist.intersection(TFb_regionlist))

	#initiate counter for how many TFa and TFb sites are in the background genes (not in the GoI set)
	TFa_sites_background = 0
	TFb_sites_background = 0

	#beginning counting the number of binding sites in the background genes for TFa and TFb 
	for gene, sitelist in BackgroundTFDict[TFa].iteritems():
		TFa_sites_background = TFa_sites_background + len(sitelist)
	for gene, sitelist in BackgroundTFDict[TFb].iteritems():
		TFb_sites_background = TFb_sites_background + len(sitelist)

	#Using the TFa and TFb background site size and the amount of open DNA in GoI and background, determine how many GoI sites are expected for TFa and TFb
	TFa_GoI_expected_sites = int(round((GoITotalDNACount * TFa_sites_background) / float(BackgroundTotalDNACount)))
	TFb_GoI_expected_sites = int(round((GoITotalDNACount * TFb_sites_background) / float(BackgroundTotalDNACount)))

	#convert gene open DNA values in prob_list to probabilities of gene open DNA over total GoI DNA	
	for x in xrange(len(prob_list)):
		prob_list[x] = prob_list[x] / float(totalDNA)

	#initiate random lists that represent the random shuffling of TFa and TFb sites within the open DNA for the GoI
	random_listA = []
	random_listB = []

	random_region_listA = []
	random_region_listB = []

	#initiate random lists that represent the random placement of sites using the background expected site counts for TFa and TFb
	background_expected_listA = []
	background_expected_listB = []

	background_expected_region_listA = []
	background_expected_region_listB = []

	#create 51 random lists of TFa and TFb sites, with selection being weighted based on the open DNA per gene
	region_value_keys = region_value_dictionary.keys()
	x = 0
	while x < 51:
		random_listA.append(set(np.random.choice(geneID_list, TFa_sites_actual, prob_list)))
		random_listB.append(set(np.random.choice(geneID_list, TFb_sites_actual, prob_list)))
		
		background_expected_listA.append(set(np.random.choice(geneID_list, TFa_GoI_expected_sites, prob_list)))
		background_expected_listB.append(set(np.random.choice(geneID_list, TFb_GoI_expected_sites, prob_list)))

		A = np.random.choice(region_value_keys, TFa_sites_actual, region_problist)
		setA = Set([])
		for number in A:
			setA.add(region_value_dictionary[number])
		random_region_listA.append(setA)		
		
		B = np.random.choice(region_value_keys, TFb_sites_actual, region_problist)
		setB = Set([])
		for number in B:
			setB.add(region_value_dictionary[number])
		random_region_listB.append(setB)

		C = np.random.choice(region_value_keys, TFa_GoI_expected_sites, region_problist)
		setC = Set([])
		for number in C:
			setC.add(region_value_dictionary[number])
		background_expected_region_listA.append(setC)
		
		D = np.random.choice(region_value_keys, TFb_GoI_expected_sites, region_problist)
		setD = Set([])
		for number in D:
			setD.add(region_value_dictionary[number])
		background_expected_region_listB.append(setD)

		x = x + 1
	#create lists to determine the intersection of TFa and TFb sites in a genes open DNA
	#list for shuffle of TFa and TFb sites and find the intersection
	intersection_counter_randomAB = []
	#list for shuffle of TFa only while TFb real sites are maintained
	intersection_counter_randomA = []
	#list for shuffle of TFb only while TFa real sites are maintinaed
	intersection_counter_randomB = []
	#list for intersection of TFa and TFb using random site placement based on expected background site counts of TFa and TFb
	background_intersect_expected_counter = []

	#Start the process of comparing gene lists to find the mean intersection for each condition
	for x in xrange(len(random_listA)):
		randomA_intersect = set(TFb_genelist).intersection(random_listA[x])
		intersection_counter_randomA.append(len(randomA_intersect))
		for y in xrange(len(random_listB)):
			randomB_intersect = set(TFa_genelist).intersection(random_listB[y])
			intersection_counter_randomB.append(len(randomB_intersect))

			Bintersect = random_listA[x].intersection(random_listB[y])
			intersection_counter_randomAB.append(len(Bintersect))

			background_intersect = background_expected_listA[x].intersection(background_expected_listB[y])
			background_intersect_expected_counter.append(len(background_intersect))
	#Values for TFa sites being shuffled while TFb sites from GoI are maintained
        randomA_intersection_mean = np.mean(intersection_counter_randomA)
        randomA_intersection_stdev = np.std(intersection_counter_randomA)

	#Values for TFb sites being shuffled while TFa sites from GoI are maintained
        randomB_intersection_mean = np.mean(intersection_counter_randomB)
        randomB_intersection_stdev = np.std(intersection_counter_randomB)

	#Values for shuffling of TFa and TFb sites based on GoI site numbers
	randomAB_intersection_mean = np.mean(intersection_counter_randomAB)
	randomAB_intersection_stdev = np.std(intersection_counter_randomAB)

	#Values for random intersection of TFa and TFb  using the frequency of TFa and TFb sites in the background region set
	background_intersection_expected_mean = np.mean(background_intersect_expected_counter)
	background_intersection_expected_stdev = np.std(background_intersect_expected_counter)

	#Values for the real intersction of TFa and TFb in the GoI set
	GoIintersect = set(TFa_genelist).intersection(TFb_genelist)
	GoIintersect_size = len(GoIintersect)

	if randomA_intersection_stdev == 0:
		Zscore_randomA = 0
	else:
		Zscore_randomA = (GoIintersect_size - randomA_intersection_mean) / float(randomA_intersection_stdev)
	if randomB_intersection_stdev == 0:
		Zscore_randomB = 0
	else:
		Zscore_randomB = (GoIintersect_size - randomB_intersection_mean) / float(randomB_intersection_stdev)
	if randomAB_intersection_stdev == 0:
		Zscore_randomAB = 0
	else:
		Zscore_randomAB = (GoIintersect_size - randomAB_intersection_mean) / float(randomAB_intersection_stdev)
	if background_intersection_expected_stdev == 0:
		Zscore_background_expected = 0
	else:
		Zscore_background_expected = (GoIintersect_size - background_intersection_expected_mean) / float(background_intersection_expected_stdev)

	Zscore_combined = Zscore_background_expected + ((Zscore_randomA + Zscore_randomB + Zscore_randomAB) / float(3))
	
#Below here is all stuff that needs to be added to gene_cooccurrence anything above is already made in the program

	intersection_region_counter_randomAB = []

	intersection_region_counter_randomA = []

	intersection_region_counter_randomB = []

	background_intersection_region_expected_counter = []

	for x in xrange(len(random_region_listA)):
        	randomA_region_intersect = set(TFb_regionlist).intersection(random_region_listA[x])
	        intersection_region_counter_randomA.append(len(randomA_region_intersect))
        	for y in xrange(len(random_region_listB)):
                	randomB_region_intersect = set(TFa_regionlist).intersection(random_region_listB[y])
	                intersection_region_counter_randomB.append(len(randomB_region_intersect))

        	        Bintersect_region = random_region_listA[x].intersection(random_region_listB[y])
                	intersection_region_counter_randomAB.append(len(Bintersect_region))

	                background_intersect_region = background_expected_region_listA[x].intersection(background_expected_region_listB[y])
        	        background_intersection_region_expected_counter.append(len(background_intersect_region))

	randomA_region_intersection_mean = np.mean(intersection_region_counter_randomA)
	randomA_region_intersection_stdev = np.std(intersection_region_counter_randomA)

	#Values for TFb sites being shuffled while TFa sites from GoI are maintained
	randomB_region_intersection_mean = np.mean(intersection_region_counter_randomB)
	randomB_region_intersection_stdev = np.std(intersection_region_counter_randomB)

	#Values for shuffling of TFa and TFb sites based on GoI site numbers
	randomAB_region_intersection_mean = np.mean(intersection_region_counter_randomAB)
	randomAB_region_intersection_stdev = np.std(intersection_region_counter_randomAB)

	#Values for random intersection of TFa and TFb  using the frequency of TFa and TFb sites in the background region set
	background_region_intersection_expected_mean = np.mean(background_intersection_region_expected_counter)
	background_region_intersection_expected_stdev = np.std(background_intersection_region_expected_counter)

	#Values for the real intersction of TFa and TFb in the GoI set
	GoIintersect_region = set(TFa_regionlist).intersection(TFb_regionlist)
	GoIintersect_region_size = len(GoIintersect_region)

	if randomA_region_intersection_stdev == 0:
	        Zscore_randomA_region = 0
	else:
	        Zscore_randomA_region = (GoIintersect_region_size - randomA_region_intersection_mean) / float(randomA_region_intersection_stdev)
	if randomB_region_intersection_stdev == 0:
        	Zscore_randomB_region = 0
	else:
        	Zscore_randomB_region = (GoIintersect_region_size - randomB_region_intersection_mean) / float(randomB_region_intersection_stdev)
	if randomAB_region_intersection_stdev == 0:
	        Zscore_randomAB_region = 0
	else:
	        Zscore_randomAB_region = (GoIintersect_region_size - randomAB_region_intersection_mean) / float(randomAB_region_intersection_stdev)
	if background_region_intersection_expected_stdev == 0:
	        Zscore_background_expected_region = 0
	else:
        	Zscore_background_expected_region = (GoIintersect_region_size - background_region_intersection_expected_mean) / float(background_region_intersection_expected_stdev)

	Zscore_combined_region = Zscore_background_expected_region + ((Zscore_randomA_region + Zscore_randomB_region + Zscore_randomAB_region) / float(3))

	#Start proximity calculations
	RegionSizeList = []
	DistDict = {}
	
	within50bp_count = 0
	outside50bp_count = 0
	
	forward_min_distance = 0
	forward_middle_distance = 0
	forward_farthest_distance = 0

	reverse_min_distance = 0
	reverse_middle_distance = 0
	reverse_farthest_distance = 0

	motiflengthA = 0
	motiflengthB = 0
	maxoverlap = 0

	#print "Region Cooccurrence: ", GoIintersect_region_size
	if GoIintersect_region_size == 0:
                proximity_within50bp_pvalue, forward_min_distance_pvalue, forward_middle_distance_pvalue, forward_farthest_distance_pvalue, reverse_min_distance_pvalue, reverse_middle_distance_pvalue, reverse_farthest_distance_pvalue = 1, 1, 1, 1, 1, 1, 1

	else:
		for region in GoIintersect_region:
			regionSize = region[2] - region[1] + 1
			motiflengthA = TFaRegion_SiteDict[region][0][2] - TFaRegion_SiteDict[region][0][1] + 1
			motiflengthB = TFbRegion_SiteDict[region][0][2] - TFbRegion_SiteDict[region][0][1] + 1

			#Determine how much of one TF site is allowed to overlap with a different TF site
			maxoverlap = min(((motiflengthA/2)-1), ((motiflengthB/2)-1))
			#Entries for each set of TF-TF spacings
			DistDict['spacings'] = 100 + 2*maxoverlap

			for siteA in TFaRegion_SiteDict[region]:
				for siteB in TFbRegion_SiteDict[region]:
					RegionSizeList.append(regionSize)
					middleA = siteA[1] + (motiflengthA/2)
					middleB = siteB[1] + (motiflengthB/2)
					#print "middleA, middleB:", middleA, middleB
					#This section checks if the TF-TF sites fall within 50 bp of each other and doesn't overlap too much
        	                        if siteA[5] == 'p':
                	                        if (middleA - middleB) <= 0:
                        	                        dist = '+'
                                	                spacing = siteB[1] - siteA[2]
							#print "spacing: ", spacing
                                        	        if spacing >= -1*maxoverlap and spacing <= 50:
                                                	        if spacing < 11:
	                                                                forward_min_distance = forward_min_distance + 1
        	                                                elif spacing < 26:
                	                                                forward_middle_distance = forward_middle_distance + 1
                        	                                else:
                                	                                forward_farthest_distance = forward_farthest_distance + 1


	                                                        if (spacing, dist) not in DistDict:
        	                                                        DistDict[(spacing, dist)] = 1
                	                                        else:
                        	                                        DistDict[(spacing, dist)] += 1
                                	                        within50bp_count = within50bp_count + 1
                                        	        else:
                                                	        outside50bp_count = outside50bp_count + 1
	                                                        continue
        	                                elif (middleA - middleB) > 0:
                	                                dist = '-'
                        	                        spacing = siteA[1] - siteB[2]
							#print "spacing: ", spacing
                                	                if spacing >= -1*maxoverlap and spacing <= 50:
                                        	                if spacing < 11:
                                                	                reverse_min_distance = reverse_min_distance + 1
	                                                        elif spacing < 26:
        	                                                        reverse_middle_distance = reverse_middle_distance + 1
                	                                        else:
                        	                                        reverse_farthest_distance = reverse_farthest_distance + 1
	
        	                                                if (spacing, dist) not in DistDict:
                	                                                DistDict[(spacing, dist)] = 1
                        	                                else:
                                	                                DistDict[(spacing, dist)] += 1
	                                                        within50bp_count = within50bp_count + 1
        	                                        else:
                	                                        outside50bp_count = outside50bp_count + 1
                        	                                continue
	
	                                elif siteA[5] == 'm':
        	                                if (middleA - middleB) < 0:
                	                                dist = '-'
                        	                        spacing = siteB[1] - siteA[2]
							#print "spacing: ", spacing
                                	                if spacing >= -1*maxoverlap and spacing <= 50:

	                                                       if spacing < 11:
        	                                                       reverse_min_distance = reverse_min_distance + 1
                	                                       elif spacing < 26:
                        	                                       reverse_middle_distance = reverse_middle_distance + 1
                                	                       else:
	                                                               reverse_farthest_distance = reverse_farthest_distance + 1
	
        	                                               if (spacing, dist) not in DistDict:
                	                                               DistDict[(spacing, dist)] = 1
                        	                               else:
                                	                               DistDict[(spacing, dist)] += 1
                                        	               within50bp_count = within50bp_count + 1
	                                                else:
        	                                        	outside50bp_count = outside50bp_count + 1
                	                                        continue
                        	                elif (middleA - middleB) >= 0:
                                	                dist = '+'
                                        	        spacing = siteA[1] - siteB[2]
							#print "spacing: ", spacing
	                                                if spacing >= -1*maxoverlap and spacing <= 50:
	
        	                                                if spacing < 11:
                	                                                forward_min_distance = forward_min_distance + 1
                        	                                elif spacing < 26:
                                	                                forward_middle_distance = forward_middle_distance + 1
                                        	                else:
                                                	                forward_farthest_distance = forward_farthest_distance + 1


	                                                        if (spacing, dist) not in DistDict:
        	                                                        DistDict[(spacing, dist)] = 1
                	                                        else:
                        	                                        DistDict[(spacing, dist)] += 1
                                	                        within50bp_count = within50bp_count + 1
                                        	        else:
	                                                        outside50bp_count = outside50bp_count + 1
        	                                                continue
	
	
        	number = 0

	        random_within50bp_count = 0
        	random_outside50bp_count = 0

	        random_forward_min_distance = 0
	        random_forward_middle_distance = 0
	        random_forward_farthest_distance = 0

        	random_reverse_min_distance = 0
	        random_reverse_middle_distance = 0
	        random_reverse_farthest_distance = 0

        	while number < 1001:
	                number = number + 1
        	        for RegionSize in RegionSizeList:
                	        #Determines how much space within the open region was available for a potential interaction
	                        Available_Space = RegionSize - motiflengthA - motiflengthB + maxoverlap

        	                #Determines a random spacing and orientation of TFa and TFb within the given available space of an open region
                	        random_spacing = random.randint(-1*maxoverlap, Available_Space)
                        	random_orientation = random.choice(['-', '+'])

	                        if random_spacing < 51:
        	                        random_within50bp_count = random_within50bp_count + 1

                	                if random_spacing < 11:
                        	                if random_orientation == '+':
                                	                random_forward_min_distance = random_forward_min_distance + 1
                                        	else:
	                                                random_reverse_min_distance = random_reverse_min_distance + 1
        	                        elif random_spacing < 26:
                	                        if random_orientation == '+':
                        	                        random_forward_middle_distance = random_forward_middle_distance + 1
                                	        else:
	                                                random_reverse_middle_distance = random_reverse_middle_distance + 1

        	                        elif random_spacing < 51:
                	                        if random_orientation == '+':
                        	                        random_forward_farthest_distance = random_forward_farthest_distance + 1
                                	        else:
                                        	        random_reverse_farthest_distance = random_reverse_farthest_distance + 1
	                        else:
        	                        random_outside50bp_count = random_outside50bp_count + 1
		#print "within 50bp, random within 50bp, outside 50 bp, random outside 50bp: ", within50bp_count, random_within50bp_count, outside50bp_count, random_outside50bp_count
		if (random_within50bp_count == 0) or (random_outside50bp_count ==0):
			proximity_within50bp_pvalue = 1			        
		else:
			proximity_within50bp_pvalue_nonlog = scipy.stats.chi2_contingency([[within50bp_count, outside50bp_count], [random_within50bp_count, random_outside50bp_count]], lambda_="log-likelihood")[1]
		        if proximity_within50bp_pvalue_nonlog == 0:
        		        proximity_within50bp_pvalue =0
		        else:
        		        proximity_within50bp_pvalue = math.log10(proximity_within50bp_pvalue_nonlog)


	        totalcount = within50bp_count + outside50bp_count
	        random_totalcount = random_within50bp_count + random_outside50bp_count

		if (totalcount - forward_min_distance == 0) or (random_totalcount - random_forward_min_distance == 0) or random_forward_min_distance == 0:
			forward_min_distance_pvalue = 1
		else:
			#print forward_min_distance, totalcount - forward_min_distance, random_forward_min_distance, random_totalcount - random_forward_min_distance
	        	forward_min_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[forward_min_distance, totalcount - forward_min_distance], [random_forward_min_distance, random_totalcount - random_forward_min_distance]], lambda_="log-likelihood")[1]
		        if forward_min_distance_pvalue_nonlog == 0:
        		        forward_min_distance_pvalue = 0
		        else:
        		        forward_min_distance_pvalue = math.log10(forward_min_distance_pvalue_nonlog)

		if (totalcount - forward_middle_distance == 0) or (random_totalcount - random_forward_middle_distance == 0) or random_forward_middle_distance == 0:
			forward_middle_distance_pvalue = 1
		else:
			#print forward_middle_distance, totalcount - forward_middle_distance, random_forward_middle_distance, random_totalcount - random_forward_middle_distance
		        forward_middle_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[forward_middle_distance, totalcount - forward_middle_distance], [random_forward_middle_distance, random_totalcount - random_forward_middle_distance]], lambda_="log-likelihood")[1]
        		if forward_middle_distance_pvalue_nonlog == 0:
	        	        forward_middle_distance_pvalue = 0
	        	else:
		                forward_middle_distance_pvalue = math.log10(forward_middle_distance_pvalue_nonlog)

		if (totalcount - forward_farthest_distance == 0) or (random_totalcount - random_forward_farthest_distance == 0) or random_forward_farthest_distance == 0:
			forward_farthest_distance_pvalue = 1
		else:
			#print forward_farthest_distance, totalcount - forward_farthest_distance, random_forward_farthest_distance, random_totalcount - random_forward_farthest_distance
        		forward_farthest_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[forward_farthest_distance, totalcount - forward_farthest_distance], [random_forward_farthest_distance, random_totalcount - random_forward_farthest_distance]], lambda_="log-likelihood")[1]
		        if forward_farthest_distance_pvalue_nonlog == 0:
        		        forward_farthest_distance_pvalue = 0
	        	else:
	        	        forward_farthest_distance_pvalue = math.log10(forward_farthest_distance_pvalue_nonlog)

		if (totalcount - reverse_min_distance == 0) or (random_totalcount - random_reverse_min_distance == 0) or random_reverse_min_distance == 0:
			reverse_min_distance_pvalue = 1
		else:
			#print reverse_min_distance, totalcount - reverse_min_distance, random_reverse_min_distance, random_totalcount - random_reverse_min_distance
		        reverse_min_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[reverse_min_distance, totalcount - reverse_min_distance], [random_reverse_min_distance, random_totalcount - random_reverse_min_distance]], lambda_="log-likelihood")[1]
	        	if reverse_min_distance_pvalue_nonlog == 0:
        	        	reverse_min_distance_pvalue = 0
	        	else:
	        	        reverse_min_distance_pvalue = math.log10(reverse_min_distance_pvalue_nonlog)

		if (totalcount - reverse_middle_distance == 0) or (random_totalcount - random_reverse_middle_distance == 0) or random_reverse_middle_distance == 0:
			reverse_middle_distance_pvalue = 1
		else:
			#print reverse_middle_distance, totalcount - reverse_middle_distance, random_reverse_middle_distance, random_totalcount - random_reverse_middle_distance
		        reverse_middle_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[reverse_middle_distance, totalcount - reverse_middle_distance], [random_reverse_middle_distance, random_totalcount - random_reverse_middle_distance]], lambda_="log-likelihood")[1]
        		if reverse_middle_distance_pvalue_nonlog == 0:
		                reverse_middle_distance_pvalue = 0
	        	else:
        	        	reverse_middle_distance_pvalue = math.log10(reverse_middle_distance_pvalue_nonlog)

		if (totalcount - reverse_farthest_distance == 0) or (random_totalcount - random_reverse_farthest_distance == 0) or random_reverse_farthest_distance == 0:
			reverse_farthest_distance_pvalue = 1
		else:
			#print reverse_farthest_distance, totalcount - reverse_farthest_distance, random_reverse_farthest_distance, random_totalcount - random_reverse_farthest_distance
		        reverse_farthest_distance_pvalue_nonlog = scipy.stats.chi2_contingency([[reverse_farthest_distance, totalcount - reverse_farthest_distance], [random_reverse_farthest_distance, random_totalcount - random_reverse_farthest_distance]], lambda_="log-likelihood")[1]
        		if reverse_farthest_distance_pvalue_nonlog == 0:
                		reverse_farthest_distance_pvalue = 0
		        else:
        		        reverse_farthest_distance_pvalue = math.log10(reverse_farthest_distance_pvalue_nonlog)

	return TFa_sites_actual, TFb_sites_actual, GoIintersect_size, randomA_intersection_mean, randomB_intersection_mean, randomAB_intersection_mean, background_intersection_expected_mean, Zscore_randomA, Zscore_randomB, Zscore_randomAB, Zscore_background_expected, GoIintersect, Zscore_combined, GoIintersect_region_size, randomA_region_intersection_mean, randomB_region_intersection_mean, randomAB_region_intersection_mean, background_region_intersection_expected_mean, Zscore_randomA_region, Zscore_randomB_region, Zscore_randomAB_region, Zscore_background_expected_region, GoIintersect_region, Zscore_combined_region, within50bp_count, outside50bp_count, DistDict, proximity_within50bp_pvalue, forward_min_distance_pvalue, forward_middle_distance_pvalue, forward_farthest_distance_pvalue, reverse_min_distance_pvalue, reverse_middle_distance_pvalue, reverse_farthest_distance_pvalue

def TSS_TF_parser(TSSFileName, BackgroundDNADict, GeneListFileName, in_directory, upstream, downstream):
	print TSSFileName
	TFlist = []
	#File that has all of the used TSS locations
        TSSfile = open(TSSFileName, "r")
	#File that has all of the Genes to be considered as Genes of Interest
        GeneListFile = open(GeneListFileName, "r")
        BackgroundDict = {}
	GenesOfInterestDict = {}

	BackgroundTFDict = {}
	GenesOfInterestTFDict = {}

	GenesOfInterest = Set([])
	BackgroundGenes = Set([])
        InputGeneCount = 0

	#Counters for # of genes in Genes of Interest and Background
        MatchGeneCount = 0
	BackgroundGeneCount = 0

	#Counters for the total amount of open DNA in Genes of Interest or Background
	BackgroundDNACount = 0
	GoIDNACount = 0

	tempset = Set([])

	#Generates the temporary Genes of Interest list,
	#need to see if there is any open DNA at the genes before fully considering it in the analysis
        for geneline in GeneListFile:
                disregard = geneline.strip()
                GeneLine = disregard.split()
                GeneName = str(GeneLine[0])
                Gene = GeneName.upper()
		tempset.add(Gene)

	#This loop looks generates the Gene of Interest Set and Background Set, 
	#BackgroundDNADict identifies which Genes have open DNA within the provided upstream and downstream regions
	#If there is no open DNA then the gene is not added to either list
	for gene, DNAcount in BackgroundDNADict.iteritems():
		if gene in tempset:
			GenesOfInterest.add(gene)
			GoIDNACount = GoIDNACount + DNAcount
		else:
			BackgroundGenes.add(gene)
			BackgroundDNACount = BackgroundDNACount + DNAcount


	print "Start making the TSSDict"	
	TSSDict = {'chr1':[], 'chr2':[], 'chr3':[], 'chr4':[], 'chr5':[], 'chr6':[], 'chr7':[], 'chr8':[], 'chr9':[], 'chr10':[], 'chr11':[], 'chr12':[], 'chr13':[], 'chr14':[], 'chr15':[], 'chr16':[], 'chr17':[], 'chr18':[], 'chr19':[], 'chr20':[], 'chr21':[], 'chr22':[], 'chrX':[], 'chrY':[]}

	geneindict = 0
	#Create a Dictionary that contains all of the TSS of genes in this run
	#The keys are chromosome number and the values are [start of regulatory region, end of regulatory region, TSS value, and gene name]
	for TSSline in TSSfile:
		remove = TSSline.strip()
		TSS = remove.split()
		TSS_gene = str(TSS[0])
		if TSS_gene not in BackgroundGenes and TSS_gene not in GenesOfInterest:
			continue
		else:
			geneindict = geneindict+1
			TSS_chromosome = str(TSS[1])
			if TSS_chromosome not in TSSDict.keys():
				TSSDict.update({TSS_chromosome:[]})
			TSS_location = int(TSS[2])
			TSS_strand = str(TSS[3])
                	BackgroundGeneCount = BackgroundGeneCount + 1
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
			TSSDict[TSS_chromosome].append((RegulatoryRange[1], RegulatoryRange[2], TSS_location, TSS_gene))
	print "Finish making TSSDict"
	for chrom in TSSDict.keys():
		print chrom, len(TSSDict[chrom])
	print "Genes in TSS dictionary: ", geneindict
	
	genecount = 0
	#Sorty the TSS dictionary values from smallest to largest to help with time efficiency in later checks
	for TSS_chromosome, sitelist in TSSDict.iteritems():
		sortedlist = sorted(sitelist, key = lambda x: x[1])
		TSSDict[TSS_chromosome] = sortedlist
		genecount = genecount + len(sitelist)
	print "The total genecount is: ", genecount
	
 	TFDict = {}
	BackgroundTFDict = {}
	GenesOfInterestTFDict = {}	
	TFList = []
	NoGoIList = []
	#Open the transcription factor putative site directory and loop over each transcription factor
	for TF in os.listdir(in_directory):
		filelocation = in_directory + TF
		#Find the overlaps betwen the transcription factor sites and the regulatory regions around the TSS
		#Returns a dictionary with TFs as keys and values as a dictionary with Genes as key and site lists as values
                output_dictionary = Overlap(filelocation, TSSDict, upstream, downstream)

		if output_dictionary.keys() == []:
			print "The following file has no binding sites: ", TF
			continue
		TFList.append(TF)
		TFDict[TF] = output_dictionary
		#This loop separates the dictionaries into Background sites and Gene of Interest sites
		for BoundGene, sitelist in TFDict[TF].iteritems():

			if BoundGene not in BackgroundDict.keys() and BoundGene not in GenesOfInterestDict.keys():

				if BoundGene not in GenesOfInterest:
					BackgroundDict[BoundGene] = {TF:sitelist}

					if TF not in BackgroundTFDict.keys():
						BackgroundTFDict[TF] = {BoundGene:sitelist}

					else:
						BackgroundTFDict[TF].update({BoundGene:sitelist})

				else:
					GenesOfInterestDict[BoundGene] = {TF:sitelist}

					if TF not in GenesOfInterestTFDict.keys():
						GenesOfInterestTFDict[TF] = {BoundGene:sitelist}

					else:
						GenesOfInterestTFDict[TF].update({BoundGene:sitelist})

			else:
				if BoundGene not in GenesOfInterest:
					BackgroundDict[BoundGene].update({TF:sitelist})
					if TF not in BackgroundTFDict.keys():
						BackgroundTFDict[TF] = {BoundGene:sitelist}
					else:
						BackgroundTFDict[TF].update({BoundGene:sitelist})
				else:
					GenesOfInterestDict[BoundGene].update({TF:sitelist})
                                        if TF not in GenesOfInterestTFDict.keys():
                                                GenesOfInterestTFDict[TF] = {BoundGene:sitelist}
                                        else:
                                                GenesOfInterestTFDict[TF].update({BoundGene:sitelist})

	MatchGeneCount = len(GenesOfInterestDict.keys())
	BackgroundGeneCount = len(BackgroundDict.keys())
	print "Genes used in provided list: ", MatchGeneCount
	print "Background Genes: ", BackgroundGeneCount

	#Return all of the relevant information about the transcription factors, Gene of Interest gene # and DNA count, and Background gene # and DNA count
	return (GenesOfInterestDict, TFList, MatchGeneCount, BackgroundDict, len(GenesOfInterest), len(BackgroundGenes), BackgroundTFDict, GenesOfInterestTFDict, GoIDNACount, BackgroundDNACount)

def KL_distance_statistic(BackgroundGeneDict, BackgroundGeneList, GoISize, GoISiteList):
	#This section determines the Kullback-Leibler distance/divergence of the Genes of Interest from the Background
	#50 random data sets, of the same size as the Gene of Interest set, are generated from the Background
	#The KL-distance is generated for each random compared to random to generate a distribution of KL-distances
	#Then the KL-distance between the Genes of Interest and the random background sets is calculated
	#The average KL-distance of the GoI vs. Background is then compared to the distribution of Background vs. Background values to give a p-value
	BackgroundGeneSize = len(BackgroundGeneList)
	list_of_lists = []
	number = 0
	#Generate a list of gene lists, representing the random background gene sets
	while number < 50:
		number = number + 1
		random_selection = random.sample(range(BackgroundGeneSize), GoISize)
		BackgroundSiteList = []
		for y in random_selection:
			GeneName = BackgroundGeneList[y]
			if GeneName not in BackgroundGeneDict.keys():
				BackgroundSiteList.append(0)
			else:
				size = len(BackgroundGeneDict[GeneName])
				BackgroundSiteList.append(size)
		list_of_lists.append(BackgroundSiteList)		
	
	#This section generates the histograms that will be used to determine the KL-distance values	
	maximum = max(max(y) for y in list_of_lists)

	defined_bins = range(0, maximum+2)

	GoI_entropy_distribution = []
	Background_entropy_distribution = []
	histogram_GoI = np.histogram(GoISiteList, bins = defined_bins)
	histogram_GoI_pseudo = histogram_GoI[0]+1

	#Loops that generate the KL-distance values
	for samplea in list_of_lists:
		#This section generates the GoI vs. a random background set KL-Distance and adds it to a list of values
		histogram_samplea = np.histogram(samplea, bins=defined_bins)
                histogram_samplea_pseudo = histogram_samplea[0]+1
                kl_value_goi = scipy.stats.entropy(histogram_GoI_pseudo, histogram_samplea_pseudo)
                GoI_entropy_distribution.append(kl_value_goi)

                for sampleb in list_of_lists:
                        if samplea == sampleb:
				pass
                        else:
				#This section generates the Random Background vs. Random Background KL-Distance values
                                histogram_sampleb = np.histogram(sampleb, bins=defined_bins)
                                histogram_sampleb_pseudo = histogram_sampleb[0]+1
                                kl_value_controls = scipy.stats.entropy(histogram_samplea_pseudo, histogram_sampleb_pseudo)
                                Background_entropy_distribution.append(kl_value_controls)

	Background_distribution_size = len(Background_entropy_distribution)
	GoI_average = reduce(lambda x, y: x + y, GoI_entropy_distribution)/float(len(GoI_entropy_distribution))
        Background_average = reduce(lambda x, y: x + y, Background_entropy_distribution)/float(Background_distribution_size)

        larger_values = 0

	#compares the GoI average KL-Distance to each of the random vs. random background KL-Distances
	#this comparison is what is used to generate the KL-Distance p-value
        for value in Background_entropy_distribution:
                if value >= GoI_average:
                        larger_values = larger_values + 1
                else:
                        pass
	#Return the p-value
        return larger_values/float(Background_distribution_size)
			
def ProcessCLI(args):
        print args
	in_directory = "/N/u/jubudka/Mason/MyProgramFiles/RWPE_ExampleTFs/"
        in_filename = "hg19_TSS_list.txt"
	in_genelist = "ETS1_Down_Up_Genelist.txt"
	in_background = "hg19_background_withOpenRegion.txt"
	out_countfilename = "Test_TFcounts_1_13_17.txt"
	Gene_Cooccurfile = "Test_Gene_Cooccur_1_13_17.txt"
	Region_Cooccurfile = "Test_Region_Cooccur_1_13_17.txt"
	Proximity_file = "Test_Proximity_1_13_17.txt"
	upstream = 20000
	downstream = 10000
	similarfile = "/N/u/jubudka/Mason/matalign/Matalign_Similar_TF_Output_2.txt"

        for i in xrange(len(args)):
                if args[i] == "-i":
                        in_filename = args[i+1]
		elif args[i] == "-d":
			in_directory = args[i+1]
                elif args[i] == "-o":
                       out_countfilename = args[i+1]
		elif args[i] == "-u":
			upstream = int(args[i+1])
		elif args[i] == "-do":
			downstream = int(args[i+1])
		elif args[i] == "-g":
			in_genelist = args[i+1]
		elif args[i] == "-b":
			in_background = args[i+1]
		elif args[i] == "-gc":
			Gene_Cooccurfile = args[i+1]
		elif args[i] == "-rc":
			Region_Cooccurfile = args[i+1]
		elif args[i] == "-pc":
			Proximity_file = args[i+1]
		elif args[i] == "-s":
			similarfile = args[i+1]

	print "input filename is ", in_filename
	print "Directory location is ", in_directory
#	print "output filename is ", out_filename
	print "distance upstream of promoter ", upstream
	print "distance downstream of promoter ", downstream
	print "Gene list of interest ", in_genelist
	print "Background Gene List: ", in_background
		
	similarTFfile = open(similarfile, "r")

	similarityDict = {}

	print "Start Similarity Dictionary"
	for lines in similarTFfile:
                remove = lines.strip()
                pair = remove.split()
                if pair[0] in similarityDict.keys():
                        if pair[1] in similarityDict[pair[0]]:
                                pass
                        else:
                                similarityDict[pair[0]].append(pair[1])
                else:
                        similarityDict[pair[0]] = [pair[1]]


                if pair[1] in similarityDict.keys():
                        if pair[0] in similarityDict[pair[1]]:
                                pass
                        else:
                                similarityDict[pair[1]].append(pair[0])
                else:
                        similarityDict[pair[1]] = [pair[0]]
	print "Finish Similarity Dictionary"

	BackgroundDNA = open(in_background, "r")
	BackgroundDNADict = {}
	BackgroundRegionDict = {}

	print "Start Background File Parsing"

        position = 1
        for backgroundline in BackgroundDNA:
                if position == 1:
                        position = 2
                        continue
                disregard = backgroundline.strip()
                GeneLine = disregard.split('\t')
                GeneName = str(GeneLine[0])
                Gene = GeneName.upper()
		Regions = GeneLine[2].split('//')
		RegionList = [eval(x) for x in Regions]
		
		BackgroundDNADict[Gene] = int(GeneLine[1])
		BackgroundRegionDict[Gene] = RegionList
		
	print "Finish Background File Parsing"

	
	print "Start TSS File Parsing"
	#Read through the input files and generate the dictionaries for data analysis 
	BindingInformation = TSS_TF_parser(in_filename, BackgroundDNADict, in_genelist, in_directory, upstream, downstream)
	print "Finish TSS File Parsing"

	#Variables generated from parsing the input files
	GoIGeneDict = BindingInformation[0]
	TFlist = BindingInformation[1]
	GeneCount = BindingInformation[2]
	BackgroundGeneDict = BindingInformation[3] 	

	GoITotalGeneCount = BindingInformation[4]
	BackgroundTotalGeneCount = BindingInformation[5]

	GoITFDict = BindingInformation[7]
	BackgroundTFDict = BindingInformation[6]

	#Value for all of the open DNA surrounding the GoI
	GoITotalDNACount = int(BindingInformation[8])
	print "DNA in GeneList: ", GoITotalDNACount
	#Value for all of the open DNA surrounds the background genes
	BackgroundTotalDNACount = int(BindingInformation[9])
	print "DNA in Background: ", BackgroundTotalDNACount
	#Ratio of GoI available DNA to that of Background Genes available DNA
	DNAratio = GoITotalDNACount/float(BackgroundTotalDNACount)
	print "Ratio DNA Genelist to DNA Background: ", DNAratio

	BackgroundGeneList = BackgroundGeneDict.keys()
	GoIGeneList = GoIGeneDict.keys()
	#Generate a list of all of the genes that have some amount of open DNA
	GeneList = GoIGeneList + BackgroundGeneList
	# number of genes in GoI
	GoISize = len(GoIGeneList)
	# number of genes in Background
	BSize = len(BackgroundGeneList)

	print "Finish TSS File Parsing"

	#Keep a running list of information about all of the TFs 
	TFcountInfo = []
	
	#Open the output file for the transcription factor overrepresentation output
	a = open(out_countfilename, "w")
	a.write('TF name' + '\t' + 'motif length' + '\t' + '# sites in genelist' + '\t' + '# Genes with site' + '\t' + 'Sites/Gene' + '\t' + '# sites in background' + '\t' + '# Genes with sites in background' + '\t' + 'Sites/Gene background' + '\t' + 'Site to Gene Factor' + '\t' + 'KL_p_value' + '\t' + 'Corrected KL pval' + '\t' + 'Zscore by DNA#' + '\t' + 'p-value from Zscore' + '\t' + 'Corrected normPval' + '\t' + 'Fishers Exact p-value' + '\t' + 'Corrected Fishers' + '\t' + 'Genes in list with site' + '\n')

        b = open(Gene_Cooccurfile, "w")
        b.write('Trans Factor A' + '\t' + 'TransFactor B' + '\t' + '# sites TFa' + '\t' + '# sites TFb' + '\t' + '# Gene Co-occurrence' + '\t' + 'Intersect_randomA' + '\t' + 'Intersect_randomB' + '\t' + 'Intersect_randomAB' + '\t' + 'Intersect_Background_Expectation' + '\t' + 'Zscore_randomA' + '\t' + 'Zscore_randomB' + '\t' + 'Zscore_randomAB' + '\t' + 'Zscore_Background_Expectation' + '\t' + 'Combined Zscore' + '\t' + 'Co-occurring genes' + '\n') 
	Gene_CooccurList = []

	c = open(Region_Cooccurfile, "w")
	c.write('Trans Factor A' + '\t' + 'TransFactor B' + '\t' + '# sites TFa' + '\t' + '# sites TFb' + '\t' + '# Region Co-occurrence' + '\t' + 'Intersect_randomA' + '\t' + 'Intersect_randomB' + '\t' + 'Intersect_randomAB' + '\t' + 'Intersect_Background_Expectation' + '\t' + 'Zscore_randomA' + '\t' + 'Zscore_randomB' + '\t' + 'Zscore_randomAB' + '\t' + 'Zscore_Background_Expectation' + '\t' + 'Combined Zscore' + '\t' + 'Co-occurring regions' + '\n')
	Region_CooccurList = []

	d = open(Proximity_file, "w")
	d.write('TFa' + '\t' + 'TFb' + '\t' + 'Count_within_50bp' + '\t' + 'Count Outside 50 bp' + '\t' + 'p-value within 50 bp' + '\t' + 'p-value TFa-TFb min distance' + '\t' + 'p-value TFa-TFb middle distance' + '\t' + 'p-value TFa-TFb farthest distance' + '\t' + 'p-value TFb-TFa min distance' + '\t' + 'p-value TFb-TFa middle distance' + '\t' + 'p-value TFb-TFa farthest distance' + '\n')
	Proximity_CooccurList = []

	print "Start TF analysis"
	#Run through all the transcription factors provided	
	for TF in TFlist:	
		GoISiteList = []
		Genelist = []
		#Counter for the amount of DNA in GoI that is assigned to the TF's binding sites for all GoI
 		GoITFdna = 0
		GoIsitecount = 0
		#Counter for the amount of DNA in background that is assigned to the TF's binding sites for all background genes
		BackgroundTFdna = 0 
		BackgroundGenes = []
		#Counter for number of TF binding sites in DNA near gene in background
		Backsitecount = 0
		motiflength = 0
		#Counter for number of GoI that have at least 1 TF binding site
		numGenes = 0
		#Counter for number of background genes that have at least 1 TF binding site
		BnumGenes = 0
		Zscore = 0

		if TF in BackgroundTFDict and TF in GoITFDict:
			#Check to confirm the TF has at least 1 site in the background and GoI sets
			pass
		else:
			#If the TF isn't in both Dictionaries then there is no way it is overrepresented and the rest of the loop is unncecessary
			continue

		for gene in GeneList:
			if gene in BackgroundGeneList:
				BackgroundGenes.append(gene)
				try:
					#This confirms that this gene has at least 1 TF binding site nearby
					trial = BackgroundTFDict[TF]
				except:
					#If the gene doesn't have at least 1 TF site then skip the loop
					continue

				if gene in BackgroundTFDict[TF].keys():
					BnumGenes = BnumGenes + 1
					motiflength = BackgroundTFDict[TF][gene][0][4]
					Bsitecount = len(BackgroundTFDict[TF][gene])
					Backsitecount = Backsitecount + Bsitecount
					BackgroundTFdna = BackgroundTFdna + (motiflength * Bsitecount)
				else:
					continue
			else:
				try:
					#This confirms that this gene has at least 1 TF binding site nearby
					trial = GoITFDict[TF]
				except:
					#If the gene doesn't have at least 1 TF site then skip the loop
					GoISiteList.append(0)
					continue

				if gene not in GoITFDict[TF]:
					GoISiteList.append(0)
					continue
				else:
					Genelist.append(gene)
					sitecount = len(GoITFDict[TF][gene])
					motiflength = GoITFDict[TF][gene][0][4]
					GoISiteList.append(sitecount)
					GoITFdna = GoITFdna + (motiflength * sitecount)
					numGenes = numGenes + 1
					GoIsitecount = GoIsitecount + sitecount
		
		if len(GoITFDict[TF]) == 0:
			Zscore = 0
		else:
			#General approach borrowed from oPossum protocol, based on binomial distribution of binding site DNA differences
			mu = BackgroundTFdna * DNAratio
			P = BackgroundTFdna / float(BackgroundTotalDNACount)
			sigma = math.sqrt((GoITotalDNACount * P) * (1 - P))
			Zscore = ((GoITFdna - 0.5 - mu)/sigma)
					
		if len(BackgroundTFDict[TF]) == 0 and len(GoITFDict[TF]) == 0:
			#If the TF isn't in the dictionary then give a p-value of 1
			summary_statistic = 1.0
		else:
			#Calculation of a p-value from the Kullback-Leibler Divergence. Its an estimate for the difference between probability distributions
			summary_statistic = KL_distance_statistic(BackgroundTFDict[TF], BackgroundGenes, GoISize, GoISiteList)

                overrepresentation_information = [TF, int(motiflength), GoIsitecount, numGenes, Backsitecount, BnumGenes, Genelist, summary_statistic, Zscore]
		TFcountInfo.append(overrepresentation_information)

		#Generate co-occurrence statistics for pairs of transcription factors	
		for TFb in TFlist:
			if TFb in BackgroundTFDict and TFb in GoITFDict:
				pass
			else:
				continue

			if TF == TFb:
				continue

                        if TF in similarityDict.keys() and TFb in similarityDict[TF]:
                                continue
			else:
				if TFb in similarityDict.keys():
					similarityDict[TFb].append(TF)
				else:
					similarityDict[TFb] = [TF]

			#(TFa, TFb, BackgroundTFDict, TFDict, DNADict, BackgroundTotalDNACount, GoITotalDNACount, GoIGenes, BackgroundRegionDict)
			gene_cooccur_weighted = Gene_Cooccurrence(TF, TFb, BackgroundTFDict, GoITFDict, BackgroundDNADict, BackgroundTotalDNACount, GoITotalDNACount, GoIGeneList, BackgroundRegionDict)
			
#TFa_sites_actual, TFb_sites_actual, GoIintersect_size, randomA_intersection_mean, randomB_intersection_mean, randomAB_intersection_mean, background_intersection_expected_mean, Zscore_randomA, Zscore_randomB, Zscore_randomAB, Zscore_background_expected, GoIintersect, Zscore_combined, GoIintersect_region_size, randomA_region_intersection_mean, randomB_region_intersection_mean, randomAB_region_intersection_mean, background_region_intersection_expected_mean, Zscore_randomA_region, Zscore_randomB_region, Zscore_randomAB_region, Zscore_background_expected_region, GoIintersect_region, Zscore_combined_region, within50bp_count, outside50bp_count, DistDict, proximity_within50bp_pvalue, forward_min_distance_pvalue, forward_middle_distance_pvalue, forward_farthest_distance_pvalue, reverse_min_distance_pvalue, reverse_middle_distance_pvalue, reverse_farthest_distance_pvalue
			# GoI_intersect_size, randomA_intersectAVG, randomB_intersectAVG, randomAB_intersectAVG, Zscore_randomA, Zscore_randomB, Zscore_randomAB, Actual_Intersectlist

			Gene_CooccurList.append([TF, TFb, gene_cooccur_weighted[0], gene_cooccur_weighted[1], gene_cooccur_weighted[2], gene_cooccur_weighted[3], gene_cooccur_weighted[4], gene_cooccur_weighted[5], gene_cooccur_weighted[6], gene_cooccur_weighted[7], gene_cooccur_weighted[8], gene_cooccur_weighted[9], gene_cooccur_weighted[10], gene_cooccur_weighted[11], gene_cooccur_weighted[12]])


			Region_CooccurList.append([TF, TFb, gene_cooccur_weighted[0], gene_cooccur_weighted[1], gene_cooccur_weighted[13], gene_cooccur_weighted[14], gene_cooccur_weighted[15], gene_cooccur_weighted[16], gene_cooccur_weighted[17], gene_cooccur_weighted[18], gene_cooccur_weighted[19], gene_cooccur_weighted[20], gene_cooccur_weighted[21], gene_cooccur_weighted[22], gene_cooccur_weighted[23]])			


			Proximity_CooccurList.append([TF, TFb, gene_cooccur_weighted[24], gene_cooccur_weighted[25], gene_cooccur_weighted[27], gene_cooccur_weighted[28], gene_cooccur_weighted[29], gene_cooccur_weighted[30], gene_cooccur_weighted[31], gene_cooccur_weighted[32], gene_cooccur_weighted[33]])

	print "Finish TF analysis"

	print "Start TF and TF pair file generation" 

	number = len(TFcountInfo)
	#This section combines all of the statistical information about the Transcription Factor sets
	#It is used to decide what to output for TF overrepresentation
	for TFinfo in sorted(TFcountInfo, key=lambda occurrence: occurrence[7]):
		geneString = ", ".join(TFinfo[6])
		TFname = TFinfo[0]
		Motif = TFinfo[1]
		numSites = TFinfo[2]
		numGenes = TFinfo[3]
		if numSites == 0:
			ratio = 1
		else:
			ratio = float(numSites)/numGenes
		BnumSites = TFinfo[4]
		BnumGenes = TFinfo[5]
		if BnumSites == 0:
			Bratio = 1
		else:
			Bratio = float(BnumSites)/BnumGenes
		ratioBratio = ratio/Bratio
		SiteToGeneFactor = 1 + math.log10(ratioBratio)

		#Distribution difference between GoI and random background sets
		KLnum = TFinfo[7]
		#Bonferroni Correction for KL-Distance p-value
		correctKL = KLnum*float(number)
		#Zscore based on total DNA differences between GoI and Background
		Zscore = TFinfo[8]
		#P-value based on normal approximation to the binomial distribution (same approach as oPOSSUM)
		normPvalue = scipy.stats.norm.sf(Zscore)
		#Bonferroni Corrected p-value
		correctednormPvalue = normPvalue*float(number)
		#Fisher Exact test to compare putative TF regulated GoI genes vs Background genes 
		fisherExact = scipy.stats.chi2_contingency([[numGenes, BnumGenes], [(GoISize - numGenes), (BSize - BnumGenes)]], lambda_="log-likelihood")[1]
		#Bonferroni corrected p-value
		correctedFisherExact = fisherExact*float(number)
		
		#Create the Transcription Factor Overrepresentation Output file based on the listed criteria
		if KLnum < 0.25:
			outputString = '%s\t%d\t%d\t%d\t%.5g\t%d\t%d\t%.5g\t%.5g\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%s\n' % (TFname, Motif, numSites, numGenes, ratio, BnumSites, BnumGenes, Bratio, SiteToGeneFactor, KLnum, correctKL, Zscore, normPvalue, correctednormPvalue, fisherExact, correctedFisherExact, geneString)
			a.write(outputString)


	correction_value = len(Gene_CooccurList)

        for occurrence in sorted(Gene_CooccurList, key=lambda x: x[14], reverse = True):
                geneString = ", ".join(occurrence[13])
                outputString = "%s\t%s\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%s\n" % (occurrence[0], occurrence[1], occurrence[2], occurrence[3], occurrence[4], occurrence[5], occurrence[6], occurrence[7], occurrence[8], occurrence[9], occurrence[10], occurrence[11], occurrence[12], occurrence[14], geneString)
                b.write(outputString)


	for occurrence in sorted(Region_CooccurList, key=lambda x: x[14], reverse = True):
		regionList = []
		for region in occurrence[13]:
			string = str(region[0]) + ':' + str(region[1]) + '-' + str(region[2])
			regionList.append(string)
		regionString = ", ".join(regionList)
		outputString = "%s\t%s\t%d\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%s\n" % (occurrence[0], occurrence[1], occurrence[2], occurrence[3], occurrence[4], occurrence[5], occurrence[6], occurrence[7], occurrence[8], occurrence[9], occurrence[10], occurrence[11], occurrence[12], occurrence[14], regionString)	
		c.write(outputString)
	
	for occurrence in sorted(Proximity_CooccurList, key = lambda x: x[4]):
		pval_within50_corrected = float((math.pow(10, occurrence[4])) * correction_value)
		if pval_within50_corrected > 1:
                	pval_within50_corrected = 1
                pval_f_min_corrected = float((math.pow(10, occurrence[5])) * correction_value)
                if pval_f_min_corrected > 1:
                        pval_f_min_corrected = 1
                pval_f_mid_corrected = float((math.pow(10, occurrence[6])) * correction_value)
                if pval_f_mid_corrected > 1:
                        pval_f_mid_corrected = 1
                pval_f_far_corrected = float((math.pow(10, occurrence[7])) * correction_value)
                if pval_f_far_corrected > 1:
                        pval_f_far_corrected = 1
                pval_r_min_corrected = float((math.pow(10, occurrence[8])) * correction_value)
                if pval_r_min_corrected > 1:
                        pval_r_min_corrected = 1
                pval_r_mid_corrected = float((math.pow(10, occurrence[9])) * correction_value)
                if pval_r_mid_corrected > 1:
                        pval_r_mid_corrected = 1
                pval_r_far_corrected = float((math.pow(10, occurrence[10])) * correction_value)
                if pval_r_far_corrected > 1:
                        pval_r_far_corrected = 1
                outputString = "%s\t%s\t%d\t%d\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\t%.4g\n" % (occurrence[0], occurrence[1], occurrence[2], occurrence[3], pval_within50_corrected, pval_f_min_corrected, pval_f_mid_corrected, pval_f_far_corrected, pval_r_min_corrected, pval_r_mid_corrected, pval_r_far_corrected)
                d.write(outputString)

	print "Finished File Generation"

        a.close()
        b.close()
	c.close()
	d.close()

if __name__ == '__main__':
     ProcessCLI(sys.argv)


