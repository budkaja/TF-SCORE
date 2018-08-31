import sys
import math
import sys
import itertools
import os.path
import string
import json
import MOODS

def matrix_mapper_long(matrix, matrix_name, threshold, outputDirectory, sequences):

	matchingCount = 0
	maxScore = 0
	scoreList = []

	motifLength = len(matrix[0])

	for position in range(0, motifLength):
		Avalue = float(matrix[0][position])
		Cvalue = float(matrix[1][position])
		Gvalue = float(matrix[2][position])
		Tvalue = float(matrix[3][position])
		ratioList = [Avalue, Cvalue, Gvalue, Tvalue]
		score = max(ratioList)
		scoreList.append(score)
		maxScore += score

        name_of_file = matrix_name
        completeName = os.path.join(outputDirectory, name_of_file)
        outputFile = open(completeName, "w")

        for sequenceName in sequences.keys():
                sequence = sequences[sequenceName]
                chromosome, start, stop, strand = sequenceName.replace(':', ' ').replace('-', ' ').split()
                for fragmentStart in xrange(0, (len(sequence) - motifLength + 1)):
                        fragment = sequence[fragmentStart : fragmentStart + motifLength]
			fragmentScore = 0
			bestRemaining = sum(scoreList)
			for position in range(0, motifLength):
				nucleotide = str(fragment[position])
				currentScore = 0
				if nucleotide not in ['A', 'C', 'G', 'T']:
					continue
				else:
					if nucleotide == 'A':
						currentScore = float(matrix[0][position])
						fragmentScore = fragmentScore + currentScore
					elif nucleotide == 'C':
						currentScore = float(matrix[1][position])
						fragmentScore = fragmentScore + currentScore
					elif nucleotide == 'G':
						currentScore = float(matrix[2][position])
						fragmentScore = fragmentScore + currentScore
					elif nucleotide == 'T':
						currentScore = float(matrix[3][position])
						fragmentScore = fragmentScore + currentScore
					bestRemaining = float(fragmentScore + sum(scoreList[position+1:])) 
					if bestRemaining < threshold:
						break
					else:
						pass			

			if fragmentScore < threshold:
				continue

			else:
				if strand == 'p':
                                	outputString = "%s\t%d\t%d\t%s\t%s\t%f\n" % (chromosome, (int(start)+int(fragmentStart)-1), (int(start)+int(fragmentStart)+int(motifLength)-1), fragment, strand, fragmentScore)
                                else:
                                	outputString = "%s\t%d\t%d\t%s\t%s\t%f\n" % (chromosome, (int(stop)-(int(fragmentStart)+int(motifLength))), (int(stop)-int(fragmentStart)), fragment, strand, fragmentScore)

                                outputFile.write(outputString)
        print "Finished: ", matrix_name

def matrix_mapper(matrix, matrix_name, threshold, outputDirectory, sequences):
        matchingCount = 0
	maxScore = 0
        scoreList = []
        
	A = matrix[0]
	C = matrix[1]
	G = matrix[2]
	T = matrix[3]

	motifLength = len(A)
	for position in range(0, motifLength):
                Avalue = float(matrix[0][position])
	        Cvalue = float(matrix[1][position])
              	Gvalue = float(matrix[2][position])
                Tvalue = float(matrix[3][position])
		ratioList = [Avalue, Cvalue, Gvalue, Tvalue]
	        score = max(ratioList)
        	scoreList.append(score)
	        maxScore += score

        maxScoreHolder = [maxScore]
	ScorerList = maxScoreHolder + scoreList
	
	List = list(itertools.product('ACGT', repeat = 5))
	newList = []

	for item in List:
        	newList.append(''.join(item))

        diff = maxScore - threshold
	count = 0

        storePassedSegments = {}
	ranges = {}
        for x in range(0, (motifLength/5)):

		storePassedSegments[x] = {}
                startpos = x * 5
                rangeScores = ScorerList[(startpos + 1):(startpos + 6)]
	        rangeMaxScore = sum(rangeScores)
                ranges[x] = rangeMaxScore
	        nicecount = 0
                badcount = 0
        	for seq in newList:
	        	score = rangeMaxScore
                	for y in range(0,5):
                        	seqdiff = 0
        	        	currentPos = startpos + y
                                currentNuc = seq[y]
                    		if str(currentNuc) == 'A':
                        		currentScore = float(matrix[0][currentPos])
                    		elif str(currentNuc) == 'C':
        	        		currentScore = float(matrix[1][currentPos])
	            		elif str(currentNuc) == 'G':
                       			currentScore = float(matrix[2][currentPos])
                    		elif str(currentNuc) == 'T':
                        		currentScore = float(matrix[3][currentPos])

		                seqdiff = rangeScores[y] - currentScore

        			score = score - seqdiff

				if (rangeMaxScore - score) > diff:
                        		badcount = badcount + 1
		                        break
            		if (rangeMaxScore - score) > diff:
        			badcount = badcount + 1
				pass
	    		else:
                		nicecount = nicecount + 1
                		storePassedSegments[x][seq] = score
                
        	count = count + 5

	remaining = motifLength - count
        rangeScores = ScorerList[(count+1):(count+remaining+1)]
        rangeMaxScore = sum(rangeScores)

        shortList = list(itertools.product('ACGT', repeat = remaining))
        newShortList = []
        
	if remaining == 0:
        	pass 
        else:
                ranges[(count/5)] = rangeMaxScore
	        storePassedSegments[(count/5)] = {}
	        for item in shortList:
                	newShortList.append(''.join(item))
	        for sequence in newShortList:
                	score = rangeMaxScore
	        	for num in range(0, remaining):
                        	seqdiff = 0
        	    		currentPos = count + num
	            		currentNuc = sequence[num]
                    		if str(currentNuc) == 'A':
                        		currentScore = float(A[currentPos])
                    		elif str(currentNuc) == 'C':
        	        		currentScore = float(matrix[1][currentPos])
	                        elif str(currentNuc) == 'G':
                        		currentScore = float(G[currentPos])
                    		elif str(currentNuc) == 'T':
                        		currentScore = float(T[currentPos])

				seqdiff = rangeScores[num] - currentScore
	            		score = score - seqdiff
                    		if (rangeMaxScore - score) > diff:
                        		break
                	if (rangeMaxScore - score) > diff:
                    		pass
	        	else:
                    		storePassedSegments[(count/5)][sequence] = score

		count = count + remaining

	fragmentNameList = []			        
	for fragmentName, seqdict in storePassedSegments.iteritems():
                fragmentNameList.append(fragmentName)

        largest = max(fragmentNameList)

	countera = 0
        counterb = 1
	scorediff = 0
	while countera < largest:
               	smallseqdict = {}
                for sequence1, score1 in storePassedSegments[countera].iteritems():
	                for sequence2, score2 in storePassedSegments[counterb].iteritems():
                               	first = ranges[countera]
                       	        second = ranges[counterb]
               	                combinedMaxScore = (float(first) + float(second))
                                currentseq = sequence1 + sequence2
	                        currentscore = score1 + score2
				scorediff = combinedMaxScore - currentscore
                       	        if scorediff > diff:
               	                        pass
                                else:
	                                smallseqdict[currentseq] = float(currentscore)
               	del storePassedSegments[countera]
                first1 = ranges[countera]
	        second1 = ranges[counterb]
               	del ranges[countera]
                ranges[counterb] = float(first1) + float(second1)
	        storePassedSegments[counterb] = smallseqdict
               	countera = countera + 1
                counterb = counterb + 1
	kmerDict = storePassedSegments[countera]
	name_of_file = matrix_name
        completeName = os.path.join(outputDirectory, name_of_file)
        outputFile = open(completeName, "w")
        for sequenceName in sequences.keys():
        	sequence = sequences[sequenceName]
                chromosome, start, stop, strand = sequenceName.replace(':', ' ').replace('-', ' ').split()
                for fragmentStart in range(0, len(sequence) - motifLength + 1):
                	fragment = sequence[fragmentStart : fragmentStart + motifLength]
                    	if fragment in kmerDict:
                        	if strand == 'p':
                        		outputString = "%s\t%d\t%d\t%s\t%s\t%f\n" % (chromosome, (int(start)+int(fragmentStart)-1), (int(start)+int(fragmentStart)+int(motifLength)-1), fragment, strand, kmerDict[fragment])
                        	else:
		                       	outputString = "%s\t%d\t%d\t%s\t%s\t%f\n" % (chromosome, (int(stop)-(int(fragmentStart)+int(motifLength))), (int(stop)-int(fragmentStart)), fragment, strand, kmerDict[fragment])

		                outputFile.write(outputString)
	print "Finished: ", matrix_name

def reverse_complement(sequencea):
	reverse = sequencea[::-1]

	baseComplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
	bases = list(reverse)
	complement = ''.join([baseComplement[base] for base in bases])

	return complement

def ProcessCLI(args):
    
    outputDirectory = '/N/u/jubudka/Mason/BindingFiles/'
    weightMatrixDirectory = '/N/u/jubudka/Mason/PWMsmall/'
    sequencesFileName = 'FASTA_All_Merged_Encode.fasta'
    p_val = 0.0001

    print args
    for i in xrange(len(args)):
        if args[i] == "-f":
	    sequencesFileName = args[i+1]
	    print "Fasta file is: ", sequencesFileName
	elif args[i] == "-p":
	    weightMatrixDirectory = args[i+1]
	    print "PWM file is: ", weightMatrixDirectory
	elif args[i] == "-o":
	    outputDirectory = args[i+1]
	    print "Output file is: ", outputDirectory
	elif args[i] == "-t":
	    p_val = float(args[i+1])

    if not os.path.exists(outputDirectory):
	os.makedirs(outputDirectory)
	

    # file for saving average score stuff
    # load position weight matrices
    # order is A C G T
    sequences = {}
    seqIDs = []
    current_sequence = ''
    sequencesFile = open(sequencesFileName)

    aCount = 0
    cCount = 0
    gCount = 0
    tCount = 0
    totalLength = 0

    for lines in sequencesFile:
	line = lines.strip()
	if line == '':
	    continue
	if (line[0].startswith('>')):
	    seqIDs.append(line[1:])
	    #add previous sequence to dictionary
	    #create the reverse complement and add to dictionary
	    #perform nucleotide counting
	    #reset sequence to '' for next fasta sequence
	    if (len(current_sequence) > 0):
		upper_current_sequence = current_sequence.upper()
		seqID = seqIDs.pop(0)
		sequences[seqID + ' ' + 'p'] = upper_current_sequence
		reverseSequence = reverse_complement(upper_current_sequence)
		sequences[seqID + ' ' + 'm'] = reverseSequence
		aCount = aCount + upper_current_sequence.count('A')
                cCount = cCount + upper_current_sequence.count('C')
                gCount = gCount + upper_current_sequence.count('G')
                tCount = tCount + upper_current_sequence.count('T')
                totalLength = totalLength + len(current_sequence)

	    current_sequence = ''
	else:
	    current_sequence += line

    upper_current_sequence = current_sequence.upper()
    seqID = seqIDs.pop(0)
    sequences[seqID + ' ' + 'p'] = upper_current_sequence
    reverseSequence = reverse_complement(upper_current_sequence)
    sequences[seqID + ' ' + 'm'] = reverseSequence

    aCount = aCount + upper_current_sequence.count('A')
    cCount = cCount + upper_current_sequence.count('C')
    gCount = gCount + upper_current_sequence.count('G')
    tCount = tCount + upper_current_sequence.count('T')
    totalLength = totalLength + len(current_sequence)

    aContent = aCount/float(totalLength)
    cContent = cCount/float(totalLength)
    gContent = gCount/float(totalLength)
    tContent = tCount/float(totalLength)
  
    backgroundScores = {'A':aContent, 'C':cContent, 'G':gContent, 'T':tContent}
    bg = [backgroundScores['A'], backgroundScores['C'], backgroundScores['G'], backgroundScores['T']]
    print bg

    matrix_names = [filename for filename in os.listdir(weightMatrixDirectory) if filename[-4:] == '.pfm']
    pseudocount = 1

    matrices = [MOODS.load_matrix(weightMatrixDirectory + filename) for filename in matrix_names]

    matrices = [MOODS.count_log_odds(matrix, bg, pseudocount) for matrix in matrices]

    thresholds = [MOODS.threshold_from_p(matrix, bg, p_val) for matrix in matrices]


    for (matrix, matrix_name, threshold) in zip(matrices, matrix_names, thresholds):
	    motifLength = len(matrix[0])
	    if motifLength >= 18:
		matrix_mapper_long(matrix, matrix_name, threshold, outputDirectory, sequences)
		continue
	    else:
		matrix_mapper(matrix, matrix_name, threshold, outputDirectory, sequences)

    print "Finished"		

if __name__ == '__main__':
        ProcessCLI(sys.argv)
