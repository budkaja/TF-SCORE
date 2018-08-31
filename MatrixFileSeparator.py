import sys
import itertools
import os.path
import string

weightMatricesFileName = '/N/u/jubudka/Mason/JASPARhumanPWM.txt'
 
weightMatrices = {}

weightMatricesFile = open(weightMatricesFileName)

weightMatricesData = filter(None, weightMatricesFile.read().split(">"))

for matrix in weightMatricesData:
	matrixData = matrix.strip().split("\r\n")
	weightMatrices[matrixData[0]] = [x.split("\t") for x in matrixData[1:]]

samplelist = open('SampleList.txt', "w")

for MatrixName, MatrixList in weightMatrices.iteritems():
	name = MatrixName.split(" ")
	Name = "_".join(name)+".matrix"
	print Name

	samplelist.write(Name + '\n')
	
	a = "\t".join(MatrixList[0])
	c = "\t".join(MatrixList[1])
	g = "\t".join(MatrixList[2])
	t = "\t".join(MatrixList[3])

	file = open(Name, "w")
	file.write('A' + '\t' + '|' + '\t' + a + '\n')
        file.write('C' + '\t' + '|' + '\t' + c + '\n')
        file.write('G' + '\t' + '|' + '\t' + g + '\n')
        file.write('T' + '\t' + '|' + '\t' + t)
	
	file.close
samplelist.close
