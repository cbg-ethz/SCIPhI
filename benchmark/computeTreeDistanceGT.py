import sys
import math
import random

random.seed(42)

inFileNameTruth = sys.argv[1]
inFileNameInf = sys.argv[2]

def computeCellDist(values, gts):
    numCells = len(gts)
    for i in range(numCells):
        for j in range(i + 1, numCells):
            values[i][j] += abs(gts[i] - gts[j])

def normalize(values, numMuts):
    numCells = len(values)
    for i in range(numCells):
        for j in range(i + 1, numCells):
            values[i][j] = float(values[i][j])/float(numMuts)

def computeDist(valuesTruth, valuesInf):
    numCells = len(valuesTruth)
    result = 0
    for i in range(numCells):
        for j in range(i + 1, numCells):
            result += pow(valuesTruth[i][j] - valuesInf[i][j], 2)
    return math.sqrt(float(result/((numCells - 1) * (numCells) / 2)))
    #return (float(result/((numCells - 1) * (numCells) / 2)))

inFileTruth = open(inFileNameTruth, 'r')
numCells = len(inFileTruth.readline().strip().split("\t")) -1
truthValues = [[0 for i in range(numCells)] for j in range(numCells)]
numMuts = 0
for line in inFileTruth:
    numMuts += 1
    lineSplit = line.strip().split("\t")[1:]
    gts = list(map(lambda x:int(x != '0'), lineSplit))
    computeCellDist(truthValues, gts)
normalize(truthValues, numMuts)


#for i in range(len(truthValues)):
#    print(truthValues[i])

def shuffleGT(gts, randomOrder):
    temp = gts
    for i in range(len(temp)):
        temp[i] = gts[randomOrder[i]]
    return temp

inFileInf = open(inFileNameInf, 'r')
inFileInf.readline()
infValues = [[0 for i in range(numCells)] for j in range(numCells)]
numMuts = 0
randomOrder = [i for i in range(numCells)]
#print(randomOrder)
random.shuffle(randomOrder)
#print(randomOrder)
for line in inFileInf:
    numMuts += 1
    #lineSplit = line.strip().split("\t")[2:]
    #gts = list(map(lambda x: float(x.split("|")[2]), lineSplit))
    lineSplit = line.strip().split("\t")[1:]
    gts = list(map(lambda x: int(x != '0'), lineSplit))
    gts = shuffleGT(gts, randomOrder)
    computeCellDist(infValues, gts)
normalize(infValues, numMuts)

#for i in range(len(infValues)):
#    print(infValues[i])

print(computeDist(truthValues, infValues))
