import sys
import math
import random

random.seed(42)

inFileNameTruth = sys.argv[1]
inFileNameInf = sys.argv[2]
randomize = sys.argv[3]

def computeCellDist(cell1, cell2):
    dist = []
    for i in range(len(cell1)):
        dist.append(abs(cell1[i] - cell2[i]))
    return dist

def computeSquareDist(dist1, dist2):
    dist = 0
    for i in range(len(dist1)):
        dist += pow(dist1[i] - dist2[i], 2)
    return float(dist / len(dist1))

def computeTreeDist(distValues):
    numCells = len(distValues)
    dist = 0
    for i in range(numCells):
        for j in range(i + 1, numCells):
            dist += distValues[i][j]
    return (float(math.sqrt(dist/((numCells - 1) * (numCells) / 2))))

def shuffleGT(gts, randomOrder):
    temp = gts
    for i in range(len(temp)):
        temp[i] = gts[randomOrder[i]]
    return temp

inFileInf = open(inFileNameInf, 'r')
numCells = len(inFileInf.readline().strip().split("\t")) - 2
if(randomize == '1'):
    numCells += 1

infValues = [[0 for i in range(numCells)] for j in range(numCells)]
infGT = [ [] for i in range(numCells)]
infDict = {}
numMuts = 0
randomOrder = [i for i in range(numCells)]
random.shuffle(randomOrder)
for line in inFileInf:
    numMuts += 1

    if(randomize == '1'):
        lineSplit = line.strip().split("\t")
        lineSplit[0] = str(numMuts)
        gts = list(map(lambda x:float(x != '0'), lineSplit[1:]))
        gts = shuffleGT(gts, randomOrder)

    else:
        lineSplit = line.strip().split("\t")[1:]
        gts = list(map(lambda x: float(x.split("|")[2]), lineSplit[1:]))
    
    infDict[lineSplit[0]] = gts

for key, value in sorted(infDict.items()):
    for i in range(numCells):
        infGT[i].append(value[i])

inFileTruth = open(inFileNameTruth, 'r')
truthGT = [ [] for i in range(numCells)]
truthValues = [[0 for i in range(numCells)] for j in range(numCells)]
truthDict = {}
numMuts = 0
for line in inFileTruth:
    numMuts += 1
    lineSplit = line.strip().split("\t")
    if(randomize == '1'):
        lineSplit[0] = str(numMuts)
    gts = list(map(lambda x:float(x != '0'), lineSplit[1:]))
    if lineSplit[0] in infDict:
        truthDict[lineSplit[0]] = gts     
        
for key, value in sorted(truthDict.items()):
    for i in range(numCells):
        truthGT[i].append(value[i])

distValues = [[0 for i in range(numCells)] for j in range(numCells)]
for i in range(numCells):
    for j in range(i + 1, numCells):
        dist1 = computeCellDist(truthGT[i], truthGT[j])
        dist2 = computeCellDist(infGT[i], infGT[j])
        distValues[i][j] = computeSquareDist(dist1, dist2 )

print(computeTreeDist(distValues))
