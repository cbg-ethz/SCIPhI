import sys
import re

def goLeft(childVec, pos):
    if childVec[pos][0] != -1:
        return childVec[pos][0]
    return pos
def goRight(childVec, pos):
    if childVec[pos][1] != -1:
        return childVec[pos][1]
    return pos
def goUp(parentVec, pos):
    if parentVec[pos] != -1:
        return parentVec[pos]
    return pos

def createSet(parentVec, childVec, nameVec):
    labels = [""] * len(parentVec)
    visits = [0] * len(parentVec)

    root = parentVec.index(-1)
    idx = root
    while (visits[root] != 3):
        visits[idx] += 1
        if visits[idx] == 1:
            idx = goLeft(childVec, idx)
        elif visits[idx] == 2:
            idx = goRight(childVec, idx)
        else:
            if nameVec[idx] != '':
                labels[idx] = labels[idx] + nameVec[idx]
            newIdx = goUp(parentVec, idx)
            if labels[newIdx] == '':
                labels[newIdx] = labels[idx]
            else:
                if newIdx != root:
                    labels[newIdx] = labels[newIdx] + "-" + labels[idx]
            idx = newIdx

    labelSet = set(())
    for label in labels:
        if len(label.split("-")) > 1:
            sortLabel = label.split("-")
            sortLabel.sort()
            labelSet.add('-'.join(sortLabel))
            #print("21",idx)
    return labelSet
        

inGTName = sys.argv[1]
inITName = sys.argv[2]

inGT = open(inGTName, 'r')
numNodes = 0
for line in inGT:
    if "->" in line:
        nodeId = int(re.split("\t | |-", line)[0]) + 1
        if (nodeId > numNodes):
            numNodes = nodeId

parentVecGT = [-1] * numNodes
nameVecGT = [""] * numNodes
childVecGT = [[-1,-1] for i in range(numNodes)] 

#print("parentVecGT: ", parentVecGT)
#print("childVecGT: ", childVecGT)

inGT.seek(0)
for line in inGT:
    if "->" in line:
        lineSplit = re.split("\t | |- |> |\[", line)
        #print(lineSplit)
        parentVecGT[int(lineSplit[0])] = int(lineSplit[2])
        if childVecGT[int(lineSplit[2])][0] == -1:
            childVecGT[int(lineSplit[2])][0] = int(lineSplit[0])
        else:
            childVecGT[int(lineSplit[2])][1] = int(lineSplit[0])
        sampleId = int(lineSplit[0]) - int((numNodes)/2)
        if sampleId >= 0:
            nameVecGT[int(lineSplit[0])] = str(sampleId)
#print(nameVecGT)
#print(parentVecGT)
#print(childVecGT)

setGT = createSet(parentVecGT, childVecGT, nameVecGT)



parentVecIT = [-1] * numNodes
nameVecIT = [""] * numNodes
childVecIT = [[-1,-1] for i in range(numNodes)] 

#print("parentVecIT: ", parentVecIT)
#print("childVecIT: ", childVecIT)

inIT = open(inITName, 'r')
for line in inIT:
    if "->" in line:
        lineSplit = re.split("-|>| ", line)
        parentVecIT[int(lineSplit[2])] = int(lineSplit[0])
        if childVecIT[int(lineSplit[0])][0] == -1:
            childVecIT[int(lineSplit[0])][0] = int(lineSplit[2])
        else:
            childVecIT[int(lineSplit[0])][1] = int(lineSplit[2])
    elif ".bam" in line:
        lineSplit = re.split("\[|\"|\.", line)
        nameVecIT[int(lineSplit[0])] = lineSplit[2]
#print(nameVecIT)
#print(parentVecIT)
#print(childVecIT)

setIT = createSet(parentVecIT, childVecIT, nameVecIT)

print(len(setGT ^ setIT))


