import sys

inFileName = sys.argv[1]

inFile = open(inFileName, 'r')
inFile.readline()

def getXCord(pix, x):
    return (x + len(pix[0])) % len(pix[0])

def getYCord(pix, y):
    return (y + len(pix)) % len(pix)

def computeDist(pix):
    dist = 0.0
    for x in range(len(pix[0])):
        for y in range(len(pix)):
            lDist = 0
            numN = -1
            for xPos in range(x-1,x+2):
                for yPos in range(y-1,y+2):
                    if pix[y][x] != 'NaN':
                        neighbor = pix[getYCord(pix,yPos)][getXCord(pix,xPos)]
                        if neighbor != 'NaN':
                            lDist += abs(float(pix[y][x]) - float(neighbor))
                            numN += 1
            if numN > 0:
                dist += lDist/numN
    return dist

pix = []
for line in inFile:
    lineSplit = line.strip().split("\t")[1:]
    pix.append(lineSplit)

print(computeDist(pix) / (len(pix)*len(pix[0])))
#def computeDist(pix):
#    dist = 0.0
#    for x in range(len(pix[0])):
#        for y in range(len(pix)):
#            for xPos in range(x-1,x+2):
#                for yPos in range(y-1,y+2):
#                    #print(x, y, xPos, yPos, pix[getYCord(pix,yPos)][getXCord(pix,xPos)])
#                    dist += abs(pix[y][x] - pix[getYCord(pix,yPos)][getXCord(pix,xPos)])
#    return dist
#
#pix = []
#for line in inFile:
#    lineSplit = line.strip().split("\t")[1:]
#    #for i in range(len(lineSplit)):
#    #    if lineSplit[i] == 'NaN':
#    #        lineSplit[i] = '0.5'
#    pix.append([float(value) for value in lineSplit])
#
#print(computeDist(pix) / (8 * len(pix)*len(pix[0])))
