import sys

inFileName = sys.argv[1]
minCov = int(sys.argv[2])

inFile = open(inFileName, 'r')
for line in inFile:
    if line.startswith("#"):
        continue

    lineSplit = line.strip().split("\t")
    for value in lineSplit[9:-1]:
        if value != './.':
            cellSplit = value.split(":")
            cov = int(cellSplit[2])
            alt = int(cellSplit[1].split(",")[1])
            if cov >= minCov:
                print(alt/cov)
