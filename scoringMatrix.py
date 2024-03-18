import numpy as np


def matrix(numStrings, match, mismatch, indel):
    matrix = np.zeros((numStrings+1,numStrings+1))
    inum = 0
    jnum = 0
    for inum in range(len(matrix)): 
        row = matrix[inum]
        for jnum in range (len(row)): 
            if inum == jnum: 
                row[jnum]= match
            elif inum == numStrings or jnum== numStrings: 
                row[jnum] = indel
            else:
                row[jnum] = mismatch
            jnum+=1
        inum+=1
    return matrix
