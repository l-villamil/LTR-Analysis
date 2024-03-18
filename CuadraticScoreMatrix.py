
import numpy as np 
import dataReaderFASTASeondVersion as dr
import scoringMatrix as scores

def matchNumber (i_chain, j_chain, scoringMatrix): 
    #i y j son las cadenas i y j-ésima 

    matrix = np.zeros ((len(i_chain)+1,len(j_chain)+1))

    initialCost = 0
    matrix[0][0] = initialCost

    for row in range (len(i_chain)+1): 
        for col in range (len(j_chain)+1): 
            if (row == 0 and col == 0):
                pass
            elif (row == 0 and col != 0):
                matrix[row][col] =  matrix[row][col-1]+ cost(i_chain[row],'-',scoringMatrix)
            elif (col == 0 and row != 0): 
                matrix[row][col] = matrix [row-1][col]+cost('-',j_chain[col],scoringMatrix)
            else: 
                matrix[row][col] = max (
                                         matrix[row-1][col-1]+cost(i_chain[row-1],j_chain[col-1],scoringMatrix), 
                                         matrix[row][col-1]+cost(i_chain[row-1],'-',scoringMatrix), 
                                         matrix[row-1][col]+cost('-',j_chain[col-1],scoringMatrix)
                                         )

    return matrix




def cost (let1, let2, scoringMatrix):

    heading = ['A','C','G','T','-']
    let1Pos = heading.index(let1)
    let2Pos = heading.index(let2)
    return scoringMatrix[let1Pos][let2Pos]
    


#Backtracking para hallar nuevas cadenas con huecos 
def backTracking (matrix, i_chain, j_chain,scoringMatrix): 
    gapNumber=0

    i = len (matrix)-1
    j = len (matrix[0])-1
    matches = 0 
    lcs = "" #Longest common subsequence
    while j > 0 or  i > 0: 
        if i == 0 and j !=0: 
            i_chain =i_chain[:i]+"-"+i_chain[i:] 
            j-=1
        elif j ==0 and i !=0: 
            j_chain = j_chain[:j]+"-"+j_chain[j: ]
            gapNumber+=1
            i-=1
        elif  matrix[i][j]== matrix[i-1][j-1]+ cost(i_chain[i-1],j_chain[j-1],scoringMatrix): #I prefer to leave it with no spaces if is the same letter or assume the mismatch
            if cost(i_chain[i-1],j_chain[j-1],scoringMatrix) ==1:
                matches+=1
                lcs += i_chain[i-1]
            j-=1
            i-=1
        elif matrix[i][j] == matrix[i-1][j]+cost('-',j_chain[j-1],scoringMatrix):
            j_chain = j_chain[:j]+"-"+j_chain[j: ]
            gapNumber+=1
            i-=1
        elif matrix[i][j] == matrix [i][j-1]+cost(i_chain[i-1],'-',scoringMatrix): 
            i_chain =i_chain[:i]+"-"+i_chain[i:] 
            gapNumber+=1
            j-=1
            

    return matches, gapNumber,lcs[::-1],i_chain, j_chain
    

def pairwiseSequencing( i_chain, j_chain,scoringMatrix): 

    matchMatrix= matchNumber(i_chain,j_chain,scoringMatrix)
    #print(matchMatrix)

    newChains = backTracking(matchMatrix, i_chain, j_chain,scoringMatrix)    
    
    return newChains



match,mismatch,indel =1,0,0




"""i_chain ="TTTGACGATG"
j_chain = "TTCGAGGGATA"
a = pairwiseSequencing(i_chain,j_chain,scores.matrix(4,1,0,0))
print (a)
"""
"""genesDict = dr.getGenes("./data/exampleFLOGenes.fa")
maxChainLen = []
timesAlg1= []
spaceAlg1= []

timesAlg2= []
spaceAlg2= []

twoCompleteGenes = dr.getTwoRandom(genesDict)


minLen = min(len(twoCompleteGenes[0]), len(twoCompleteGenes[1]))
size = 10
twoFractionGenes = dr.getFractionFixedSize(twoCompleteGenes,size,size//10)
print(twoFractionGenes[0],twoFractionGenes[1])

"""


#print("The number of matches: {} \nThe number of gaps added: {}\nThe longest common subsequence {}\nThe chains: \n {}\n {}".format(result[0],result[1],result[2],result[3],result[4]))



"""
i_chain = cg.generateGenome(10)
j_chain = cg.generateGenome(10)

result = pairwiseSequencing(i_chain,j_chain,scores.matrix(4,1,-1,0))
print("The number of matches: {} \nThe number of gaps added: {}\nThe longest common subsequence {}\nThe chains: \n {}\n {}".format(result[0],result[1],result[2],result[3],result[4]))
"""
"""
def multipleIterations(iterations, maxLength): 
    fig, axs = plt.subplots(iterations)
    for e in range (iterations): 
        iLen = random.randint(0, maxLength)
        jLen = iLen+random.randint(0, 0)

        i_chain = cg.generateGenome(iLen)
        j_chain = cg.generateGenome(jLen)

        resultGraph = pairwiseSequencing(i_chain, j_chain)
        result = resultGraph[0]
        graph = result[1]

        print("The number of matches: {} \nThe number of gaps added: {}\nThe chains: \n {}\n {}".format(result[0],result[1],result[2],result[3]))

        #axs[e].imshow(graph, interpolation='nearest', cmap=cm.Greys_r)
    #fig.show()

        #print("The number of matches: {} \nThe number of gaps added: {}\nThe chains: \n {}\n {}".format(result[0],result[1],result[2],result[3]))
# First number is number of matches, second includes the new pair in sequences and last number is number og gaps added

#multipleIterations(1,10)
"""
