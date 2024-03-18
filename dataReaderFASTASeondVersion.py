
import random

def getGenes (path:str): 
    f = open(path, "r")

    genes = {}
    i = 0 

    line = f.readline()
    while line!=None and len(line)!=0 and len(genes)<20000:
        if line[0] == ">": 
            name = line[1:]
            gene_i = ""
            line = f.readline()

            while len(line)!=0 and line[0] != ">": 
                gene_i+= line.strip("\n")
                line = f.readline()
                
            genes[name] = gene_i 
            
        
    return genes




def getLTRTransposons(transposons): 
    ltrs = {}
    for tName in transposons: 
        type = tName.split("#")
        if type[1][:3] == "LTR": 
            ltrs[tName] = transposons[tName]
    return ltrs

def getNoLTRTransposons(transposons): 
    ltrs = {}
    for tName in transposons: 
        type = tName.split("#")

        if type[1][:3]  == "LTR":  
            ltrs[tName] = transposons[tName]
    return ltrs




def getTwoRandom(genes:dict): 
    
    keys = list(genes.keys())
    randKey1 = random.choice(keys)
    gene1 = genes[randKey1]
    #To obtain different genes
    del genes[randKey1]
    keys = list(genes.keys())

    gene2 = genes[random.choice(keys)]

    return gene1,gene2


def getFraction(twoGenes, minSize, maxSize,maxdiffSize): 
    #Important: MaxSize must be  the min size (avoid problems)
    gene1 = twoGenes[0]
    gene2 =twoGenes[1]

    size1 = random.randint(minSize,maxSize)
    #Get random size of size 2 respect to a maz dif passed by parameter

    if (random.randint(0,1) == 0): 
        size2 = size1-random.randint(0,maxdiffSize)
    else: 
        size2 = size1+random.randint(0,maxdiffSize)



    if len(gene1)<len(gene2):
       start = random.randint(0,len(gene1)-size1)
    else:
       start = random.randint(0,len(gene2)-size1)

    if len(gene1)<len(gene2):
        result = gene1[start:start+size1],gene2[start:start+size2]
    else: 
        result = gene1[start:start+size2],gene2[start:start+size1]


       
    return result

def getFractionFixedSize(twoGenes, size, maxdiffSize): 
    gene1 = twoGenes[0]
    gene2 =twoGenes[1]

    
    if (random.randint(0,1) == 0): 
        size2 = size-random.randint(0,maxdiffSize)
    else: 
        size2 = size+random.randint(0,maxdiffSize)
    size = int(size)
    size2= int(size2)

    if len(gene1)<len(gene2):
       start = random.randint(0,len(gene1)-size)
    else:
       start = random.randint(0,len(gene2)-size)

    if len(gene1)<len(gene2):
        result = gene1[start:start+size],gene2[start:start+size2]
    else: 
        result = gene1[start:start+size2],gene2[start:start+size]

    return result

    









