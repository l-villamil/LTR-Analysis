import localAlignment as la


# Filtra la DB y selecciona los genes con tamaño mayor o igual minlength
def FIlterGenByLen(genes,minlength):
    filteredGenes={}
    for gen in genes:
        if minlength<=len(genes[gen]):
           filteredGenes[gen]=genes[gen]
    return filteredGenes


#Toma los extremos de una secuencia
def getExtremeSequence(gen,window=1000):
    return [gen[0:window], gen[-window:]]


# Toma los extremos de todos los genes y los devuelve en un diccionario
def multipleExtremes(genes):
    extremes = {}
    for t in genes: 
        extremes[t] = getExtremeSequence(genes[t])
    return extremes



# Retorna la cantidad de matches exactos en los dos extremos

def getTotalPercentageMatches(extreme1, extreme2):
    realMatch = 0
    repeatedLenght= len(extreme1)
    for i in range (0, repeatedLenght): 
        if extreme1[i] == extreme2[i]:
            realMatch += 1
    return [(realMatch/repeatedLenght)*100,(1-(realMatch/repeatedLenght))*100]



def getMatches(extremes):#Extrems bein the dictionary with two extremes of each LTR

    #dataPoints = np.linspace(0,len(extremes),len(extremes))

    dictMatches = {}
    matchPercentage = []
    mismatchPercentage = []
    for name in extremes: 
        extreme1 = extremes[name][0]
        extreme2= extremes[name][1]

        realMatchPercentage, realUnmatchPercentage = getTotalPercentageMatches(extreme1, extreme2)
        if realMatchPercentage>80: 
            print(name)
        

        """extremes[name].append(realMatch)
        extremes[name].append(repeatedLenght-realMatch)"""
        dictMatches[name] = [realMatchPercentage, realUnmatchPercentage]
        matchPercentage.append(realMatchPercentage)
        mismatchPercentage.append(realUnmatchPercentage)

    return dictMatches


def getRealMatchesUnmatchesList(dict):

    realMatch = []
    realUnmatch = []
    for name in dict: 
        realMatch.append(dict[name][0])
        realUnmatch.append(dict[name][1])
    return realMatch, realUnmatch




def getAlignedLocal(extremes):
    alignedMatches = []
    ltrSize=[]
    coords=[]
    i=1
    results="results.txt"
    try:

        with open(results,'w') as archivo:
            

            for t in extremes: 
                seq1 = extremes[t][0]
                seq2 =extremes[t][1]
                localResult=la.main(seq1,seq2)
                origin_coords=localResult[4]
                max_coords=localResult[3]
                sim_tup = localResult[2]
                seq2_new=localResult[1]
                seq1_new=localResult[0]
                sim_tup=sim_tup[1]
                # print(max_coords,origin_coords,sim_tup)

                coords.append([max_coords,origin_coords])
                alignedMatches.append(sim_tup)
                ltrSize.append(len(seq1_new))
                archivo.write(''.join(seq1_new)+ '\n')
                archivo.write(''.join(seq2_new)+ '\n')
                archivo.write(str(sim_tup)+ '\n')
                archivo.write(str(max_coords)+ '\n')
                archivo.write(str(origin_coords)+ '\n')



                # if(i>10):
                #    break
                
                print(i)
                i+=1
    except Exception as e:
        print("Error al escribir el archivo:", str(e))
    return alignedMatches,ltrSize,coords