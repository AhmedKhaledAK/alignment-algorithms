Q = "AGTCAT"
R = "AAGTATCGA"

symbols = ['A', 'C', 'G', 'T']


def getKMers(Q, k):
    kmers = [Q[i : i + k] for i in range(len(Q) - k + 1)]
    print(kmers)
    return kmers


def findSynonyms(kmers):
    synons = []
    for i in range(len(kmers)):
        kmer = kmers[i]
        synons.append((kmer, i))
        kmer = list(kmer)
        for j in range(len(kmer)):
            orgChar = kmer[j]
            
            ######
            for k in range(len(symbols)):
                if symbols[k] != orgChar:
                    kmer[j] = symbols[k]
                    synons.append(("".join(kmer), i))
                
            kmer[j] = orgChar
    return synons            


def calKey(synons, k):
    for i in range(len(synons)):
        kmer = synons[i][0]
        key = 0
        for j in range(len(kmer)):
            if kmer[j] == 'C':
                key += (4 ** (k-j-1)) * 1
            elif kmer[j] == 'G':
                key += (4 ** (k-j-1)) * 2
            elif kmer[j] == 'T':
                key += (4 ** (k-j-1)) * 3

        synons[i] = (kmer, synons[i][1], key)

    synons = sorted(synons, key=lambda x: x[2])
    return synons

def main():
    kmers = getKMers(Q, 3)
    synons = findSynonyms(kmers)
    synons = calKey(synons, 3)
    print(synons)

    ### DB sequence
    print("-----")
    kmersDB = getKMers(R, 3)
    for i in range(len(kmersDB)):
        kmersDB[i] = (kmersDB[i], i)
        
    kmersDB = calKey(kmersDB, 3)
    print(kmersDB)


main()