Q = "GACAGC"
R = "ACGGATTCCATAT"

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


def ungappedExtend(qKmer, dbKmer, qkmers, k):
    qWord = qkmers[qKmer[1]]
    print(qWord)

    dbWord = dbKmer[0]

    startq = qKmer[1]
    startdb = dbKmer[1]
    threshold = 1
    score = 0


    startidx_q = startq
    lastidx_q = startq + k
    startidx_db = startdb
    lastidx_db = startdb + k

    ##seed score 
    i = 0
    while i < len(qWord):
        if qWord[i] == dbWord[i]:
            score += 1
        else:
            score -= 1
        i += 1


    ##right extension

    i, j = startq + k, startdb + k

    while i < len(Q) and j < len(R) and score >= threshold:
        if Q[i] == R[j]:
            score += 1
        else:
            score -= 1
        i += 1
        j += 1

    lastidx_q = i - 1
    lastidx_db = j - 1


    ##left extension
    
    i, j = startq - 1, startdb - 1

    while i >= 0 and j >= 0 and score > threshold:
        if Q[i] == R[j]:
            score += 1
        else:
            score -= 1
        i -= 1
        j -= 1  

    if i < 0 or j < 0:
        i += 1
        j += 1

    startidx_q = i
    startidx_db = j

    return (Q[startidx_q : lastidx_q + 1] , R[startidx_db : lastidx_db + 1])

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


    print("----")
    for kmer in kmersDB:
        print(kmer)
        output = [item for item in synons if item[2] == kmer[2]]
        print(output)
        if len(output) == 1:
            alignment = ungappedExtend(output[0], kmer, kmers, 3)
            print("alginment: ")
            print(alignment)
        print("----")

main()



'''
def ungapped_extension(query_seq, db_seq, start_q, start_db, k):
  score = 0
  threshold = 1

  startidx_q = start_q
  startidx_db = start_db
  endidx_q = start_q + k
  endidx_db = start_db + k

  print("\n-----UNGAPPED EXTENSION------")
  print("Start query index: ", start_q)
  print("Start db index: ", start_db)


  # calculate score for length of seed
  for i in range(k):
    if query_seq[start_q+i] == db_seq[start_db+i]:
      score = score + 1
    else:
      score = score - 1
  print("seed score: ", score)


  # extend to the right
  query_idx = start_q + k
  db_idx = start_db + k
  
  while(1):
    # if out of range, break
    if query_idx >= len(query_seq) or db_idx >= len(db_seq):
      break

    # if match or mismatch
    if query_seq[query_idx] == db_seq[db_idx]:
      score = score + 1
    else:
      score = score - 1
    
    if score < threshold:
      break

    query_idx = query_idx + 1
    db_idx = db_idx + 1

  endidx_db = db_idx - 1
  endidx_q = query_idx - 1

  print("score after aligning to the right: ", score)


  # extend to the left
  query_idx = start_q - 1
  db_idx = start_db - 1
  print("qidx: ", query_idx)
  print("didx: ", db_idx)
  
  while(1):
    # query_idx = query_idx - 1
    # db_idx = db_idx - 1

    if query_idx < 0:
      query_idx = 0
      db_idx = db_idx + 1
      break

    if db_idx < 0:
      db_idx = 0
      query_idx = query_idx + 1
      break

    if query_seq[query_idx] == db_seq[db_idx]:
      score = score + 1
    else:
      score = score - 1
    
    if score < threshold:
      print("under threshold!")
      query_idx = query_idx + 1
      db_idx = db_idx + 1
      break
    
    query_idx = query_idx - 1
    db_idx = db_idx - 1

  startidx_db = max(0, db_idx)
  startidx_q = max(0, query_idx)

  print("Query sequence: ", query_seq[startidx_q:endidx_q+1],startidx_q, endidx_q)
  print("DB sequence: ", db_seq[startidx_db:endidx_db+1], startidx_db, endidx_db)
'''