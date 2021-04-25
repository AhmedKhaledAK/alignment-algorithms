
query_seq = 'GACAGC'
db_seq = 'ACGGATTCCATAT'



def getSubMatrix(match, mismatch, seq1, seq2):
    H, W = len(seq1), len(seq2)
    matrix = [[0 for i in range(W)] for y in range(H)]
    
    i = 0
    while i < H:
        j = 0
        while j < W:
            if seq1[i] == seq2[j]:
                matrix[i][j] = match
            else:
                matrix[i][j] = mismatch
            j += 1
        i += 1

    return matrix


def getAlignment(scoringMatrix, pointerMatrix, startPoint, seq1, seq2):
    alignments = []

    for sp in range(len(startPoint)):
        str1 = str2 = ""
        i, j = startPoint[sp][0], startPoint[sp][1]
        while scoringMatrix[i][j] != 0:
            if pointerMatrix[i][j][2] == 0:
                str1 += seq1[i - 1]
                str2 += seq2[j - 1]
            elif pointerMatrix[i][j][2] == 1:
                str1 += seq1[i - 1]
                str2 += "-"
            elif pointerMatrix[i][j][2] == 2:
                str1 += "-"
                str2 += seq2[j - 1]
            hold = i
            i = pointerMatrix[i][j][0]
            j = pointerMatrix[hold][j][1]
        alignments.append((str1[::-1], str2[::-1]))
    return alignments



def alignUsingSW(seq1, seq2, match=1, mismatch=-1, gapPenalty=1):
    subMatrix = getSubMatrix(match, mismatch, seq1, seq2)

    H, W = len(seq1) + 1, len(seq2) + 1
    scoringMatrix = [[0 for i in range(W)] for y in range(H)]
    pointerMatrix = [[(-1,-1,-1, -1) for i in range(W)] for y in range(H)]

    maxScore = -1
    startPoint = []

    i = 0
    while i < H:
        j = 0
        while j < W:
            if i == 0 or j == 0:
                scoringMatrix[i][j] = 0
            else:
                scoringMatrix[i][j] = max(0, scoringMatrix[i - 1][j - 1] + subMatrix[i - 1][j - 1], scoringMatrix[i - 1][j] - gapPenalty, scoringMatrix[i][j - 1] - gapPenalty)
                
                if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + subMatrix[i - 1][j - 1]:
                    pointerMatrix[i][j] = (i - 1, j - 1, 0, scoringMatrix[i][j])
                elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] - gapPenalty:
                    pointerMatrix[i][j] = (i - 1, j, 1, scoringMatrix[i][j])
                elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] - gapPenalty:
                    pointerMatrix[i][j] = (i, j - 1, 2, scoringMatrix[i][j])

            if maxScore < scoringMatrix[i][j]:
                maxScore = scoringMatrix[i][j]
                #startPoint = (i, j)
            
            j += 1
        i += 1


    i = 0
    while i < H:
        j = 0
        while j < W:
            if scoringMatrix[i][j] == maxScore:
                startPoint.append((i, j))

            j += 1
        i += 1

    alignment = getAlignment(scoringMatrix, pointerMatrix, startPoint, seq1, seq2)
    return (scoringMatrix, pointerMatrix, maxScore, alignment)




def find_kmers(k, str):
  kmers = [str[i : i+k] for i in range(len(str)-k+1)]
  print("K-mers in sequence:")
  print(kmers)
  return kmers


def helper_func(kmer, i, kmer_index):
  symbols = ['A', 'C', 'G', 'T']
  words = []
  words.append((kmer, kmer_index))
  kmer = list(kmer)

  char = kmer[i]

  for symbol in symbols:
    if char != symbol:
      kmer[i] = symbol

      words.append(("".join(kmer), kmer_index))

  print(words)
  return words



def find_synonyms(kmers, k):
  print("Finding synonyms to each k-mer in query sequence")
  neighbors = []
  for i in range(len(kmers)):
    print("\nCurrent Word: " + kmers[i])
    for j in range(k):
      words = helper_func(kmers[i], j, i)
      neighbors += words

  print
  print("Neighbors")
  print(neighbors)
  return neighbors


def calculate_keys(words, k):
  words_keys = []

  for word_tuple in words:
    word = word_tuple[0]
    key = 0

    # iterate char by char through word
    for i in range(len(word)):
      if word[i] == 'A':
        weight = 0
      elif word[i] == 'C':
        weight = 1
      elif word[i] == 'G':
        weight = 2
      elif word[i] == 'T':
        weight = 3

      key = key + weight * 4**(k-i-1)
    
    words_keys.append((word_tuple[0], word_tuple[1], key))

  return words_keys


def calculate_hsp_score(query_seq, db_seq, qidx, dbidx):
  if query_seq[qidx] == 'A' and db_seq[dbidx] == 'A':
    return 1
  elif (query_seq[qidx] == 'A' and db_seq[dbidx] == 'T') or (query_seq[qidx] == 'T' and db_seq[dbidx] == 'A'):
    return -1
  elif (query_seq[qidx] == 'A' and db_seq[dbidx] == 'G') or (query_seq[qidx] == 'G' and db_seq[dbidx] == 'A'):
    return -0.5
  elif (query_seq[qidx] == 'A' and db_seq[dbidx] == 'C') or (query_seq[qidx] == 'C' and db_seq[dbidx] == 'A'):
    return -1
  elif (query_seq[qidx] == 'T' and db_seq[dbidx] == 'T'):
    return 1
  elif (query_seq[qidx] == 'T' and db_seq[dbidx] == 'G') or (query_seq[qidx] == 'G' and db_seq[dbidx] == 'T'):
    return -1
  elif (query_seq[qidx] == 'T' and db_seq[dbidx] == 'C') or (query_seq[qidx] == 'C' and db_seq[dbidx] == 'T'):
    return -0.5
  elif (query_seq[qidx] == 'G' and db_seq[dbidx] == 'G'):
    return 1
  elif (query_seq[qidx] == 'G' and db_seq[dbidx] == 'C') or (query_seq[qidx] == 'C' and db_seq[dbidx] == 'G'):
    return -1
  elif (query_seq[qidx] == 'C' and db_seq[dbidx] == 'C'):
    return 1



def ungapped_extension(query_seq, db_seq, start_q, start_db, k):
  score = 0
  threshold = 1

  startidx_q = start_q
  startidx_db = start_db
  endidx_q = start_q + k
  endidx_db = start_db + k

  print("\n-----UNGAPPED EXTENSION between DB k-mer and Query k-mer------")
  print("Start query index: ", start_q)
  print("Start db index: ", start_db)


  # calculate score for length of seed
  for i in range(k):
    score = score + calculate_hsp_score(query_seq, db_seq, start_q+i, start_db+i)

  print("seed score: ", score)


  # extend to the right
  query_idx = start_q + k
  db_idx = start_db + k
  
  while(query_idx < len(query_seq) and db_idx < len(db_seq)):
    # if match or mismatch
    score = score + calculate_hsp_score(query_seq, db_seq, query_idx, db_idx)
    
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
  
  while query_idx >= 0 and db_idx >= 0:
    score = score + calculate_hsp_score(query_seq, db_seq, query_idx, db_idx)
    
    if score < threshold:
      # print("under threshold!")
      query_idx = query_idx + 1
      db_idx = db_idx + 1
      break
    
    query_idx = query_idx - 1
    db_idx = db_idx - 1

  if query_idx < 0 or db_idx < 0:
    query_idx += 1
    db_idx += 1

  startidx_db = db_idx
  startidx_q = query_idx

  print("Query sequence: ", query_seq[startidx_q:endidx_q+1],startidx_q, endidx_q)
  print("DB sequence: ", db_seq[startidx_db:endidx_db+1], startidx_db, endidx_db)
  return query_seq[startidx_q : endidx_q + 1], db_seq[startidx_db : endidx_db + 1], score



print("Query Sequence: " + query_seq)

# find all 3-mers of the query sequence
kmers = find_kmers(3, query_seq)

# find synonyms to all query sequence 3-mers
words = find_synonyms(kmers, 3)

# should only pick 3-mers that have a score above a certain threshold according to blosum matrix

words_keys = calculate_keys(words, 3)

query_bst = sorted(words_keys, key=lambda x: x[2])
print
print("Print each synonym and its calculated key")
print(query_bst)
print



#find all 3-mers in db sequence
kmers_db = find_kmers(3, db_seq)
for i in range(len(kmers_db)):
  kmers_db[i] = (kmers_db[i], i)
print("\nEach DB sequence k-mer and its position")
print(kmers_db)

kmers_db_keys = calculate_keys(kmers_db, 3)
print("\n\nEach DB sequence k-mer and its key")
print(kmers_db_keys)


print("\n\nFor each DB seq k-mer")
results = []
for kmer in kmers_db_keys:
  print
  # from db seq
  print("DB seq k-mer: ")
  print(kmer)
  output = [item for item in query_bst if item[2] == kmer[2]]

  print("\nCorresponding seeds in query seq:")
  # from query seq
  print(output)

  for op in output:
    str1, str2, score = ungapped_extension(query_seq, db_seq, op[1], kmer[1], 3)
    (scoringMat, pointerMat, maxScore, alignment) = alignUsingSW(str1, str2, 1, -1, 1)
    print("RESULT: ", alignment)
    results.append(alignment)

print("Resulting alignments:")
print(results)