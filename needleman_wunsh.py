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



def alignUsingNW(seq1, seq2, match=1, mismatch=-1, gapPenalty=1):

    subMatrix = getSubMatrix(match, mismatch, seq1, seq2)

    H, W = len(seq1) + 1, len(seq2) + 1
    scoringMatrix = [[0 for i in range(W)] for y in range(H)]
    pointerMatrix = [[(-1,-1,-1, -1) for i in range(W)] for y in range(H)]

    i = 0
    while i < H:
        j = 0
        while j < W:
            if i == 0 and j == 0:
                continue
            if i == 0:
                scoringMatrix[i][j] = scoringMatrix[i][j - 1] - gapPenalty
            elif j == 0:
                scoringMatrix[i][j] = scoringMatrix[i - 1][j] - gapPenalty
            else:
                pass