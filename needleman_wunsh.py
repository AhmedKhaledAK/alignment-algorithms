def alignUsingNW(seq1, seq2, match=1, mismatch=-1, gapPenalty=1):

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