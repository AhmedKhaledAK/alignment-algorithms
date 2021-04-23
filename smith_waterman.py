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

            if maxScore <= scoringMatrix[i][j]:
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

def main():
    (scoringMat, pointerMat, maxScore, alignment) = alignUsingSW("ACAGC", "ACGGA", 1, -1, 1)

    print(scoringMat)
    print("------")
    print(pointerMat)
    print("------")
    print(alignment)
    print("------")
    print(maxScore)


main()