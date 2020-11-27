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
    str1 = str2 = ""
    i, j = startPoint[0], startPoint[1]
    while i != 0 and j != 0:
        if pointerMatrix[i][j][2] == 0:
            str1 += seq1[i - 1]
            str2 += seq2[j - 1]
            print(scoringMatrix[i][j])
        elif pointerMatrix[i][j][2] == 1:
            str1 += seq1[i - 1]
            str2 += "-"
            print(scoringMatrix[i][j])
        elif pointerMatrix[i][j][2] == 2:
            str1 += "-"
            str2 += seq2[j - 1]
            print(scoringMatrix[i][j])
        print(pointerMatrix[i][j])
        hold = i
        i = pointerMatrix[i][j][0]
        j = pointerMatrix[hold][j][1]
    return (str1[::-1], str2[::-1])

def alignUsingNW(seq1, seq2, match=1, mismatch=-1, gapPenalty=-1):

    subMatrix = getSubMatrix(match, mismatch, seq1, seq2)

    H, W = len(seq1) + 1, len(seq2) + 1
    scoringMatrix = [[0 for i in range(W)] for y in range(H)]
    pointerMatrix = [[(-1,-1,-1, -1) for i in range(W)] for y in range(H)]

    i = 0
    while i < H:
        j = 0
        while j < W:
            if i == 0 and j == 0:
                j += 1
                continue
            if i == 0:
                scoringMatrix[i][j] = scoringMatrix[i][j - 1] + gapPenalty
            elif j == 0:
                scoringMatrix[i][j] = scoringMatrix[i - 1][j] + gapPenalty
            else:
                scoringMatrix[i][j] = max(scoringMatrix[i - 1][j - 1] + subMatrix[i - 1][j - 1], scoringMatrix[i - 1][j] + gapPenalty, scoringMatrix[i][j - 1] + gapPenalty)

            if scoringMatrix[i][j] == scoringMatrix[i - 1][j - 1] + subMatrix[i - 1][j - 1]:
                pointerMatrix[i][j] = (i - 1, j - 1, 0, scoringMatrix[i][j])
            elif scoringMatrix[i][j] == scoringMatrix[i - 1][j] + gapPenalty:
                pointerMatrix[i][j] = (i - 1, j, 1, scoringMatrix[i][j])
            elif scoringMatrix[i][j] == scoringMatrix[i][j - 1] + gapPenalty:
                pointerMatrix[i][j] = (i, j - 1, 2, scoringMatrix[i][j])

            
            j += 1
        i += 1


    alignment = getAlignment(scoringMatrix, pointerMatrix, (H - 1, W - 1), seq1, seq2)
    return (scoringMatrix, pointerMatrix, alignment)

def main():
    (scoringMatrix, pointerMatrix, alignment) = alignUsingNW("GATTACA", "GCATGCU")
    
    print(scoringMatrix)
    print("-------")
    print(pointerMatrix)
    print("-------")
    print(alignment)

main()