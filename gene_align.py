# fills in dp matrix and finds score
def dp_align(dp_matrix, sub_matrix, gap_penalty, idx1, idx2): 
    for i in range (1,idx1+1):
        for j in range (1,idx2+1):
            match = dp_matrix[i-1][j-1] + sub_matrix[seq1[i-1]][seq2[j-1]]
            gap1 = dp_matrix[i-1][j] + gap_penalty
            gap2 = dp_matrix[i][j-1] + gap_penalty
            dp_matrix[i][j] = max(match,gap1,gap2)
    return dp_matrix[idx1][idx2]

# uses dp matrix to find alignment
def traceback(dp_matrix, sub_matrix, gap_penalty, idx1, idx2):
    score = 0
    alignment1 = ""
    alignment2 = ""
    while idx1 > 0 or idx2 > 0:
        val = dp_matrix[idx1][idx2]
        if idx1-1 >= 0 and idx2-1 >= 0 and (not dp_matrix[idx1-1][idx2-1] == None) and val == dp_matrix[idx1-1][idx2-1] + sub_matrix[seq1[idx1-1]][seq2[idx2-1]]:
            if seq1[idx1-1] == seq2[idx2-1]:
                score+=1
            else:
                score -=1
            alignment1 = str(seq1[idx1-1]) + alignment1
            alignment2 = str(seq2[idx2-1]) + alignment2
            idx2 -= 1
            idx1 -= 1
        elif idx1-1 >= 0 and (not dp_matrix[idx1-1][idx2] == None) and val == dp_matrix[idx1-1][idx2] + gap_penalty:
            score -=1
            alignment1 = str(seq1[idx1-1]) + alignment1
            alignment2 = "-" + alignment2
            idx1 -= 1 
        elif idx2-1 >= 0 and (not dp_matrix[idx1][idx2-1] == None) and val == dp_matrix[idx1][idx2-1] + gap_penalty:
            score -=1
            alignment2 = str(seq2[idx2-1]) + alignment2
            alignment1 = "-" + alignment1
            idx2 -= 1
    return score, alignment1, alignment2

# read in sequences to align from fasta files
def read_sequences(seq_file1, seq_file2):
    with open(seq_file1) as f1:
        with open(seq_file2) as f2:
            f1.readline()
            f2.readline()
            seq1 = ""
            seq2 = ""
            for line in f1:
                 seq1 += line.strip()
            for line in f2:
                seq2 += line.strip()
            #i = 0
            #length = 10000
            #while i < length:
            #    c1 = f1.read(1)
            #    c2 = f2.read(1)
             #   if c1 == '' or c2 == '':
             #       i = length
             #   if c1.upper() == "A" or c1.upper() == "C" or c1.upper() == "T" or c1.upper() == "G":
              #      if c2.upper() == "A" or c2.upper() == "C" or c2.upper() == "T" or c2.upper() == "G":
              #          i+=1
               #         seq1 += c1.upper()
               #         seq2 += c2.upper()
    return seq1, seq2
 
# create initial dynamic programming matrix 
def matrix_init(seq1, seq2, gap_penalty):
    dp_matrix = [[None for y in range(len(seq2)+1)] for x in range(len(seq1)+1)]
    dp_matrix[0][0] = 0
    for i in range(1, len(seq2)+1):
        dp_matrix[0][i] = dp_matrix[0][i-1] + gap_penalty
    for i in range(1, len(seq1)+1):
        dp_matrix[i][0] = dp_matrix[i-1][0] + gap_penalty
    return dp_matrix

if __name__ == "__main__":
    sub_matrix_file = "subs.txt"
    gap_penalty = -1
    length = 10

    # read in substituion matrix from file
    with open(sub_matrix_file) as f:
        sub_matrix = {} 
        letters = f.readline().split()
        for letter in letters:
            sub_matrix[letter] = {}
        for line in f:
            vals = line.split()
            for idx, val in enumerate(vals[1:]):
                sub_matrix[letters[idx]][vals[0]] = int(val)

    seq_file1 = "Data/ND4L_Canis_Lupus_Familiaris.fna"
    seq_file2 = "Data/ND4L_Panthera_Tigris.fna"
    seq1, seq2 = read_sequences(seq_file1, seq_file2)

    dp_matrix = matrix_init(seq1, seq2, gap_penalty)
    print "\nCANIS LUPUS FAMILIARIS VS. PANTHERA TIGRIS"
    dp_align(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))
    score1, a11, a12 = traceback(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))
    print "The optimal alignment between given sequences has score", score1

    seq_file1 = "Data/ND4L_Felis_Catus.fna"
    seq_file2 = "Data/ND4L_Panthera_Tigris.fna"
    seq1, seq2 = read_sequences(seq_file1, seq_file2)

    dp_matrix = matrix_init(seq1, seq2, gap_penalty)
    print "\nFELIS CATUS VS. PANTHERA TIGRIS"
    dp_align(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))
    score2, a21, a22 = traceback(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))
    print "The optimal alignment between given sequences has score", score2
 
    #write alignment results to a file
    with open('ND4L_results.txt', 'w') as f:
        f.write("CANIS LUPUS FAMILIARIS VS. PANTHERA TIGRIS")
        f.write("\nThe optimal alignment between given sequences has score " + str(score1))
        f.write("\n" + a11)
        f.write("\n" + a12)
        f.write("\n\nFELIS CATUS VS. PANTHERA TIGRIS")
        f.write("\nThe optimal alignment between given sequences has score " + str(score2))
        f.write("\n" + a21)
        f.write("\n" + a22)
  
  
