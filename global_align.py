import sys         # for reading args from command line 

# fills in dp matrix and finds score
def dp_align(dp_matrix, sub_matrix, gap_penalty, idx1, idx2): 
    for i in range (1,idx1+1):
        for j in range (1,idx2+1):
            match = dp_matrix[i-1][j-1] + sub_matrix[seq1[i-1]][seq2[j-1]]
            gap1 = dp_matrix[i-1][j] + gap_penalty
            gap2 = dp_matrix[i][j-1] + gap_penalty
            dp_matrix[i][j] = max(match,gap1,gap2)
    return dp_matrix[idx1][idx2]

def traceback(dp_matrix, sub_matrix, gap_penalty, idx1, idx2):
    alignment1 = ""
    alignment2 = ""
    while idx1 > 0 or idx2 > 0:
        val = dp_matrix[idx1][idx2]
        if idx1-1 >= 0 and idx2-1 >= 0 and (not dp_matrix[idx1-1][idx2-1] == None) and val == dp_matrix[idx1-1][idx2-1] + sub_matrix[seq1[idx1-1]][seq2[idx2-1]]:
            alignment1 = str(seq1[idx1-1]) + alignment1
            alignment2 = str(seq2[idx2-1]) + alignment2
            idx2 -= 1
            idx1 -= 1
        elif idx1-1 >= 0 and (not dp_matrix[idx1-1][idx2] == None) and val == dp_matrix[idx1-1][idx2] + gap_penalty:
            alignment1 = str(seq1[idx1-1]) + alignment1
            alignment2 = "-" + alignment2
            idx1 -= 1 
        elif idx2-1 >= 0 and (not dp_matrix[idx1][idx2-1] == None) and val == dp_matrix[idx1][idx2-1] + gap_penalty:
            alignment2 = str(seq2[idx2-1]) + alignment2
            alignment1 = "-" + alignment1
            idx2 -= 1
    return alignment1, alignment2

if __name__ == "__main__":
    # retrieve needed filenames and gap-penalty from the commandline args
    seq_file1 = sys.argv[1]
    seq_file2 = sys.argv[2]
    sub_matrix_file = sys.argv[3]
    gap_penalty = int(sys.argv[4])

    # read in sequences to align from fasta files 
    with open(seq_file1) as f:
        lines = f.readlines()
        seq1 = lines[1]

    with open(seq_file2) as f:
        lines = f.readlines()
        seq2 = lines[1]
    
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

    # create initial dynamic programming matrix 
    dp_matrix = [[None for y in range(len(seq2)+1)] for x in range(len(seq1)+1)]
    dp_matrix[0][0] = 0
    for i in range(1, len(seq2)+1):
        dp_matrix[0][i] = dp_matrix[0][i-1] + gap_penalty
    for i in range(1, len(seq1)+1):
        dp_matrix[i][0] = dp_matrix[i-1][0] + gap_penalty

    # creates dp matrix and finds alignment score
    print "The optimal alignment between given sequences has score", dp_align(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))

    # finds and prints alignment
    a1, a2 = traceback(dp_matrix, sub_matrix, gap_penalty, len(seq1), len(seq2))
    print a1
    print a2
  
