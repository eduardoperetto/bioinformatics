import os
import datetime

class SequenceReader:
    def __init__(self, file_path):
        self.file_path = file_path

    def read_sequences(self):
        sequences = {}
        with open(self.file_path, 'r') as file:
            identifier = ''
            sequence = ''
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if identifier:
                        sequences[identifier] = sequence
                        sequence = ''
                    identifier = line[1:]
                else:
                    sequence += line
            if identifier:
                sequences[identifier] = sequence
        return sequences


class AlignmentResult:
    def __init__(self, alignment_seq1, alignment_seq2, score):
        self.alignment_seq1 = alignment_seq1
        self.alignment_seq2 = alignment_seq2
        self.score = score

    def calculate_identity(self):
        matches = sum(1 for a, b in zip(self.alignment_seq1, self.alignment_seq2) if a == b)
        return (matches / len(self.alignment_seq1)) * 100


class NeedlemanWunschAligner:
    def __init__(self, match_score=2, mismatch_score=-2, gap_penalty=-4):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_penalty = gap_penalty

    def align(self, seq1, seq2):
        m, n = len(seq1), len(seq2)
        score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

        for i in range(m + 1):
            score_matrix[i][0] = self.gap_penalty * i
        for j in range(n + 1):
            score_matrix[0][j] = self.gap_penalty * j

        for i in range(1, m + 1):
            for j in range(1, n + 1):
                match = score_matrix[i - 1][j - 1] + (self.match_score if seq1[i - 1] == seq2[j - 1] else self.mismatch_score)
                delete = score_matrix[i - 1][j] + self.gap_penalty
                insert = score_matrix[i][j - 1] + self.gap_penalty
                score_matrix[i][j] = max(match, delete, insert)

        alignment_seq1, alignment_seq2 = '', ''
        i, j = m, n
        while i > 0 and j > 0:
            current_score = score_matrix[i][j]
            diagonal_score = score_matrix[i - 1][j - 1]
            up_score = score_matrix[i - 1][j]
            left_score = score_matrix[i][j - 1]

            if current_score == diagonal_score + (self.match_score if seq1[i - 1] == seq2[j - 1] else self.mismatch_score):
                alignment_seq1 = seq1[i - 1] + alignment_seq1
                alignment_seq2 = seq2[j - 1] + alignment_seq2
                i -= 1
                j -= 1
            elif current_score == up_score + self.gap_penalty:
                alignment_seq1 = seq1[i - 1] + alignment_seq1
                alignment_seq2 = '-' + alignment_seq2
                i -= 1
            else:
                alignment_seq1 = '-' + alignment_seq1
                alignment_seq2 = seq2[j - 1] + alignment_seq2
                j -= 1

        while i > 0:
            alignment_seq1 = seq1[i - 1] + alignment_seq1
            alignment_seq2 = '-' + alignment_seq2
            i -= 1
        while j > 0:
            alignment_seq1 = '-' + alignment_seq1
            alignment_seq2 = seq2[j - 1] + alignment_seq2
            j -= 1

        final_score = score_matrix[m][n]
        return AlignmentResult(alignment_seq1, alignment_seq2, final_score)


def save_alignment_output(result, identity):
    if not os.path.exists('output'):
        os.makedirs('output')
    filename = datetime.datetime.now().strftime("%H_%M_%d_%m.txt")
    output_path = os.path.join('output', filename)
    
    with open(output_path, 'w') as output_file:
        output_file.write("\nAlignment Results:\n")
        output_file.write(f"Score: {result.score}\n")
        output_file.write(f"Identity: {identity:.2f}%\n")
        output_file.write(f"Sequence 1 Alignment:\n{result.alignment_seq1}\n")
        output_file.write(f"Sequence 2 Alignment:\n{result.alignment_seq2}\n")


def main():
    file_path = 'input.fasta'
    reader = SequenceReader(file_path)
    sequences = reader.read_sequences()

    ids = list(sequences.keys())
    if len(ids) < 2:
        save_alignment_output("Not enough sequences to align.", None)
        return

    seq1, seq2 = sequences[ids[0]], sequences[ids[1]]
    aligner = NeedlemanWunschAligner()
    result = aligner.align(seq1, seq2)
    identity = result.calculate_identity()

    save_alignment_output(result, identity)


if __name__ == "__main__":
    main()
