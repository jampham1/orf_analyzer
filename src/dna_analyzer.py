
CODON_TABLE = {
    "ATA":"I","ATC":"I","ATT":"I","ATG":"M",
    "ACA":"T","ACC":"T","ACG":"T","ACT":"T",
    "AAC":"N","AAT":"N","AAA":"K","AAG":"K",
    "AGC":"S","AGT":"S","AGA":"R","AGG":"R",
    "CTA":"L","CTC":"L","CTG":"L","CTT":"L",
    "CCA":"P","CCC":"P","CCG":"P","CCT":"P",
    "CAC":"H","CAT":"H","CAA":"Q","CAG":"Q",
    "CGA":"R","CGC":"R","CGG":"R","CGT":"R",
    "GTA":"V","GTC":"V","GTG":"V","GTT":"V",
    "GCA":"A","GCC":"A","GCG":"A","GCT":"A",
    "GAC":"D","GAT":"D","GAA":"E","GAG":"E",
    "GGA":"G","GGC":"G","GGG":"G","GGT":"G",
    "TCA":"S","TCC":"S","TCG":"S","TCT":"S",
    "TTC":"F","TTT":"F","TTA":"L","TTG":"L",
    "TAC":"Y","TAT":"Y","TAA":"*","TAG":"*",
    "TGC":"C","TGT":"C","TGA":"*","TGG":"W",
}

# Retrieves fasta file and assigns value to "sequence" variable
def read_fasta(filename):
    sequences = {}
    header = None

    with open(filename) as file:
        for line in file:
            line = line.strip()

            if line.startswith(">"):
                header = line[1:]
                sequences[header] = ""
            else:
                sequences[header] += line.upper()

    return sequences

# Generates the complementary strand to original sequence
def reverse_complement(sequence):
    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complement[base] for base in reversed(sequence))

#Analyzes ORFs for start/end position, reading frame, AA count w/ names
def find_orfs(sequence):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}

    orfs = []

    for frame in range(3):

        for i in range(frame, len(sequence) - 2):

            codon = sequence[i:i+3]

            if codon == start_codon:

                for j in range(i, len(sequence) - 2, 3):

                    codon = sequence[j:j+3]

                    if codon in stop_codons:

                        start = i
                        end = j + 3

                        orf_seq = sequence[start:end]
                        protein = translate_dna(orf_seq)
                        aa_length = len(protein)

                        orfs.append((start, end, frame, aa_length, protein))

    return orfs

# Analyzes both original and complementary sequences; adds to ORF list
def analyze_sequence(name, sequence, output_file):

    strands = {
        "forward": sequence,
        "reverse": reverse_complement(sequence)
    }

    for strand_name, seq in strands.items():

        orfs = find_orfs(seq)

        for start, end, frame, aa, protein in orfs:

            output_file.write(
                f"{name}\t{strand_name}\t{frame}\t{start}\t{end}\t{aa}\t{protein}\n"
            )

# Analyzes input sequence, generates output in a new file called "orf_results.tsv"
def analyze_fasta(input_file, output_file="orf_results.tsv"):

    sequences = read_fasta(input_file)

    with open(output_file, "w") as out:

        out.write("Sequence\tStrand\tFrame\tStart\tEnd\tAminoAcids\tProtein\n")

        for name, seq in sequences.items():
            analyze_sequence(name, seq, out)

    print(f"Results saved to {output_file}")

# Converts codons directs to AA name
def translate_dna(dna_sequence):
    protein = ""

    for i in range(0, len(dna_sequence) - 2, 3):
        codon = dna_sequence[i:i+3]
        amino_acid = CODON_TABLE.get(codon, "X")  # X = unknown

        if amino_acid == "*":  # stop codon
            break

        protein += amino_acid

    return protein


# run the script
analyze_fasta("../data/dna.fasta")

# AFTER THIS POINT, USES OUTPUT FILE DATA FROM "orf_results.tsv"
def load_orf_results(filename):
    orf_lengths = []
    orf_frames = []

    with open(filename) as f:
        next(f)  # skip header

        for line in f:
            parts = line.strip().split("\t")

            aa_length = int(parts[5])
            frame = int(parts[2])

            orf_lengths.append(aa_length)
            orf_frames.append(frame)

    return orf_lengths, orf_frames


import matplotlib.pyplot as plt

def plot_orf_lengths(plot_lengths):

    plt.figure()
    plt.hist(plot_lengths, bins=20)

    plt.xlabel("Amino Acid Length")
    plt.ylabel("Frequency")
    plt.title("ORF Length Distribution")

    plt.savefig("orf_length_distribution.png")  # save for GitHub
    plt.show()


def plot_frame_distribution(plot_frames):

    frame_counts = {0: 0, 1: 0, 2: 0}

    for f in plot_frames:
        frame_counts[f] += 1

    plt.figure()
    plt.bar(frame_counts.keys(), frame_counts.values())

    plt.xlabel("Reading Frame")
    plt.ylabel("Number of ORFs")
    plt.title("ORF Distribution by Frame")

    plt.savefig("frame_distribution.png")
    plt.show()


lengths, frames = load_orf_results("../results/orf_results.tsv")

plot_orf_lengths(lengths)
plot_frame_distribution(frames)

plt.yscale("log")