import random

def generate_dna(length=2000):
    bases = ["A", "T", "G", "C"]
    return "".join(random.choice(bases) for _ in range(length))


def write_fasta(sequence, filename="output.fasta", header=">sequence_1"):

    with open(filename, "w") as f:
        f.write(header + "\n")

        # wrap sequence every 60 characters (standard FASTA format)
        for i in range(0, len(sequence), 60):
            f.write(sequence[i:i+60] + "\n")


# generate sequence
dna_sequence = generate_dna(2000)

# write to fasta file
write_fasta(dna_sequence, "dna.fasta")

print("FASTA file created: dna.fasta")