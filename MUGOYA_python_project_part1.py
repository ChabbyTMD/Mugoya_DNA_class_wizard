import operator as op

standard_code = {
    "UUU": "F",
    "UUC": "F",
    "UUA": "L",
    "UUG": "L",
    "UCU": "S",
    "UCC": "S",
    "UCA": "S",
    "UCG": "S",
    "UAU": "Y",
    "UAC": "Y",
    "UAA": "*",
    "UAG": "*",
    "UGA": "*",
    "UGU": "C",
    "UGC": "C",
    "UGG": "W",
    "CUU": "L",
    "CUC": "L",
    "CUA": "L",
    "CUG": "L",
    "CCU": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAU": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGU": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "AUU": "I",
    "AUC": "I",
    "AUA": "I",
    "AUG": "M",
    "ACU": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAU": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGU": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GUU": "V",
    "GUC": "V",
    "GUA": "V",
    "GUG": "V",
    "GCU": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAU": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGU": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


aa_mol_weights = {
    "A": 89.09,
    "C": 121.15,
    "D": 133.1,
    "E": 147.13,
    "F": 165.19,
    "G": 75.07,
    "H": 155.16,
    "I": 131.17,
    "K": 146.19,
    "L": 131.17,
    "M": 149.21,
    "N": 132.12,
    "P": 115.13,
    "Q": 146.15,
    "R": 174.2,
    "S": 105.09,
    "T": 119.12,
    "V": 117.15,
    "W": 204.23,
    "X": 0,
    "Y": 181.19,
}


# Declaring seq parent class
class seq:
    # Call instance attributes
    def __init__(self, name, organism, sequence, type):
        self.name = name
        self.organism = organism
        self.sequence = sequence
        self.type = type

    # define the info function
    def info(self):
        print("{},{},{},{}".format(self.name, self.organism, self.type, self.sequence))

    # define the length function
    def length(self):
        print("{}".format(len(self.sequence)))

    # define the fasta_out function.

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "\n"
            + self.sequence
        )
        f.close()


# Protein Child Class.


class protein(seq):
    def __init__(self, name, organism, sequence, type, size):
        # Declare size attribute for protein child class
        self.size = size
        super().__init__(name, organism, sequence, type)

    def mol_weight(self):
        """
        Iterate through protein sequence and calculate molecular weight using reference weights from the aa_mol_weights dictionary and print the result.
        """
        total_mol_weight = 0
        for amino in self.sequence:
            total_mol_weight += aa_mol_weights[amino]
        print(total_mol_weight)

    def fasta_out(self):
        f = open("{}.fa".format(self.name), "w")
        f.write(
            ">"
            + self.name
            + "_"
            + self.organism
            + "_"
            + self.type
            + "_"
            + self.size
            + "\n"
            + self.sequence
        )
        f.close()


# Nucleotide Child Class
class nucleotide(seq):
    def __init__(self, name, organism, sequence, type):
        super().__init__(name, organism, sequence, type)

    def gc_content(self):
        """
        Using the countOf method in the operator python module, count all occurences of G and C in the nucleotide sequence, divided by the total length of the sequence and multiplied by 100 to give the percentage of GC content.
        """
        print(
            "{}%".format(
                (op.countOf(self.sequence, "G") + op.countOf(self.sequence, "C"))
                / len(self.sequence)
                * 100
            )
        )


# DNA Grand-Child Class
## NOTE: (reverse_complement method must be ran first before six_frames method)
class DNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        # Initialising the reverse complemented string to none
        self.rev_comp = None
        super().__init__(name, organism, sequence, type)

    def transcribe(self):
        print("{}".format(self.sequence.replace("T", "U")))
        return self.sequence.replace("T", "U")

    def reverse_complement(self):
        reverse_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
        # Turn nucleotide string into list of nucleotides
        dna_list = list(self.sequence)

        # reverse the list of nucleotides
        dna_list.reverse()

        rev_comp_list = []
        # Append to a list the complementary nucleotide to reverse complement list
        for nuc in dna_list:
            rev_comp_list.append(reverse_dict[nuc])

        # Concatenate all nucleotides in the reverse complement list into a string.
        rev_comp_seq = "".join(rev_comp_list)
        # Print and return the reverse complemented string.
        print(rev_comp_seq)
        self.rev_comp = rev_comp_seq

    def six_frames(self):
        frameF = self.sequence
        frameR = self.rev_comp
        print("Forward strand reading frames for sequence {} \n".format(self.sequence))

        # Using string slicing to get the 3 reading frames for the forward and reverse strands
        # Forward reading frames.
        for slice_indexF in range(0, 3, 1):
            print(frameF[slice_indexF:])

        # Reverse reading frames
        print("\n")
        print("Reverse strand reading frames for sequence {} \n".format(self.rev_comp))
        for slice_indexR in range(0, 3, 1):
            print(frameR[slice_indexR:])


# RNA Grand-Child Class
## NOTE:(start method must be called first before translate method)
class RNA(nucleotide):
    def __init__(self, name, organism, sequence, type):
        # Initialize start index with None
        self.start_index = None
        super().__init__(name, organism, sequence, type)

    def start(self):
        print("{}".format(self.sequence.find("AUG")))
        self.start_index = self.sequence.find("AUG")

    def translate(self):
        """
        Loop through the RNA transcript string in 3 character chunks starting from the first occurence of start codon AUG and create a list of codons upto a stop codon defined in the standard_code dictionary.
        Subsequently, the next loop translates the codons to individual amino acids that are concatenated into one protein sequence.
        """
        transcript = self.sequence
        codon_list = []
        protein_sequence = ""
        for i in range(self.start_index, len(transcript), 3):
            codon = transcript[i : i + 3]
            if len(codon) == 3:
                codon_list.append(codon)

        for codon in codon_list:
            if standard_code[codon] == "*":
                break
            else:
                protein_sequence += standard_code[codon]

        print(protein_sequence)
        return protein_sequence


# Project Part1

# Part B:

# DNA Section: Providing sample name, DNA sequence, organism and DNA sequence type to create fasta file containing these data, return a reverse complement of the sequence as well as six reading frame.
uidA_DNA = DNA(
    name="uidA_DNA",
    sequence="CGCATGTTACGTCCTGTAGAAACCCCAACCCGTGAAATCAAAAAA",
    organism="Bacteria",
    type="DNA",
)

uidA_DNA.fasta_out()
print("\n")
print("uidA Reverse complement sequence:")
uidA_DNA.reverse_complement()
print("\n")
print("uidA Reading frames:")
print("\n")
uidA_DNA.six_frames()
print("\n")

# RNA Section: Utilise the transcribe method from DNA grand-child class to create an RNA sequence. Create an RNA object that takes name, RNA sequence, organism and RNA nucleotide type and writes a fasta file containing these data as well as translating the RNA and returning a protein sequence.
print("RNA Transcript:")
uidA_RNA = uidA_DNA.transcribe()

uiDA_RNA = RNA(name="uidA_RNA", sequence=uidA_RNA, organism="Bacteria", type="RNA")

uiDA_RNA.fasta_out()
print("Index of the start codon:")
uiDA_RNA.start()
print("\n")
print("Protein Sequence:")
uidA_protein = uiDA_RNA.translate()

# Protein Section: Accept protein sequence as input generated from RNA grand child class and create a protein child class object with protein name, sequence, organism, sequence type and protein sequence length. Subseuently create a fasta file containing these data in addition to calculating the molecular weight of the protein.
uiDA_protein = protein(
    name="uidA_protein",
    sequence=uidA_protein,
    organism="Bacteria",
    type="protein",
    size=str(len(uidA_protein)),
)

uiDA_protein.fasta_out()
print("Molecular weight:")
uiDA_protein.mol_weight()
