def read_dna_sequence_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            dna_sequence = file.read().strip()
            return dna_sequence
    except FileNotFoundError:
        print(f"Error: File not found at {file_path}")
        return None

def transcribe_to_rna(dna_sequence):
    # Replace 'T' with 'U' to transcribe to RNA
    rna_sequence = dna_sequence.replace('T', 'U')
    return rna_sequence

def translate_to_protein(rna_sequence):
    # Define a dictionary mapping RNA codons to amino acids
    codon_table = {
        'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
        'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
        'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S',
        'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
        'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W',
        'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }

    protein_sequence = ""
    
    # Iterate over the RNA sequence, reading three bases at a time
    for i in range(0, len(rna_sequence), 3):
        codon = rna_sequence[i:i+3]
        
        # Stop translation if a stop codon is encountered
        if codon in ('UAA', 'UAG', 'UGA'):
            break

        # Lookup the amino acid for the current codon
        amino_acid = codon_table.get(codon, 'X')  # 'X' for unknown amino acid
        
        # Append the amino acid to the protein sequence
        protein_sequence += amino_acid

    return protein_sequence

# Example usage:
if __name__ == "__main__":
    # Input the file path containing the DNA sequence
    file_path = input("Enter the file path containing the DNA sequence: ")
    
    # Read the DNA sequence from the file
    dna_sequence = read_dna_sequence_from_file(file_path)
    
    if dna_sequence:
        # Transcribe the DNA sequence to RNA
        rna_sequence = transcribe_to_rna(dna_sequence)
        print(f"RNA Sequence: {rna_sequence}")
        
        # Translate the RNA sequence to protein
        protein_sequence = translate_to_protein(rna_sequence)
        print(f"Protein Sequence: {protein_sequence}")
