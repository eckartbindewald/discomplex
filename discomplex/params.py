

# Size/Volume: Small, Medium, Large (1, 2, 3)
# Polarity: Non-polar, Polar (0, 1)
# Charge: Acidic, Neutral, Basic (-1, 0, 1)
# Hydrophobicity: Hydrophobic, Neutral, Hydrophilic (-1, 0, 1)
# Aromaticity: Non-aromatic, Aromatic (0, 1)
aa_encoding= {
    'A': [1, 0, 0, -1, 0],  # Alanine
    'R': [3, 1, 1, 0, 0],   # Arginine
    'N': [2, 1, 0, 0, 0],   # Asparagine
    'D': [2, 1, -1, 0, 0],  # Aspartic acid
    'C': [2, 1, 0, -2, 0],  # Cysteine, more hydrophobic
    'Q': [2, 1, 0, 0, 0],   # Glutamine
    'E': [2, 1, -1, 0, 0],  # Glutamic acid
    'G': [1, 0, 0, 0, 0],   # Glycine
    'H': [3, 1, 0, 0, 1],   # Histidine
    'I': [3, 0, 0, -2, 0],  # Isoleucine, more hydrophobic
    'L': [3, 0, 0, -1, 0],  # Leucine
    'K': [3, 1, 1, 1, 0],   # Lysine
    'M': [3, 0, 0, -2, 0],  # Methionine, more hydrophobic
    'F': [3, 0, 0, -2, 1],  # Phenylalanine, more hydrophobic
    'P': [2, 0, 0, 0, 0],   # Proline
    'S': [2, 1, 0, 0, 0],   # Serine
    'T': [2, 1, 0, -1, 0],  # Threonine, slightly hydrophobic
    'W': [3, 1, 0, -2, 1],  # Tryptophan, more hydrophobic
    'Y': [3, 1, 0, -1, 1],  # Tyrosine, less hydrophobic
    'V': [3, 0, 0, -1, 0]   # Valine
}