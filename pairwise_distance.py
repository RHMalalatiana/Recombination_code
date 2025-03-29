import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations

def compute_pairwise_identity(sequences, step_size):
    """
    Compute pairwise identity between sequences in a given interval.

    Parameters:
    - sequences (list of str): Aligned sequences.
    - step_size (int): Length of subsequence to consider.

    Returns:
    - pairwise_identity (dict): Dictionary with pairs as keys and list of pairwise identity.
    """
    num_seqs = len(sequences)
    alignment_length = len(sequences[0])
    pairs = list(combinations(range(num_seqs), 2))  # Generate all unique pairs of sequence indices
    pairwise_identity = {pair: [] for pair in pairs}
    steps = []

    # Slide the window across the alignment
    for start in range(0, alignment_length-step_size+1, step_size):
        steps.append(start+1)
        for (i, j) in pairs:
            segment1 = sequences[i][start:start + step_size]
            segment2 = sequences[j][start:start + step_size]
            
            # Calculate pairwise identity as proportion of matches
            ident = sum(1 for a, b in zip(segment1, segment2) if a == b)
            pair_id = ident/step_size
            pairwise_identity[(i, j)].append(pair_id)
    
    return steps, pairwise_identity

def plot_pairwise_identity(steps, pairwise_identity, sequence_labels):
    """
    Plot pairwise identity across alignment for each pair of sequences.

    Parameters:
    - steps (list of int): Start positions of each window.
    - pairwise_identity (dict): Dictionary with pairs as keys and list of pairwise identity as values.
    - sequence_labels (list of str): Labels for sequences.
    """
    plt.figure(figsize=(12, 6))
    
    # Plot each pair's identity
    for (i, j), identity in pairwise_identity.items():
        label = f"{sequence_labels[i]}-{sequence_labels[j]}"
        plt.plot(steps, identity, marker='*',label=label)
    
    # Label the plot
    plt.title("Pairwise Identity Across Alignment for Each Sequence Pair")
    plt.xlabel("Position along alignment")
    plt.ylabel("Pairwise identity")
    plt.legend(title="Sequence Pairs")
    plt.grid(True)
    plt.ylim(0, 1.5)
    plt.show()
import random

def generate_dna_sequences(num_sequences, length):
    """
    Generate a specified number of random DNA sequences of a given length.

    Parameters:
    - num_sequences (int): Number of DNA sequences to generate.
    - length (int): Length of each DNA sequence.

    Returns:
    - sequences (list of str): List of generated DNA sequences.
    """
    nucleotides = ['A', 'T', 'G', 'C']
    sequences = []
    
    for _ in range(num_sequences):
        sequence = ''.join(random.choice(nucleotides) for _ in range(length))
        sequences.append(sequence)
    
    return sequences

""" # Example usage
sequences = generate_dna_sequences(3,2000) # Example sequences
sequence_labels = ["Seq1", "Seq2", "Seq3"]         # Sequence names                   
step_size = 20                                    # Set step size

# Compute pairwise identity
steps, pairwise_identity = compute_pairwise_identity(sequences, step_size)

# Plot the pairwise identity
plot_pairwise_identity(steps, pairwise_identity, sequence_labels)
 """