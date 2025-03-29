# SIMULATED DATA GENERATION
import copy
import math as M
import random
from pairwise_distance import *


def jukesCantor(genome,gen_time):
  """
  Input: genome, a character string consisting of a, c, g, and t
         gen_time, genetic time
  Output: a randomly modified string according to the Juke's cantor model
    The following coupled ODEs are integrated from t=0 to t=gen_time

    dot a = -a + (c + g + t) /3
    dot c = -a + (g + t + a) /3
    dot g = -g + (t + a + c) /3
    dot t = -t + (a + c + g) /3
  
    using the solution 
 
    a(t)=(3/4)*exp(-t)*(a(0)-(c(0)+g(0)+t(0))/3)
        +(1/4)*(a(0)+c(0)+g(0)+t(0))

    and similarly for c(t), g(t), and t(t).

    In other words, if 'a' is the initial base, its persistence
    probability after a time t is 

    3/4 exp(-t) + 1/4

    and likewise for the other three base inputs.

  """
  if type(gen_time) !=  float :
    raise Exception("Time must be a float")
  for b in genome:
    if not ( b=='a' or b=='c' or b=='g' or b=='t' ):
       raise Exception("Invalid input genome")
  persistence_prob=(3.*M.exp(-gen_time)+1.)/4.
  new_genome=[]
  for base in genome:
    p=random.random()
    if p < persistence_prob:
      new_base=base
    else: 
      if base=='a':
        seq=('c','g','t')
      if base=='c':
        seq=('a','g','t')
      if base=='g':
        seq=('a','c','t')
      if base=='t':
        seq=('a','c','g')
      new_base=random.choice(seq)
    new_genome.append(new_base)
  return(new_genome)

# The following routine, which applies mutations randomly, chooses the Jukes-Cantor algorithm,
# but some other model can be used by changing the function below.

def evolve(genome,time):
   return(jukesCantor(genome,time))

def generateHelper(initialGenome,listIn):
 time=listIn[0]
 newGenome=evolve(initialGenome,time)
 listIn[0]=newGenome
 if not ( len(listIn) == 1 or len(listIn) == 3):
    raise Exception("Lists must have one or three elements") 
 if len(listIn) == 1:
   return
 else:
   generateHelper(newGenome,listIn[1])
   generateHelper(newGenome,listIn[2])
   return 

def generateDriver(initialGenome,listIn):
  """
  This algorithm assumes a tree with input of the form:

  [ genTime, leftBranch, rightBranch]

  where the branches can either be tripleton lists of the format above
  or leaves represtented as a single nonnegative number enclosed in a list.
  genTime above is a non-negative real number. A deep copy of the input list is made,
  and the genetic times are replaced with the stochastically evolved genomes.

  An initial genone at the root of the binary tree is required as input, and this
  must be in a format of a list of base charcters, i.e., 'a', 'c', 'g', or 't'.

  """
  listInBis=copy.deepcopy(listIn)
  generateHelper(initialGenome,listInBis)
  return(listInBis)
             
def extract_genomes(tree):
 """
 This function extract the genome on the output of generateDriver
 """
 if len(tree)==1 :
   return(tree)
 else :
   return(extract_genomes(tree[1])+extract_genomes(tree[2]))
 
def rec(parentA,parentB,b):
  """
  This function create a recombinant genome from parentA and parentB
  """
  return(parentA[:b]+parentB[b:])


def simulate_recombination(seq1, seq2, recombination_rate):
    """
    Simulates recombination along two sequences based on a given recombination rate.
    
    Args:
        seq1 (str): First DNA sequence (parent 1).
        seq2 (str): Second DNA sequence (parent 2).
        recombination_rate (float): Probability of recombination at each position (between 0 and 1).
    
    Returns:
        str: Recombinant sequence.
        list: List of recombination points.
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length.")

    length = len(seq1)
    recombinant_seq = []
    current_parent = seq1  # Start from first parent
    recombination_points = []

    for i in range(length):
        # Decide whether to switch parent sequences based on recombination rate
        if random.random() < recombination_rate:
            current_parent = seq2 if current_parent == seq1 else seq1
            recombination_points.append(i)  # Store recombination position

        recombinant_seq.append(current_parent[i])  # Append current base

    return "".join(recombinant_seq)

def hamming_distance(seq1,seq2):
    """
    This function compute the difference between two sequences
    """
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    sum=0
    for i in range(len(seq1)):
        if seq1[i]!=seq2[i]:
            sum+=1
    return sum