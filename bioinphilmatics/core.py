import collections
import re

CODON_TABLE = { "UUU": "F",      "CUU": "L",      "AUU": "I",      "GUU": "V",
                "UUC": "F",      "CUC": "L",      "AUC": "I",      "GUC": "V",
                "UUA": "L",      "CUA": "L",      "AUA": "I",      "GUA": "V",
                "UUG": "L",      "CUG": "L",      "AUG": "M",      "GUG": "V",
                "UCU": "S",      "CCU": "P",      "ACU": "T",      "GCU": "A",
                "UCC": "S",      "CCC": "P",      "ACC": "T",      "GCC": "A",
                "UCA": "S",      "CCA": "P",      "ACA": "T",      "GCA": "A",
                "UCG": "S",      "CCG": "P",      "ACG": "T",      "GCG": "A",
                "UAU": "Y",      "CAU": "H",      "AAU": "N",      "GAU": "D",
                "UAC": "Y",      "CAC": "H",      "AAC": "N",      "GAC": "D",
                "UAA": "Stop",   "CAA": "Q",      "AAA": "K",      "GAA": "E",
                "UAG": "Stop",   "CAG": "Q",      "AAG": "K",      "GAG": "E",
                "UGU": "C",      "CGU": "R",      "AGU": "S",      "GGU": "G",
                "UGC": "C",      "CGC": "R",      "AGC": "S",      "GGC": "G",
                "UGA": "Stop",   "CGA": "R",      "AGA": "R",      "GGA": "G",
                "UGG": "W",      "CGG": "R",      "AGG": "R",      "GGG": "G"}

def count_nucleotides(dna: str) -> list:
    counted = collections.Counter(dna)
    return [counted["A"], counted["C"], counted["G"], counted["T"]]

def transcription(dna: str) -> str:
    return dna.replace('T', 'U')

def complement(dna: str) -> str:
      c_replace = dna.replace('C', 'g')
      g_replace = c_replace.replace('G', 'c')
      a_replace = g_replace.replace('A', 't')
      t_replace = a_replace.replace('T', 'a')
      return t_replace[::-1].upper()

def fibonacci(n: int, k: int) -> int:
    fibonacci_cache = {}
    if n in fibonacci_cache:
        return fibonacci_cache[n]
     
    if n == 1:
        value = 1
    elif n == 2:
        value = 1
    elif n > 2:
        value = fibonacci(n-1, k) + k*fibonacci(n-2, k)

    fibonacci_cache[n] = value
    return value

def split_fasta(input: str) -> dict:
    split = input.split('>')[1:]
    return_dict = {}
    for entry in split:
        split_entry = entry.split('\n', 1)
        return_dict[split_entry[0]] = split_entry[1].replace('\n', '')
    return return_dict

def calculate_gc_content(input: str) -> float:
    total = len(input)
    nucleotides = count_nucleotides(input)
    return (nucleotides[1] + nucleotides[2])/total

def hamming_distance(str1: str, str2: str) -> int:
    return sum(string1 != string2 for string1, string2 in zip(str1, str2))

def mendellian_inheritance(k: int, m: int, n: int) -> float:
    # k = homozygous dominant
    # m = heterozygous
    # n = homozygous recessive
    # probability of dominant gene = sum(probability of match * probabilty match produces dominant gene)
    total = k + m + n
    dominant_dominant = (k/total)*((k-1)/(total-1))
    dominant_hetero = 2*(k/total)*(m/(total-1))
    dominant_recessive = 2*(k/total)*(n/(total-1))
    hetero_hetero = (m/total)*((m-1)/(total-1))
    hetero_recessive = 2*(m/total)*(n/(total-1))
    return 0.75*hetero_hetero + dominant_dominant + dominant_hetero + 0.5*hetero_recessive + dominant_recessive

def translation(rna: str) -> str:
    codons = [rna[i:i+3] for i in range(0, len(rna), 3)]
    amino_acids = [CODON_TABLE[codon] for codon in codons]
    return ''.join(amino_acids).split('Stop')[0]

def count_motif(dna: str, motif: str) -> list:
    positions = []
    index = dna.find(motif)
    while index != -1:
        positions.append(index+1)
        index = dna.find(motif, index + 1)
    return positions
   
def consensus_string(fasta: str) -> dict:
    dna_strings = split_fasta(fasta)
    dna_length = len([*dna_strings.values()][0])
    profile_dict = {"A": [0]*dna_length, "C": [0]*dna_length, "G": [0]*dna_length, "T": [0]*dna_length}
    for i in range(dna_length):
        for dna in dna_strings:
            profile_dict[dna_strings[dna][i]][i] += 1
    result_string = ""
    for i in range(dna_length):
        max_char = max(profile_dict, key = lambda k: profile_dict[k][i])
        result_string += max_char
    return {"consensus string": result_string, "matrix": profile_dict}

def mortal_fibonacci(n: int, m: int, cache=None) -> int:
    population = [0] * m
    population[0] = 1           #Base case, starting with 1 pair of rabbits

    for month in range(1, n):
        # Babies will be all of the adult rabbits minus the ones that will die this month
        babies = sum(population[1:])
        # Shift the ages over
        for age in range(m-1, 0, -1):
            population[age] = population[age-1]
        # Babies are now 1 month old
        population[0] = babies
    
    return sum(population)

def build_overlap_graph(fasta: str, k: int) -> list:
    adjacency_list = []
    fasta_dict = split_fasta(fasta)
    for s_label, s_entry in fasta_dict.items():
        s_suffix = s_entry[-k:]
        for t_label, t_entry in fasta_dict.items():
            t_prefix = t_entry[:k]
            if s_suffix == t_prefix and s_label != t_label:
                adjacency_list.append((s_label, t_label))
    return adjacency_list

