import re
import os
from collections import defaultdict

import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner


def read_from_fsa(fsa_file_path):
    """Reads a FASTA file and returns the description (header) and the sequence"""
    if not fsa_file_path:
        return '',''
    try:
        with open(fsa_file_path) as f:
            name = f.readline()[1:-1]
            sequence = f.read().replace('\n', '')
            print(name)
            print(sequence)
            print('Length:', len(sequence))
            return name, sequence.upper()
    except Exception as e:
        print('Unable to read:', fsa_file_path,e)

def read_enzyme_list(enzymefile):
    """Reads in and processes the enzyme list"""
    with open(enzymefile) as file:
        lines=file.readlines()
        enzymes=[line.rstrip() for line in lines]
        
    rsitelist = []
    enamelist = []
    if (len(enzymes) % 2) == 1:
        print('Enzyme list has an odd number of lines. Enzyme list should be a list of enzyme names and restriction sites in each line.')
    else:
        for linei in range(1, len(enzymes), 2):
            rsite = enzymes[linei]
            ename = enzymes[linei-1]
            # check length of restriction site
            if len(rsite) == 6 or len(rsite) == 8:
                # check whether any non-ACGT characters
                if len(re.sub('[ACGT]', '', rsite)) == 0:
                    #check whether palindromic so can focus on only one strand
                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    reverse_complement = "".join(complement.get(base, base) for base in reversed(rsite))
                    if rsite==reverse_complement:
                        # add name only if enzyme already in the list
                        rsiteexists=0 
                        for linej in range(len(rsitelist)):
                            if rsitelist[linej] == rsite:
                                rsiteexists=1
                                enamelist[linej] = enamelist[linej]+', '+ename
                                break
                        #add new entry if not yet there
                        if rsiteexists==0:
                            rsitelist.append(rsite)
                            enamelist.append(ename)
    return (np.array(rsitelist), np.array(enamelist))

def read_enzyme_list(enzymefile):
    """Reads in and processes the enzyme list"""
    with open(enzymefile) as file:
        lines=file.readlines()
        enzymes=[line.rstrip() for line in lines]
        
    rsitelist = []
    enamelist = []
    if (len(enzymes) % 2) == 1:
        print('Enzyme list has an odd number of lines. Enzyme list should be a list of enzyme names and restriction sites in each line.')
    else:
        for linei in range(1, len(enzymes), 2):
            rsite = enzymes[linei]
            ename = enzymes[linei-1]
            # check length of restriction site
            if len(rsite) == 6 or len(rsite) == 8:
                # check whether any non-ACGT characters
                if len(re.sub('[ACGT]', '', rsite)) == 0:
                    #check whether palindromic so can focus on only one strand
                    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
                    reverse_complement = "".join(complement.get(base, base) for base in reversed(rsite))
                    if rsite==reverse_complement:
                        # add name only if enzyme already in the list
                        rsiteexists=0 
                        for linej in range(len(rsitelist)):
                            if rsitelist[linej] == rsite:
                                rsiteexists=1
                                enamelist[linej] = enamelist[linej]+', '+ename
                                break
                        #add new entry if not yet there
                        if rsiteexists==0:
                            rsitelist.append(rsite)
                            enamelist.append(ename)
    return (np.array(rsitelist), np.array(enamelist))

def read_gene_plus_string(Gene_plus):
    """Splits gene sequence into 3 pieces: ORF, 1000 bp upstream ("Left_of_gene") and 1000 bp downstream ("Right_of_gene")"""
    # # Import Gene -+1000 base pairs - clean in the same way
    # with open(genefile, 'r') as file:
    #     Gene_plus = file.read().replace('\n', '')

    # Clean away numbers and spaces
    Gene_plus = Gene_plus.replace(' ', '') # spaces
    remove_digits = str.maketrans('', '', '0123456789')
    Gene_plus = Gene_plus.translate(remove_digits)   # numbers
    Left_of_gene = Gene_plus[0:1000] # 1000 bp's to the left of the Gene
    Right_of_gene = Gene_plus[(len(Gene_plus)-1000):] # 1000 bp's to the right of the Gene (includes 3' UTR segment)
    Gene = Gene_plus[1000:(len(Gene_plus)-1000)]
    
    return (Left_of_gene, Gene, Right_of_gene)

def dna_to_protein(dna):
     """Translates DNA into a Protein. It truncates the 3' tail that doesn't make a full codon.
     \n STOP is denoted by '*'"""
     
     genetic_code = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }
     protein = ''
    
     if(len(dna)%3 != 0):
        dna = dna[: -(len(dna)%3)]

     for i in range(0, len(dna), 3):
         code = dna[i:i+3]
         if(code in genetic_code.keys()):
             protein += genetic_code[code]
         else:
             protein = ''
             break

     return protein

def read_overlapping_genes(path):
  """Reads lines from Overlapping_genes.txt."""
  overlapping_genes = defaultdict(list)
  try:
    with open(path, 'r') as f:
      lines = f.readlines()
      for line in lines:
        if "Chromosome" in line:
          continue
        line = line.strip()
        split = line.split()
        gene_1 = f"{split[3]} {split[4]}"
        gene_2 = f"{split[8]} {split[9]}"
        overlapping_genes[gene_1].append(gene_2)
        overlapping_genes[gene_2].append(gene_1)
    return overlapping_genes
  except FileNotFoundError:
    print(f"File '{path}' not found.")
    return {}
  except Exception as e:
    print(f"An error occurred: {e}")
    return {}

def get_full_gene_names_mapping():
    mapping = {}
    for filename in os.listdir('/content/individual_genes/'):
      if not filename.endswith('.fasta') and not filename.endswith('.fsa') and not filename.endswith('.fa'):
        continue
      firsthalf, secondhalf = filename.replace('.fasta','').split('_', 1)
      mapping[firsthalf] = f"{firsthalf} {secondhalf}"
      mapping[secondhalf] = f"{firsthalf} {secondhalf}"
    return mapping

def read_fasta_files(directory):
     # List all fasta files in the directory
     fasta_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".fasta") or f.endswith(".fa")]
     sequences = []

     # Iterate over each fasta file
     for fasta_file in fasta_files:
         #print(f"Processing file: {fasta_file}")  # Log the file being processed
         # Parse the sequences in the fasta file
         records = list(SeqIO.parse(fasta_file, "fasta"))

         # Check if the fasta file contains exactly one sequence
         if len(records) == 1:
             record = records[0]  # Get the single sequence
             sequence_str = str(record.seq)
             # Only keep sequences of length >= 2000
             if len(sequence_str) >= 2000:
                 #print(f"Adding sequence {record.id} of length {len(sequence_str)}")  # Log the sequence being added
                 sequences.append((record.id, sequence_str))

     return sequences

def align_query_to_all(query_seq, sequences):
     best_match = None
     best_score = -float('inf')
     aligner = PairwiseAligner()
     aligner.mode = "global" #faster than local

     # Align query sequence with each sequence
     for seq_id, seq in sequences:
         # Check if the absolute difference in length is less than 100 bp

         if abs(len(query_seq) - len(seq)) < 50:
             alignments = aligner.align(query_seq, seq)
             score = alignments.score
             if score > best_score:
                 best_score = score
                 best_match = (seq_id, seq, score)

     return best_match
