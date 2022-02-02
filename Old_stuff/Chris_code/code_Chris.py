# Load Gene sequence

with open('C:\\Users\\Christopher\\Desktop\\clb2.txt', 'r') as file:

    Gene = file.read().replace('\n', '')
# Clean away numbers and spaces
Gene=Gene.replace(' ', '') # spaces
from string import digits
remove_digits=str.maketrans('', '', digits)
Gene=Gene.translate(remove_digits)   # numbers

# Import Gene -+1000 base pairs - clean in the same way
with open('C:\\Users\\Christopher\\Desktop\\clb2 +- 1000.txt', 'r') as file:
    Gene_plus = file.read().replace('\n', '')

# Clean away numbers and spaces
Gene_plus=Gene_plus.replace(' ', '') # spaces
from string import digits
remove_digits=str.maketrans('', '', digits)
Gene_plus=Gene_plus.translate(remove_digits)   # numbers
Left_of_Gene=Gene_plus[0:1000] # 1000 bp's to the left of the Gene
Right_of_Gene=Gene_plus[int(len(Gene)):int(len(Gene))+1000] # 1000 bp's to the right of the Gene (includes 3' UTR segment)

# Importing Enzyme restriction-site Database

with open('C:\\Users\\Christopher\\Desktop\\Raw enzyme data.txt') as file:
    lines=file.readlines()
    database=[line.rstrip() for line in lines]  
# convert to pandas dataframe

import pandas as pd
import numpy as np
database=pd.DataFrame(database)  
# Sort through restriction site database to discard restriction sites with either non 'A', 'C', 'G', 'T' letters and/or are less than 6 base pairs
# split database into enzyme names and restriction sites
enzyme_names=database[::2]
db = database.iloc[1: , :]
restriction_sites=db[::2]
restrict_vect=[]

for i in range(len(restriction_sites)):
    if len(restriction_sites.iloc[i,0])<6 or len(restriction_sites.iloc[i,0])>8 or len(restriction_sites.iloc[i,0])==7:
        restrict_vect.append(i)

ind=list(set(restrict_vect))

# Now delete these indexed cut sites from just the restriction cut sites table
restriction_sites=pd.DataFrame(restriction_sites)
Viables=restriction_sites.drop(restriction_sites.index[[ind]])
fake_viables=[]



for i in range(len(Viables)):
    string_to_test = Viables.iloc[i,0]
    chars_to_check = ["R", "Y", "N", "W", "K", "M", "V", "H"]
    for char in chars_to_check:
        if char in string_to_test:
            fake_viables.append(i)


Cleaned_Enzymes=Viables.drop(Viables.index[[fake_viables]]) # List of all viable cut-site sequences
Cleaned_Enzymes=list(Cleaned_Enzymes.iloc[:,0])
Cleaned_Enzymes=set(Cleaned_Enzymes)
Cleaned_Enzymes=pd.DataFrame(Cleaned_Enzymes)

# We are now eliminating non-palindromic Enzymes - as the unique ones are unique on top and bottom strand in the plasmid
Palindromic_clean_enzymes=[]

for i in range(len(Cleaned_Enzymes)):
    # Creating reverse compliment of Enzyme cut
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = "".join(complement.get(base, base) for base in reversed(Cleaned_Enzymes.iloc[i,0]))

    if Cleaned_Enzymes.iloc[i,0]==reverse_complement:
        Palindromic_clean_enzymes.append(Cleaned_Enzymes.iloc[i,0])

       
# TODO
# For each Palindromic clean enzyme, we now must go from 100 bps from the end of the Gene: CUT || DISTANCE FROM END (count in bp's) - then re-order so its in ascending order or close-ness to 100 bp's end
Cuts_in_Gene=[]
distance_of_cut_from_end_of_gene=[]

for i in range(len(Palindromic_clean_enzymes)):
    for j in range(len(Gene)):
        count=len(Gene)-j # counting down the base pairs from the end of the gene len(gene) to 0.
        if Palindromic_clean_enzymes[i]==Gene[j:j+len(Palindromic_clean_enzymes[i])] and count>=100:
            Cuts_in_Gene.append(Palindromic_clean_enzymes[i])
            distance_of_cut_from_end_of_gene.append(count)

              
Cuts_in_Gene=pd.DataFrame(Cuts_in_Gene)
distance_of_cut_from_end_of_gene=pd.DataFrame(distance_of_cut_from_end_of_gene)
Table_of_cuts_in_Gene=pd.concat([Cuts_in_Gene, distance_of_cut_from_end_of_gene], axis=1)
Table_of_cuts_in_Gene.columns=['Cuts in Gene', 'Distance from the end of the Gene']   
Sorted_Table_of_cuts_in_Gene=Table_of_cuts_in_Gene.sort_values(by=['Distance from the end of the Gene'])
Sorted_Gene_cuts=list(Sorted_Table_of_cuts_in_Gene.iloc[:,0].values.tolist())
Sorted_Gene_cut_distances=list(Sorted_Table_of_cuts_in_Gene.iloc[:,1].values.tolist())

 
linker='GGTGCTTCTGTTGGTGCTTCTGTTTCTGTTGGTCCGC'
mCherry="atggtgagcaagggcgaggaggataacatggccatcatcaaggagttcatgcgcttcaaggtgcatatggagggctccgtgaacggccacgagttcgagatcgagggcgagggcgagggccgcccctacgagggcacccagaccgccaagctgaaggtgaccaagggtggccccctgcccttcgcctgggacatcctgtcccctcagttcatgtacggctccaaggcctacgtgaagcaccccgccgacatccccgactacttgaagctgtccttccccgagggcttcaagtgggagcgcgtgatgaacttcgaggacggcggcgtggtgaccgtgacccaggactcctccctgcaggacggcgagttcatctacaaggtgaagctgcgcggcaccaacttcccctccgacggccccgtaatgcagaagaagaccatgggctgggaggcctcctccgagcggatgtaccccgaggacggcgccctgaagggcgagatcaagcagaggctgaagctgaaggacggcggccactacgacgctgaggtcaagaccacctacaaggccaagaagcccgtgcagctgcccggcgcctacaacgtcaacatcaagttggacatcacctcccacaacgaggactacaccatcgtggaacagtacgaacgcgccgagggccgccactccaccggcggcatggacgagctgtactag"
mCherry=mCherry.upper()

# For the given gene cut, in order of distance from gene, calculate 2X, 3' UTR, and if the gene cut appears anywhere in the linker, RFP, 2X and 3' UTR
# Creates a table having the cuts, the 2X and the distance from the gene end

gene_increment=10 # amount of base pairs to the left of the gene cut (MUST BE GREATER THAN 8 BP'S)
UTR_increment=20 # UTR must be at least the amount fo the gene cut
Non_unique_gene_cuts=[]
Unique_Gene_cuts=[]
two_Xs=[]
UTRs=[]

for i in range(len(Sorted_Gene_cuts)):
    a=0
    b=0
    c=0
    d=0
    # Calculating 2X
    two_X=Gene[(len(Gene)-(Sorted_Gene_cut_distances[i]+gene_increment)):(len(Gene))]
    # Calculating 3' UTR
    UTR=Right_of_Gene[0:(len(two_X)+UTR_increment)]

    # Checking uniqueness of the cut in the RFP, linker, 2X and 3' UTR
    # Linker uniqueness:
    for j in range(len(linker)):
        if Sorted_Gene_cuts[i]==linker[j:j+len(Sorted_Gene_cuts[i])]:
            a=1
            Non_unique_gene_cuts.append(Sorted_Gene_cuts[i])
            break
    # Checking RFP uniqueness (mCherry here)
    for k in range(len(mCherry)):
        if Sorted_Gene_cuts[i]==mCherry[k:k+len(Sorted_Gene_cuts[i])]:
            b=1
            Non_unique_gene_cuts.append(Sorted_Gene_cuts[i])
            break  

    # Checking 2X
    for n in range(len(two_X)):
        if Sorted_Gene_cuts[i]==two_X[n:n+len(Sorted_Gene_cuts[i])]:
            c=c+1
            Non_unique_gene_cuts.append(Sorted_Gene_cuts[i])
            break

    # Checking 3' UTR
    for m in range(len(UTR)):
        if Sorted_Gene_cuts[i]==UTR[m:m+len(Sorted_Gene_cuts[i])]:
            d=1
            Non_unique_gene_cuts.append(Sorted_Gene_cuts[i])
            break

    if a==0 and b==0 and c==1 and d==0:
        two_Xs.append(two_X)
        UTRs.append(UTR)
        Unique_Gene_cuts.append(Sorted_Gene_cuts[i])

if not Unique_Gene_cuts:
    print("Error, no unique gene cuts found in using specified gene and utr increments - please re-adjust")

# this lets us fetch the 2X's of the unique gene cuts later in the program by referring to the gene cut itself
Gene_cuts_and_2_X=pd.concat([pd.DataFrame(Unique_Gene_cuts), pd.DataFrame(two_Xs)], axis=1)
Gene_cuts_and_2_X.columns=(["Gene cuts", "2X for Gene cut"]) 

# Now that we have gathered the unique Gene cuts in the insertion sequence we can begin to examine the plasmid
# import plasmid sequence

with open('C:\\Users\\Christopher\\Desktop\\pet_ura3.txt', 'r') as file:
    plasmid = file.read().replace('\n', '')
    plasmid = plasmid.upper()


# Search through plasmid for each unique Gene cut, however we should keep a list of both the cuts found in the plasmid and those not found in the plasmid
Unique_gene_cuts_in_plasmid=[]
# Checking if cuts are in plasmid
for i in range(len(Unique_Gene_cuts)):
    for j in range(len(plasmid)):
        if Unique_Gene_cuts[i]==plasmid[j:j+len(Unique_Gene_cuts[i])]:
            Unique_gene_cuts_in_plasmid.append(Unique_Gene_cuts[i])


# We need to keep it SIMPLE
# You have the unique gene cuts
# You have the unique cuts in the MCS
# Now, for each gene cut you:
    # 1. check its not in the backbone.
    # 2. check to see if its in the MCS.
    # 3. if its in the MCS - take the unique pairs of MCS cuts and check if they are in the insertion sequence of that gene cut in the MCS.
    # 4. when you find a pair that is unique in the insertion sequence, check if, in the MCS, the gene cut lies between the pair.
    # 5. if it does - we have the ends of the insertion sequence as that pair - TODO
    # 6. if it doesn't - we move onto the next pair - if no pairs satisfy that criteria then we cannot use the gene cut.


# 1. check gene cut is not in backbone
# here, the backbone is just the plasmid sequence without the MCS
backbone_1="TCGCGCGTTTCGGTGATGACGGTGAAAACCTCTGACACATGCAGCTCCcGGACTCTAATTTGTGAGTTTAGTATACATGCATTTACTTATAATACAGTTTTTTAGTTTTGCTGGCCGCATCTTCTCAAATATGCTTCCCAGCCTGCTTTTCTGTAACGTTCACCCTCTACCTTAGCATCCCTTCCCTTTGCAAATAGTCCTCTTCCAACAATAATAATGTCAGATCCTGTAGAGACCACATCATCCACGGTTCTATACTGTTGACCCAATGCGTCTCCCTTGTCATCTAAACCCACACCGGGTGTCATAATCAACCAATCGTAACCTTCATCTCTTCCACCCATGTCTCTTTGAGCAATAAAGCCGATAACAAAATCTTTGTCGCTCTTCGCAATGTCAACAGTACCCTTAGTATATTCTCCAGTAGATAGGGAGCCCTTGCATGACAATTCTGCTAACATCAAAAGGCCTCTAGGTTCCTTTGTTACTTCTTCTGCCGCCTGCTTCAAACCGCTAACAATACCTGGGCCCACCACACCGTGTGCATTCGTAATGTCTGCCCATTCTGCTATTCTGTATACACCCGCAGAGTACTGCAATTTGACTGTATTACCAATGTCAGCAAATTTTCTGTCTTCGAAGAGTAAAAAATTGTACTTGGCGGATAATGCCTTTAGCGGCTTAACTGTGCCCTCCATGGAAAAATCAGTCAAGATATCCACATGTGTTTTTAGTAAACAAATTTTGGGACCTAATGCTTCAACTAACTCCAGTAATTCCTTGGTGGTACGAACATCCAATGAAGCACACAAGTTTGTTTGCTTTTCGTGCATGATATTAAATAGCTTGGCAGCAACAGGACTAGGATGAGTAGCAGCACGTTCCTTATATGTAGCTTTCGACATGATTTATCTTCGTTTCCTGCAGGTTTTTGTTCTGTGCAGTTGGGTTAAGAATACTGGGCAATTTCATGTTTCTTCAACACTACATATGCGTATATATACCAATCTAAGTCTGTGCTCCTTCCTTCGTTCTTCCTTCTGTTCGGAGATTACCGAATCAAAAAAATTTCAAGGAAACCGAAATCAAAAAAAAGAATAAAAAAAAAATGATGAATTGAAAAGGTGGTATGGTGCACTCTCAGTACTCCgGGAGACGGTCACAGCTTGTCTGTAAGCGGATGCCGGGAGCAGACAAGCCCGTCAGGGCGCGTCAGCGGGTGTTGGCGGGTGTCGGGGCTGGCTTAACTATGCGGCATCAGAGCAGATTGTACTGAGAGTGCACCATATCGACTACGTCGTTAAGGCCGTTTCTGACAGAGTAAAATTCTTGAGGGAACTTTCACCATTATGGGAAATGGTTCAAGAAGGTATTGACTTAAACTCCATCAAATGGTCAGGTCATTGAGTGTTTTTTATTTGTTGTATTTTTTTTTTTTTAGAGAAAATCCTCCAATATATAAATTAGGAATCATAGTTTCATGATTTTCTGTTACACCTAACTTTTTGTGTGGTGCCCTCCTCCTTGTCAATATTAATGTTAAAGTGCAATTCTTTTTCCTTATCACGTTGAGCCATTAGTATCAATTTGCTTACCTGTATTCCTTTACATCCTCCTTTTTCTCCTTCTTGATAAATGTATGTAGATTGCGTATATAGTTTCGTCTACCCTATGAACATATTCCATTTTGTAATTTCGTGTCGTTTCTATTATGAATTTCATTTATAAAGTTTATGTACAAATATCATAAAAAAAGAGAATCTTTTTAAGCAAGGATTTTCTTAACTTCTTCGGCGACAGCATCACCGACTTCGGTGGTACTGTTGGAACCACCTAAATCACCAGTTCTGATACCTGCATCCAAAACCTTTTTAACTGCATCTTCAATGGCCTTACCTTCTTCAGGCAAGTTCAATGACAATTTCAACATCATTGCAGCAGACAAGATAGTGGCGATAGGGTCAACCTTATTCTTTGGCAAATCTGGAGCAGAACCGTGGCATGGTTCGTACAAACCAAATGCGGTGTTCTTGTCTGGCAAAGAGGCCAAGGACGCAGATGGCAACAAACCCAAGGAACCTGGGATAACGGAGGCTTCATCGGAGATGATATCACCAAACATGTTGCTGGTGATTATAATACCATTTAGGTGGGTTGGGTTCTTAACTAGGATCATGGCGGCAGAATCAATCAATTGATGTTGAACCTTCAATGTAGGGAATTCGTTCTTGATGGTTTCCTCCACAGTTTTTCTCCATAATCTTGAAGAGGCCAAAAGATTAGCTTTATCCAAGGACCAAATAGGCAATGGTGGCTCATGTTGTAGGGCCATGAAAGCGGCCATTCTTGTGATTCTTTGCACTTCTGGAACGGTGTATTGTTCACTATCCCAAGCGACACCATCACCATCGTCTTCCTTTCTCTTACCAAAGTAAATACCTCCCACTAATTCTCTGACAACAACGAAGTCAGTACCTTTAGCAAATTGTGGCTTGATTGGAGATAAGTCTAAAAGAGAGTCGGATGCAAAGTTACATGGTCTTAAGTTGGCGTACAATTGAAGTTCTTTACGGATTTTTAGTAAACCTTGTTCAGGTCTAACACTACCGGTACCCCATTTAGGACCACCCACAGCACCTAACAAAACGGCATCAGCCTTCTTGGAGGCTTCCAGCGCCTCATCTGGAAGTGGAACACCTGTAGCATCGATAGCAGCACCACCAATTAAATGATTTTCGAAATCGAACTTGACATTGGAACGAACATCAGAAATAGCTTTAAGAACCTTAATGGCTTCGGCTGTGATTTCTTGACCAACGTGGTCACCTGGCAAAACGACGATCTTCTTAGGGGCAGACATTAGAATGGTATATCCTTGAAATATATATATATATATTGCTGAAATGTAAAAGGTAAGAAAAGTTAGAAAGTAAGACGATTGCTAACCACCTATTGGAAAAAACAATAGGTCCTTAAATAATATTGTCAACTTCAAGTATTGTGATGCAAGCATTTAGTCATGAACGCTTCTCTATTCTATATGAAAAGCCGGTTCCGGCGCTCTCACCTTTCCTTTTTCTCCCAATTTTTCAGTTGAAAAAGGTATATGCGTCAGGCGACCTCTGAAATTAACAAAAAATTTCCAGTCATCGAATTTGATTCTGTGCGATAGCGCCCCTGTGTGTTCTCGTTATGTTGAGGAAAAAAATAATGGTTGCTAAGAGATTCGAACTCTTGCATCTTACGATACCTGAGTATTCCCACAGTTAACTGCGGTCAAGATATTTCTTGAATCAGGCGCCTTAGACCGCTCGGCCAAACAACCAATTACTTGTTGAGAAATAGAGTATAATTATCCTATAAATATAACGTTTTTGAACACACATGAACAAGGAAGTACAGGACAATTGATTTTGAAGAGAATGTGGATTTTGATGTAATTGTTGGGATTCCATTTTTAATAAGGCAATAATATTAGGTATGTAGATATACTAGAAGTTCTCCTCGACCGGTCGATATGCGGTGTGAAATACCGCACAGATGCGTAAGGAGAAAATACCGCATCAGGAAATTGTAAgCGTTAATATTTTGTTAAAATTCGCGTTAAATTTTTGTTAAATCAGCTCATTTTTTAACCAATAGGCCGAAATCGGCAAAATCCCTTATAAATCAAAAGAATAGACCGAGATAGGGTTGAGTGTTGTTCCAGTTTGGAACAAGAGTCCACTATTAAAGAACGTGGACTCCAACGTCAAAGGGCGAAAAACCGTCTATCAGGGCGATGGCCCACTACGTGAACCATCACCCTAATCAAGTTTTTTGGGGTCGAGGTGCCGTAAAGCACTAAATCGGAACCCTAAAGGGAGCCCCCGATTTAGAGCTTGACGGGGAAAGCCGGCGAACGTGGCGAGAAAGGAAGGGAAGAAAGCGAAAGGAGCGGGCGCTAGGGCGCTGGCAAGTGTAGCGGTCACGCTGCGCGTAACCACCACACCCGCCGCGCTTAATGCGCCGCTACAGGGCGCGTCCATTCGCCATTCAGGCTGCGCAACTGTTGGGAAGGGCGATCGGTGCGGGCCTCTTCGCTATTACGCCAGCTGGCGAAAGGGGGATGTGCTGCAAGGCGATTAAGTTGGGTAACGCCAGGGTTTTCCCAGTCACGACGTTGTAAAACGACGGCCAGTGAGCGCGC"
backbone_2="GCGCGCTTGGCGTAATCATGGTCATAGCTGTTTCCTGTGTGAAATTGTTATCCGCTCACAATTCCACACAACATAcGAGCCGGAAGCATAAAGTGTAAAGCCTGGGGTGCCTAATGAGTGAGcTAACTCACATTAATTGCGTTGCGCTCACTGCCCGCTTTCCAGTCGGGAAACCTGTCGTGCCAGCTGCATTAATGAATCGGCCAACGCGCGGGGAGAGGCGGTTTGCGTATTGGGCGCTCTTCCGCTTCCTCGCTCACTGACTCGCTGCGCTCGGTCGTTCGGCTGCGGCGAGCGGTATCAGCTCACTCAAAGGCGGTAATACGGTTATCCACAGAATCAGGGGATAACGCAGGAAAGAACATGTGAGCAAAAGGCCAGCAAAAGGCCAGGAACCGTAAAAAGGCCGCGTTGCTGGCGTTTTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGGACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGACCCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACATTTCCCCGAAAAGTGCCACCTGACGTCTAAGAAACCATTATTATCATGACATTAACCTATAAAAATAGGCGTATCACGAGGCCCTTTCGTC"
Unique_Gene_cuts_not_in_backbone=[]
for i in range(len(Unique_Gene_cuts)):
    counter=0
    for j in range(len(backbone_1)):
        if Unique_Gene_cuts[i]==backbone_1[j:j+len(Unique_Gene_cuts[i])]:
            print(" We cannot use Gene cut - found in backbone 1: " + Unique_Gene_cuts[i])
            counter=counter+1
    for k in range(len(backbone_2)):
        if Unique_Gene_cuts[i]==backbone_2[j:j+len(Unique_Gene_cuts[i])]:
            print(" We cannot use Gene cut - found in backbone 2: " + Unique_Gene_cuts[i])
            counter=counter+1      
    if counter==0:
        Unique_Gene_cuts_not_in_backbone.append(Unique_Gene_cuts[i])

if not Unique_Gene_cuts_not_in_backbone:
    print("Error, Gene cut found within backbone - please readjust parameters")
     
# 2 check to see if its in the MCS

MCS="GTAATACGACTCACTATAGGGCGAATTGGGTACCGGGCCCCCCCTCGAGGTCGACGGTATCGATAAGCTTGATATCGAATTCCTGCAGCCCGGGGGATCCACTAGTTCTAGAGCGGCCGCCACCGCGGTGGAGCTCCAGCTTTTGTTCCCTTTAGTGAGGGTTAATT"
Unique_Gene_cuts_only_in_mcs=[]
position_in_MCS_start=[]
position_in_MCS_end=[]   
for i in range(len(Unique_Gene_cuts_not_in_backbone)):
    for j in range(len(MCS)):
        if Unique_Gene_cuts_not_in_backbone[i]==MCS[j:j+len(Unique_Gene_cuts_not_in_backbone[i])]:
            Unique_Gene_cuts_only_in_mcs.append(Unique_Gene_cuts_not_in_backbone[i])
            position_in_MCS_start.append(j)
            position_in_MCS_end.append(j+len(Unique_Gene_cuts_not_in_backbone[i]))


# 3 we now have to check pairs of unique MCS cuts in the plasmid to see if they are unique in the insertion sequence of the gene cuts           
cuts_in_MCS=[]
MCS_cut_starts=[]
MCS_cut_end=[]
for i in range(len(Palindromic_clean_enzymes)):
    for j in range(len(MCS)):
        if Palindromic_clean_enzymes[i]==MCS[j:j+len(Palindromic_clean_enzymes[i])]:
            cuts_in_MCS.append(Palindromic_clean_enzymes[i])
            MCS_cut_starts.append(j)
            MCS_cut_end.append(j+len(Palindromic_clean_enzymes[i]))
             

cuts_in_MCS_data=pd.DataFrame(cuts_in_MCS)
MCS_cut_starts_data=pd.DataFrame(MCS_cut_starts)
MCS_cut_end_data=pd.DataFrame(MCS_cut_end)         
MCS_cut_table=pd.concat([cuts_in_MCS_data, MCS_cut_starts_data, MCS_cut_end_data], axis=1)
MCS_cut_table.columns=["MCS cuts", "Cut Start", "Cut End"]     

 
# Just searches for all cuts that are unique in MCS
plasmid_cuts=[]
for i in range(len(Palindromic_clean_enzymes)):
    for j in range(len(plasmid)):
        if Palindromic_clean_enzymes[i]==plasmid[j:j+len(Palindromic_clean_enzymes[i])]:
            plasmid_cuts.append(Palindromic_clean_enzymes[i])
            

import pandas as pd
plasmid_cuts=pd.DataFrame(plasmid_cuts)           
plasmid_counts=plasmid_cuts.iloc[:,0].value_counts()
unique_plasmid_cuts=list(plasmid_counts[plasmid_counts==1].index[:])
#all_plasmid_cuts=list(set(plasmid_cuts.iloc[:,0])) # multiple or unique - these are the cuts in the plasmid
# Want to find the cuts in the MCS that are unique in the plasmid,                
Unique_cuts_in_MCS=list(set(unique_plasmid_cuts).intersection(set(cuts_in_MCS)))

# We now have a table of all of the unique cuts in the MCS with their start and end indexes
Table_of_unique_MCS_cuts=MCS_cut_table[MCS_cut_table['MCS cuts'].isin(Unique_cuts_in_MCS)]
Sorted_table_of_unique_MCS_cuts=Table_of_unique_MCS_cuts.sort_values(by=['Cut Start'])
# We now know the order along the MCS of the unique cuts - any gene cuts we are to use shall also be in here, as if not then we would not be considering those gene cuts as viable.
# Now just make sure we have the gene cuts in the sorted MCS unique cuts identified to us
Sorted_table_of_unique_MCS_gene_cuts=Sorted_table_of_unique_MCS_cuts[Sorted_table_of_unique_MCS_cuts["MCS cuts"].isin(Unique_Gene_cuts_only_in_mcs)]
Sorted_table_of_non_gene_unique_cuts_in_MCS=Sorted_table_of_unique_MCS_cuts[~Sorted_table_of_unique_MCS_cuts["MCS cuts"].isin(Unique_Gene_cuts_only_in_mcs)]

# 4 going through each gene cut in the MCS - then finding valid pairs that have the gene cut lie in between them, then searching the insertion sequence to verify if pair is unique
Sorted_MCS_cut_sequences=Sorted_table_of_unique_MCS_cuts.iloc[:,0].values.tolist()
# We create a table here with structure: Gene cut in mcs | list of unique MCS cuts before that gene cut | list of unique MCS cuts after that gene cut
# (as the elements the table are created, create a single list representing each row with elements [ gene cut, [list of cuts before gene cut], [list of cuts after gene cut] ])
def listToStringWithoutBrackets(list1): # function is needed to remove square brackets from cut string
    return str(list1).replace('[','').replace(']','')

INFO_ON_MCS_CUTS_FOR_GENE_ROW=[]
for k in range(len(Unique_Gene_cuts_only_in_mcs)): # for all the gene cuts in the MCS
    for i in range(len(Sorted_MCS_cut_sequences)): # for all sorted cuts in the MCS
        if Sorted_MCS_cut_sequences[i]==listToStringWithoutBrackets(Unique_Gene_cuts_only_in_mcs[k]): # if the gene cut is found in the MCS
            beginning_sorted_cuts=Sorted_MCS_cut_sequences[0:i]    # take the cuts before
            end_sorted_cuts=Sorted_MCS_cut_sequences[i+1:len(Sorted_MCS_cut_sequences)] # and the cuts after the gene cut
            # We now sort through to find the unique sorted non gene cuts in the insertation sequence of the gene cut
            useable_beginning_MCS_cuts_for_gene_cut=[]
            for j in range(len(beginning_sorted_cuts)): # Going through the cuts before the gene cut in the MCS
                ########################################
                a=0
                b=0
                c=0
                d=0
                # Calculating 2X for the gene cut in the MCS
                print("Searching through beginning cuts")
                two_X=Gene[(len(Gene)-(int(Table_of_cuts_in_Gene.loc[Table_of_cuts_in_Gene['Cuts in Gene'] == Unique_Gene_cuts_only_in_mcs[k]].iloc[:,1])+gene_increment)):(len(Gene))]
                # Calculating 3' UTR from that MCS 2X
                UTR=Right_of_Gene[0:(len(two_X)+UTR_increment)]
                
                if len(two_X) > len(Right_of_Gene):
                    print("WARNING - 2X is longer than 1000 bp's to the right of gene")
                    #break

                # Checking uniqueness of the cut in the RFP, linker, 2X and 3' UTR
                # Linker uniqueness:
                for p in range(len(linker)):
                    if beginning_sorted_cuts[j]==linker[p:p+len(beginning_sorted_cuts[j])]:
                        a=1
                        break

                # Checking RFP uniqueness (mCherry here):
                for h in range(len(mCherry)):
                    if beginning_sorted_cuts[j]==mCherry[h:h+len(beginning_sorted_cuts[j])]:
                        b=1
                        break

                # Checking 2X
                for n in range(len(two_X)):
                    if beginning_sorted_cuts[j]==two_X[n:n+len(beginning_sorted_cuts[j])]:
                        c=1   
                        break

                # Checking 3' UTR
                for m in range(len(UTR)):
                    if beginning_sorted_cuts[j]==UTR[m:m+len(beginning_sorted_cuts[j])]:
                        d=1
                        break

                if a==0 and b==0 and c==0 and d==0:
                    useable_beginning_MCS_cuts_for_gene_cut.append(beginning_sorted_cuts[j])
                    # insert this list into the appropriate spot in the table
                    
                # We now do exactly the same, but this time for the end sorted cuts in the MCS
                useable_end_MCS_cuts_for_gene_cut=[]
                for j in range(len(end_sorted_cuts)): # Going through the cuts before the gene cut in the MCS
                    a=0
                    b=0
                    c=0
                    d=0
                    print("Searching through potential end cuts")

                    # Calculating 2X for the gene cut in the MCS
                    two_X=Gene[(len(Gene)-(int(Table_of_cuts_in_Gene.loc[Table_of_cuts_in_Gene['Cuts in Gene'] == Unique_Gene_cuts_only_in_mcs[k]].iloc[:,1]+gene_increment))):(len(Gene))]
                    
                    # Calculating 3' UTR from that MCS 2X
                    UTR=Right_of_Gene[0:(len(two_X)+UTR_increment)]

                    # Checking uniqueness of the cut in the RFP, linker, 2X and 3' UTR
                    # Linker uniqueness:
                    for p in range(len(linker)):

                        if end_sorted_cuts[j]==linker[p:p+len(end_sorted_cuts[j])]:
                            a=1
                            break

                    # Checking RFP uniqueness (mCherry here):
                    for h in range(len(mCherry)):
                        if end_sorted_cuts[j]==mCherry[h:h+len(end_sorted_cuts[j])]:
                            b=1 
                            break

                    # Checking 2X
                    for n in range(len(two_X)):
                        if end_sorted_cuts[j]==two_X[n:n+len(end_sorted_cuts[j])]:
                            c=1 
                            break
                       
                    # Checking 3' UTR
                    for m in range(len(UTR)):
                        if end_sorted_cuts[j]==UTR[m:m+len(end_sorted_cuts[j])]:
                            d=1    
                            break

                    if a==0 and b==0 and c==0 and d==0:
                        useable_end_MCS_cuts_for_gene_cut.append(end_sorted_cuts[j])
                    #########################################

    # can just find the last element in the ends and the first element in the beginning cuts - creates largest bracket around gene cut in MCS
    INFO_ON_MCS_CUTS_FOR_GENE_ROW.append([Unique_Gene_cuts_only_in_mcs[k], useable_beginning_MCS_cuts_for_gene_cut, useable_end_MCS_cuts_for_gene_cut])       
            
# Create dataframe of the lists stacked atop one another
MCS_Gene_cut_DataFrame=pd.DataFrame(INFO_ON_MCS_CUTS_FOR_GENE_ROW)
# After the loop, we can now use the table to find the cuts we are to use for the end points of the insertion sequence, and the cuts to use in between the linker, RFP, 2X and 3' UTR
# lets makes the output be the insertion sequence for each gene cut, we can make this a for loop process

insertion_sequences=[]

for i in range(len(Unique_Gene_cuts_only_in_mcs)):
    MCS_partition_start=MCS_Gene_cut_DataFrame.iloc[i, 1][0]
    useable_beginning_MCS_cuts_for_gene_cut.remove(MCS_partition_start) 
    MCS_partition_end=MCS_Gene_cut_DataFrame.iloc[i, 2][len(useable_end_MCS_cuts_for_gene_cut)-1]
    useable_end_MCS_cuts_for_gene_cut.remove(MCS_partition_end) 
    useable_insertion_sequence_cuts=useable_beginning_MCS_cuts_for_gene_cut+useable_end_MCS_cuts_for_gene_cut
    # Now building the insertion sequence
    insertion_start=MCS_partition_start
    insertion_2X=Gene_cuts_and_2_X[Gene_cuts_and_2_X["Gene cuts"]==Unique_Gene_cuts_only_in_mcs[i]].iloc[:,1]
    insertion_cut_1=useable_insertion_sequence_cuts[1]
    insertion_linker=linker
    insertion_cut_2=useable_insertion_sequence_cuts[2]
    insertion_RFP=mCherry
    insertion_cut_3=useable_insertion_sequence_cuts[3]
    insertion_cut_UTR=Right_of_Gene[0:(UTR_increment+len(two_X))]
    insertion_end=MCS_partition_end
    # Puts each insetion sequence into a list of insertion sequences. Direction is 5' -> 3'.
    insertion_sequence=insertion_start+insertion_2X+insertion_cut_1+insertion_linker+insertion_cut_2+insertion_RFP+insertion_cut_3+insertion_cut_UTR+insertion_end
    print("Insertion Sequence for Gene cut: " + Unique_Gene_cuts_only_in_mcs[i] +  " with start as: " + insertion_start + " and end as: " + insertion_end + " has between segment cuts of: " + insertion_cut_1 + " " + insertion_cut_2 + " " + insertion_cut_3)
    print("Insertion Start: ")
    print(insertion_start)
    print("Insertion 2X: ")
    print(list(insertion_2X))
    print("Cut between 2X and linker: ")
    print(insertion_cut_1)
    print("Linker: ")
    print(insertion_linker)
    print("Cut between linker and RFP: ")
    print(insertion_cut_2)
    print("RFP: ")
    print(insertion_RFP)
    print("Cut between RFP and 3' UTR segment: ")
    print(insertion_cut_3)
    print("3' UTR segment: ")
    print(insertion_cut_UTR)
    print("Insertion end")
    print(insertion_end)
