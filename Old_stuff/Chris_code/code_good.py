import sys
import numpy as np
import re




# Read in the gene file +/- 1000 base pairs and break it up into 3 pieces
def read_gene_string(genefile):

    # Import Gene -+1000 base pairs - clean in the same way
    with open(genefile, 'r') as file:
        Gene_plus = file.read().replace('\n', '')

    # Clean away numbers and spaces
    Gene_plus=Gene_plus.replace(' ', '') # spaces
    remove_digits=str.maketrans('', '', '0123456789')
    Gene_plus=Gene_plus.translate(remove_digits)   # numbers
    Left_of_gene=Gene_plus[0:1000] # 1000 bp's to the left of the Gene
    Right_of_gene=Gene_plus[(len(Gene_plus)-1000):] # 1000 bp's to the right of the Gene (includes 3' UTR segment)
    Gene=Gene_plus[1000:(len(Gene_plus)-1000)]
    
    return (Left_of_gene, Gene, Right_of_gene)




# Read in and process the enzyme list
def read_enzyme_list(enzymefile):

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
    return (rsitelist, enamelist)
        
        
        
        
# Main method
if __name__ == "__main__":

    genefile   = sys.argv[1]
    enzymefile = sys.argv[2]
    minhomology = 100
    
    Left_of_gene, Gene, Right_of_gene = read_gene_string(genefile)
    
    rsitelist, enamelist              = read_enzyme_list(enzymefile)
    
    
