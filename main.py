import argparse
import re
import numpy as np


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


def find_rsite_locations(sequence, rsitelist, enamelist, starting_end, min_homology=0):
    """Searches for cutsites from "rsitelist" in a given sequence. 
    Starts by ommiting "min_homology" bp from the "starting_end" (either 5 or 3 as in 5' and 3', respectively)"""
    
    sequence_cut = sequence[min_homology:len(sequence)-min_homology]
    
    if starting_end == 5:
        rsite_position_list = np.array([sequence_cut.find(rsite) for rsite in rsitelist])
        sorted_inds = np.argsort(rsite_position_list)

    elif starting_end == 3:
        rsite_position_list = np.array([sequence_cut.rfind(rsite) for rsite in rsitelist]) 
        sorted_inds = np.argsort(-rsite_position_list)

    # TODO Is sorting really needed? as we later sort sequences by the length?

    # remove not matched rsites
    rsite_position_list = rsite_position_list[sorted_inds]
    matched_inds = (rsite_position_list != -1)
    rsite_position_list = rsite_position_list[matched_inds]
    rsite_position_list += min_homology

    # apply sort and remove inds to rsitelist and enamelist
    rsitelist = rsitelist[sorted_inds][matched_inds]
    enamelist = enamelist[sorted_inds][matched_inds]

    return rsitelist, enamelist, rsite_position_list


def generate_start_end_sequences(left_chunk, right_chunk, rsite_side, rsite_position, rsite, minhomology, alpha, stop_codon_offset=0):
    """Generates X and alphaX sequences for the cases of tagging genes (at either 5' or 3' end)

    Splitted in two cases, depending on the position of the cutsite:
    rsite_side argument can have 'left' or 'right' values. 
    For 5' UTR tagging, 'left' means inside the gene starting from 5' end, 'right' means inside the right 1000 bp
    For 3' UTR tagging, 'left' means inside the left 1000 bp, 'right' means inside the gene starting from 3' end
    For gene deletion, 'left' means inside the left 1000 bp, 'right' means inside the right 1000 bp
    stop_codon_offset is used to remove the stop codon of the gene in case we are tagging the 3' end """
    

    if rsite_side not in ['left','right']:
        raise ValueError('Wrong riste_side value')

    if rsite_side == 'left':
        X = left_chunk[rsite_position-minhomology:len(left_chunk)-stop_codon_offset]
        alphaX = right_chunk[:int(alpha*len(X))]
        return X, alphaX
    else:
        X = right_chunk[:rsite_position+minhomology]
        alphaX = left_chunk[-int(alpha*len(X)):len(left_chunk)-stop_codon_offset]
        return alphaX, X


def rsite_search(Gene, rsitelist, enamelist, modality, alpha, min_homology=0,
                 left_of_Gene='', right_of_Gene='', FPGs=[], linker=''):
    """This is the central function that finds suitable cutsites for the given list of FPGs and the linker


    Modality is used to define the application of the program:
    5 for 5' tagging with given list of FPGs,
    3 for 3' taggging with given list of FPGs, and
    0 for deletion (FPGs neglected)

    Function returns a dictionary of parameters that describe the cloning """

    if modality not in [0, 3, 5]:
        raise ValueError('Wrong starting_position value')

    #5' tagging
    if modality == 5:

        Gene_and_right_of_Gene = Gene + right_of_Gene

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = find_rsite_locations(
            Gene_and_right_of_Gene, rsitelist, enamelist, 5, min_homology)
        rsitelist_5UTR, enamelist_5UTR, rsite_position_list_5TR = find_rsite_locations(
            left_of_Gene, rsitelist, enamelist, 3, min_homology)

        start_end_sequences_gene = [generate_start_end_sequences(
            left_of_Gene, Gene_and_right_of_Gene, 'right', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)]
        start_end_sequences_5UTR = [generate_start_end_sequences(
            left_of_Gene, Gene_and_right_of_Gene, 'left', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_5TR, rsitelist_5UTR)]

        full_sequences = [
            [start_seq+FPG[:-3]+linker+end_seq for FPG in FPGs] #We don't take the stop codon from FPG
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_5UTR
        ]
        real_alphas = [len(alphaX)/float(len(X)) for alphaX,X in start_end_sequences_gene] + \
            [len(alphaX)/float(len(X)) for X,alphaX in start_end_sequences_5UTR]

        return_dict = {
            'rsitelist': np.append(rsitelist_gene, rsitelist_5UTR),
            'enamelist': np.append(enamelist_gene, enamelist_5UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_5TR),
            'rsite_places': ['rsite in gene+3`UTR']*len(rsitelist_gene)+['rsite in 5`UTR']*len(rsitelist_5UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_5UTR,
            'real_alphas': real_alphas
        }

    # 3' tagging 
    elif modality == 3:

        left_of_Gene_and_Gene = left_of_Gene + Gene

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = find_rsite_locations(
            left_of_Gene_and_Gene, rsitelist, enamelist, 3, min_homology)
        rsitelist_3UTR, enamelist_3UTR, rsite_position_list_3UTR = find_rsite_locations(
            right_of_Gene, rsitelist, enamelist, 5, min_homology)

        start_end_sequences_gene = [
            generate_start_end_sequences(left_of_Gene_and_Gene, right_of_Gene, 'left', rsite_pos, rsite, min_homology, alpha, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)
        ]
        start_end_sequences_3UTR = [
            generate_start_end_sequences(left_of_Gene_and_Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, alpha, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_3UTR, rsitelist_3UTR)
        ]

        full_sequences = [
            [start_seq+linker+FPG+end_seq for FPG in FPGs]
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_3UTR
        ]
        real_alphas = [len(alphaX)/float(len(X)) for X,alphaX in start_end_sequences_gene] + \
            [len(alphaX)/float(len(X)) for alphaX,X in start_end_sequences_3UTR]

        return_dict = {
            'rsitelist': np.append(rsitelist_gene, rsitelist_3UTR),
            'enamelist': np.append(enamelist_gene, enamelist_3UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_3UTR),
            'rsite_places': ['rsite in 5`UTR+gene']*len(rsitelist_gene)+['rsite in 3`UTR']*len(rsitelist_3UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_3UTR,
            'real_alphas': real_alphas
        }

    # Deletion
    else:
        rsitelist_5, enamelist_5, rsite_position_list_5 = find_rsite_locations(left_of_Gene, rsitelist,
                                                                            enamelist, 3, min_homology)
        rsitelist_3, enamelist_3, rsite_position_list_3 = find_rsite_locations(right_of_Gene, rsitelist,
                                                                            enamelist, 5, min_homology)

        start_end_sequences_5UTR = [generate_start_end_sequences(
            left_of_Gene, right_of_Gene, 'left', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_5, rsitelist_5)]
        start_end_sequences_3UTR = [generate_start_end_sequences(
            left_of_Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_3, rsitelist_3)]

        full_sequences = [
            [start_seq+end_seq]
            for start_seq, end_seq in start_end_sequences_5UTR+start_end_sequences_3UTR
        ]
        real_alphas = [len(alphaX)/len(X) for X,alphaX in start_end_sequences_5UTR] + \
            [len(alphaX)/len(X) for alphaX,X in start_end_sequences_3UTR]

        return_dict = {
            'rsitelist': np.append(rsitelist_5, rsitelist_3),
            'enamelist': np.append(enamelist_5, enamelist_3),
            'rsite_position_list': np.append(rsite_position_list_5, rsite_position_list_3),
            'rsite_places': ['rsite in 5`']*len(rsitelist_5)+['rsite in 3`']*len(rsitelist_3),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_5UTR+start_end_sequences_3UTR,
            'real_alphas': real_alphas
        }

    full_sequences_per_rsite = return_dict['full_sequences']
    sort = np.argsort([len(x[0]) for x in full_sequences_per_rsite])
    for key, value in return_dict.items():
        return_dict[key] = list(np.array(value)[sort])

    return return_dict


def find_compatible_MCS_rsites(MCS, rsitelist, enamelist, full_sequences, backbone_no_MCS_5, backbone_no_MCS_3):
    """Searches for the cutsites that can be used for opening up the bakcbone and putting in the insert. 
    \n Performs two searches: one starting from the 5' end and one starting from the 3' end. """
    rsite1 = ''
    ename1 = ''
    rsite2 = ''
    ename2 = ''

    backbone_no_MCS = backbone_no_MCS_3 + backbone_no_MCS_5

    rsitelist_sorted_right, enamelist_sorted_right, _ = find_rsite_locations(MCS, rsitelist, enamelist, starting_end=5)
    for i, rsite in enumerate(rsitelist_sorted_right):
        if not any(rsite in full_sequence for full_sequence in full_sequences) and not rsite in backbone_no_MCS:
            rsite1 = rsite
            ename1 = enamelist_sorted_right[i]
            break

    rsitelist_sorted_left, enamelist_sorted_left, _ = find_rsite_locations(MCS, rsitelist, enamelist, starting_end=3)
    for i, rsite in enumerate(rsitelist_sorted_left):
        if not any(rsite in full_sequence for full_sequence in full_sequences) and not (rsite in backbone_no_MCS and rsite == rsite1):
            rsite2 = rsite
            ename2 = enamelist_sorted_left[i]
            break
        
    return rsite1, ename1, rsite2, ename2

def  find_additional_cutsites(plasmids, rsitelist, enamelist):
    """Finds all appropropriate cutsites that could be added between different pieces of the insert in the plasmid (e.g. between the gene and the linker or between the linker and the FPG etc).
    \n A good cutsite has to fullfill three criteria:
    1. Divisible by 3 (this could be circuimvented by accomidating the linker length)
    2. Does not cut the final plasmids (plural in case of several FPGs)
    3. Does not introduce a stop codon"""

    good_rsite_list = []
    good_enzyme_list = []
    for rsite, enzyme in zip(rsitelist, enamelist):
        if len(rsite)%3 == 0 and dna_to_protein(rsite).count('*') == 0:
            if not any(rsite in plasmid for plasmid in plasmids):
                good_rsite_list.append(rsite)
                good_enzyme_list.append(enzyme)

    return good_enzyme_list, good_rsite_list


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

def assemble_plasmid(backbone_no_MCS_5, backbone_no_MCS_3, sequence): 
    return backbone_no_MCS_5 + sequence + backbone_no_MCS_3


def main(args):

    # Read the files (1. and 2.)
    print('\nPIPOline by Stojkovic, Gligorovski, Rahi\n\n********** Reading input files\n\nVector backbone:')
    _, backbone = read_from_fsa(args.backbone_path)
    
    print('\nMultiple cloning site extracted from vector backbone:')
    MCS = backbone[args.MCS_start_ind-1:args.MCS_end_ind]
    print(MCS,'\nLength:', len(MCS))
    backbone_no_MCS_5 = backbone[0:args.MCS_start_ind-1]
    backbone_no_MCS_3 = backbone[args.MCS_end_ind:]

    print('\nLinker:')
    _, linker = read_from_fsa(args.linker_path)
    
    print('\nGene sequence +/- 1000 bps:')
    _, Gene_plus = read_from_fsa(args.Gene_path)
    left_of_Gene, Gene, right_of_Gene = read_gene_plus_string(Gene_plus)

    print('\nFluorescent protein genes:')
    FPGs = [read_from_fsa(path)[1] for path in args.FPG_paths] if len(args.FPG_paths) else []

    rsitelist, enamelist = read_enzyme_list(args.enzyme_path)
    
    if(args.modality == 5 or args.modality == 3):
        popular_rsitelist, popular_enamelist = read_enzyme_list(args.popular_enzyme_path)

    # 3.
    rsite_dict = rsite_search(Gene, rsitelist, enamelist, args.modality, args.alpha,
                              args.min_homology, left_of_Gene, right_of_Gene, FPGs, linker)

    gene_rsitelist_sorted = rsite_dict['rsitelist']
    gene_enamelist_sorted = rsite_dict['enamelist'] 
    gene_rsite_position_list_sorted = rsite_dict['rsite_position_list'] 
    rsite_places = rsite_dict['rsite_places'] 
    full_sequences_per_rsite = rsite_dict['full_sequences']
    start_end_sequences = rsite_dict['start_end_sequences']
    real_alphas = rsite_dict['real_alphas']
    print('\n',gene_rsitelist_sorted, gene_enamelist_sorted, gene_rsite_position_list_sorted, rsite_places,'\n')

    compatible_restriction_sites = []
    optimal_plasmid = ''
    MCS_rsites = []

    print("\n\n##################################################################\n\n")

    for i in range(len(gene_rsitelist_sorted)):

        # 4. and 5.
        rsite0 = gene_rsitelist_sorted[i]
        ename0 = gene_enamelist_sorted[i]
        rsite_place = rsite_places[i]
        full_sequences = full_sequences_per_rsite[i]
        start_seq, end_seq = start_end_sequences[i]
        real_alpha = real_alphas[i]
        print(str(i+1), rsite0, ename0, 'location:',rsite_place)
            
        # check if all full sequences have only one rsite0  
        if any([(full_sequence.count(rsite0) != 1) for full_sequence in full_sequences]):
            print("""Enzyme {} cannot satisfy the conditions.
            Restriction site {} is not unique in the assembled insert sequence.""".format(ename0, rsite0))
            print("##################################################################")
            
            start_name = '5` UTR' if args.modality in [5,0] else 'Gene end'
            end_name = '3` UTR' if args.modality in [3,0] else 'Gene start'
            print("\nPieces that should be used for cloning are:\n{}: {}\n{}: {}".format(
                start_name, start_seq, end_name, end_seq
            ))
            print("\nTotal length of gene chunks that should be subcloned is", str((len(start_seq)+len(end_seq))))
            
            continue
            

        # 7.
        # Search for the most outer CS from MCS that are not in 4. 
        rsite1, ename1, rsite2, ename2 = find_compatible_MCS_rsites(MCS, rsitelist, enamelist, full_sequences, backbone_no_MCS_5, backbone_no_MCS_3)
        if not rsite1 or not rsite2:
            print("""Enzyme {} cannot satisfy the conditions.
            No good cutsites in the MCS.""".format(ename0, rsite1, rsite2))
            print("\n\n##################################################################\n\n")
            continue

        # 8. 
        full_plasmid = assemble_plasmid(backbone_no_MCS_5, backbone_no_MCS_3, rsite1 + full_sequences[0] + rsite2)

        if full_plasmid.count(rsite0) > 1:
            print("""Enzyme {} cannot satisfy the conditions.
            It cuts the backbone.""".format(ename0, rsite0))
            print("\n\n##################################################################\n\n")
        else:
            if len(compatible_restriction_sites) == 0:
                optimal_plasmid = full_plasmid
                MCS_rsites = (rsite1, rsite2)
            compatible_restriction_sites.append((ename0, rsite_place, ename1, ename2))
            print("""Enzyme {} can be used for integration into the genome by cutting the {}.
            Enzymes that should be used for cloning inside the plasmid are {}, and {}.""".format(
                ename0, rsite_place, ename1, ename2
            ))
            
            start_name = '5` UTR' if args.modality in [5,0] else 'Gene end'
            end_name = '3` UTR' if args.modality in [3,0] else 'Gene start'
            print("\nPieces that should be used for cloning are:\n{}: {}\n{}: {}".format(
                start_name, start_seq, end_name, end_seq
            ))
            print("\nDesired alpha: {}, Realized alpha: {:.2}".format(args.alpha, real_alpha))
            print("\nTotal length of gene chunks that should be subcloned is", str((len(start_seq)+len(end_seq))))

            #12. Find good additonal cutsites to add on the joints between the gene chunks, linker and FPG
                #Go through the list and check if they are good for use with these sequences
            if(args.modality != 0):
                full_plasmids  = [assemble_plasmid(backbone_no_MCS_5, backbone_no_MCS_3, rsite1 + full_sequence + rsite2) for full_sequence in full_sequences]
                good_pop_enzymes, good_pop_cutsites = find_additional_cutsites(full_plasmids, popular_rsitelist, popular_enamelist)
                if(len(good_pop_enzymes)<3):
                    print("There's no enough popular cutsites that can be added to the existing insert")
                else:  
                    #print(good_pop_enzymes)
                    print("Added the cutsites of 1. {first}, 2. {second} and 3. {third} to the final sequence\n".format(first = good_pop_enzymes[0], second = good_pop_enzymes[1], third = good_pop_enzymes[2]))
                    if(args.modality == 5):
                        print("The final insert sequence with the first FPG: {}".format(rsite1 + start_seq + good_pop_cutsites[0] + FPGs[0][:-3] + good_pop_cutsites[1] + linker + good_pop_cutsites[2] + end_seq + rsite2)) #we don't take the stop of the FPG in case of 5' labeling
                    if(args.modality == 3):
                        print("The final insert sequence with the first FPG: {}".format(rsite1 + start_seq + good_pop_cutsites[0] + linker + good_pop_cutsites[1] + FPGs[0] + good_pop_cutsites[2] + end_seq + rsite2)) #stop codon is already removed from start_seq in case of 3' labeling
                    print("\nOther popular enzymes that can be used are {}".format(good_pop_enzymes[3:]))
            print("\n\n##################################################################\n\n")
    return optimal_plasmid, compatible_restriction_sites, MCS_rsites


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Smart pop out arguments")
    parser.add_argument("--backbone_path", type=str)#, default='./Test_examples/pETUL_backbone.fsa', help="")
    parser.add_argument("--MCS_start_ind", type=int)#, default=1, help="")
    parser.add_argument("--MCS_end_ind", type=int)#, default=108, help="")
    parser.add_argument("--min_homology", type=int)#, default=70, help="")
    parser.add_argument("--alpha", type=float)#, default=1.4, help="")
    parser.add_argument('--Gene_path', type=str)#, default='./Test_examples/CLB2_3p_labeling/CLB2_pm_1000.fsa', help='')
    parser.add_argument("--linker_path", type=str, default='')#, default='./Test_examples/long_linker.fsa', help="")
    parser.add_argument("--modality", type=int)#, default=3, help="")
    parser.add_argument("--enzyme_path", type=str)#, default='./Chris_code/raw_enzyme_list.txt', help='')
    parser.add_argument("--popular_enzyme_path", type=str)
    parser.add_argument("--FPG_paths", nargs="*", type=str, default=[
        # './Test_examples/CLB2_3p_labeling/mCherry_FPG.fsa',
        # './Test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa',
        # './Test_examples/CLB2_3p_labeling/ymTq2_FPG.fsa'
    ], help="")

    args = parser.parse_args()
#    print(args)

    main(args)
