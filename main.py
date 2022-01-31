import argparse
import numpy as np
import re

def read_from_fsa(fsa_file_path):
    if not fsa_file_path:
        return '',''
    try:
        with open(fsa_file_path) as f:
            name = f.readline()[1:-1]
            code = f.read().replace('\n', '')
            return name, code.upper()
    except Exception as e:
        print('Unable to read:', fsa_file_path,e)

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
    return (np.array(rsitelist), np.array(enamelist))

# Read in the gene file +/- 1000 base pairs and break it up into 3 pieces
def read_gene_plus_string(Gene_plus):

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


def find_rsite_locations(Gene, rsitelist, enamelist, starting_end, min_homology=0):

    if starting_end == 5:
        Gene_cut = Gene[min_homology:]
        rsite_position_list = np.array([Gene_cut.find(rsite) for rsite in rsitelist])
        sorted_inds = np.argsort(rsite_position_list)
        
    elif starting_end == 3:
        Gene_cut = Gene[:len(Gene)-min_homology]
        rsite_position_list = np.array([Gene_cut.rfind(rsite) for rsite in rsitelist]) 
        sorted_inds = np.argsort(-rsite_position_list)


    # remove not matched rsites
    rsite_position_list = rsite_position_list[sorted_inds]
    matched_inds = (rsite_position_list != -1)
    rsite_position_list = rsite_position_list[matched_inds]
    if starting_end == 5:
        rsite_position_list += min_homology

    # apply sort and remove inds to rsitelist and enamelist
    rsitelist = rsitelist[sorted_inds][matched_inds]
    enamelist = enamelist[sorted_inds][matched_inds]

    return rsitelist, enamelist, rsite_position_list


def generate_start_end_sequences(left_chunk, right_chunk, rsite_side, rsite_position, rsite, minhomology, alpha, stop_codon_offset=0):

    if rsite_side not in ['left','right']:
        raise ValueError('Wrong Gene_side value')

    if rsite_side == 'left':
        X = left_chunk[rsite_position-minhomology:len(left_chunk)-stop_codon_offset]
        alphaX = right_chunk[:int(alpha*len(X))]
        return X, alphaX
    else:
        X = right_chunk[:rsite_position+len(rsite)+minhomology]
        alphaX = left_chunk[-int(alpha*len(X)):len(left_chunk)-stop_codon_offset]
        return alphaX, X


def rsite_search(Gene, rsitelist, enamelist, modality, alpha, min_homology=0,
                 left_of_Gene='', right_of_Gene='', FPGs=[], linker=''):

    if modality not in [0, 3, 5]:
        raise ValueError('Wrong starting_position value')

    if modality == 5:

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = find_rsite_locations(
            Gene, rsitelist, enamelist, 5, min_homology)
        rsitelist_5UTR, enamelist_5UTR, rsite_position_list_5TR = find_rsite_locations(
            left_of_Gene, rsitelist, enamelist, 3, min_homology)

        start_end_sequences_gene = [generate_start_end_sequences(
            left_of_Gene, Gene, 'right', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)]
        start_end_sequences_5UTR = [generate_start_end_sequences(
            left_of_Gene, Gene, 'left', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_5TR, rsitelist_5UTR)]

        full_sequences = [
            [start_seq+FPG+linker+end_seq for FPG in FPGs]
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_5UTR
        ]

        return_dict = {
            'rsitelist': np.append(rsitelist_gene, rsitelist_5UTR),
            'enamelist': np.append(enamelist_gene, enamelist_5UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_5TR),
            'rsite_places': ['rsite in gene']*len(rsitelist_gene)+['rsite in 5`UTR']*len(rsitelist_5UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_5UTR 
        }

    elif modality == 3:

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = find_rsite_locations(
            Gene, rsitelist, enamelist, 3, min_homology)
        rsitelist_3UTR, enamelist_3UTR, rsite_position_list_3UTR = find_rsite_locations(
            right_of_Gene, rsitelist, enamelist, 5, min_homology)

        start_end_sequences_gene = [
            generate_start_end_sequences(Gene, right_of_Gene, 'left', rsite_pos, rsite, min_homology, alpha, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)
        ]
        start_end_sequences_3UTR = [
            generate_start_end_sequences(Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, alpha, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_3UTR, rsitelist_3UTR)
        ]

        full_sequences = [
            [start_seq+linker+FPG+end_seq for FPG in FPGs]
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_3UTR
        ]

        return_dict = {
            'rsitelist': np.append(rsitelist_gene, rsitelist_3UTR),
            'enamelist': np.append(enamelist_gene, enamelist_3UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_3UTR),
            'rsite_places': ['rsite in gene']*len(rsitelist_gene)+['rsite in 3`UTR']*len(rsitelist_3UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_3UTR
        }

    else:
        rsitelist_5, enamelist_5, rsite_position_list_5 = find_rsite_locations(left_of_Gene, rsitelist,
                                                                            enamelist, 3, min_homology)
        rsitelist_3, enamelist_3, rsite_position_list_3 = find_rsite_locations(right_of_Gene, rsitelist,
                                                                            enamelist, 5, min_homology)

        start_end_sequences_5UTR = [generate_start_end_sequences(
            left_of_Gene, Gene, 'left', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_5, rsitelist_5)]
        start_end_sequences_3UTR = [generate_start_end_sequences(
            Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, alpha) for rsite_pos, rsite in zip(rsite_position_list_3, rsitelist_3)]

        full_sequences = [
            [start_seq+end_seq]
            for start_seq, end_seq in start_end_sequences_5UTR+start_end_sequences_3UTR
        ]

        return_dict = {
            'rsitelist': np.append(rsitelist_5, rsitelist_3),
            'enamelist': np.append(enamelist_5, enamelist_3),
            'rsite_position_list': np.append(rsite_position_list_5, rsite_position_list_3),
            'rsite_places': ['rsite in 5`']*len(rsitelist_5)+['rsite in 3`']*len(rsitelist_3),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_5UTR+start_end_sequences_3UTR
        }

    full_sequences_per_rsite = return_dict['full_sequences']
    sort = np.argsort([len(x[0]) for x in full_sequences_per_rsite])
    for key, value in return_dict.items():
        return_dict[key] = list(np.array(value)[sort])

    return return_dict


def find_compatible_MCS_rsites(MCS, rsitelist, enamelist, full_sequences, backbone_no_MCS_5, backbone_no_MCS_3):

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


def assemble_plasmid(backbone_no_MCS_5, backbone_no_MCS_3, sequence): 
    # TODO CHECK if backbone_no_MCS_5 + sequence + backbone_no_MCS_3 + backbone_no_MCS_5 bad
    return backbone_no_MCS_5 + sequence + backbone_no_MCS_3


def main(args):
    # Read all the files (1. and 2.)
    _, backbone = read_from_fsa(args.backbone_path)
    # print('Backbone len:', len(backbone))

    MCS = backbone[args.MCS_start_ind-1:args.MCS_end_ind]
    backbone_no_MCS_5 = backbone[0:args.MCS_start_ind-1]
    backbone_no_MCS_3 = backbone[args.MCS_end_ind:]
    # print('MCS len', len(MCS))

    _, linker = read_from_fsa(args.linker_path)
    # print('Linker len:', len(linker))
    _, Gene_plus = read_from_fsa(args.Gene_path)
    # print('Gene plus len:', len(Gene_plus))

    left_of_Gene, Gene, right_of_Gene = read_gene_plus_string(Gene_plus)
    # print('Gene len:', len(Gene))

    FPGs = [read_from_fsa(path)[1] for path in args.FPG_paths] if len(args.FPG_paths) else []
    # print('FPG len:', len(FPGs))

    rsitelist, enamelist = read_enzyme_list(args.enzyme_path)
    print(rsitelist[:10], enamelist[:10])


    # 3.
    rsite_dict = rsite_search(Gene, rsitelist, enamelist, args.modality, args.alpha,
                              args.min_homology, left_of_Gene, right_of_Gene, FPGs, linker)

    gene_rsitelist_sorted = rsite_dict['rsitelist']
    gene_enamelist_sorted = rsite_dict['enamelist'] 
    gene_rsite_position_list_sorted = rsite_dict['rsite_position_list'] 
    rsite_places = rsite_dict['rsite_places'] 
    full_sequences_per_rsite = rsite_dict['full_sequences']
    start_end_sequences = rsite_dict['start_end_sequences']
    # print(gene_rsitelist_sorted, gene_enamelist_sorted, gene_rsite_position_list_sorted, rsite_places)

    compatible_restriction_sites = []
    optimal_plasmid = ''
    MCS_rsites = []

    # print([len(x[0]) for x in full_sequences_per_rsite])

    for i in range(len(gene_rsitelist_sorted)):

        print('\n')

        # 4. and 5.
        rsite0 = gene_rsitelist_sorted[i]
        ename0 = gene_enamelist_sorted[i]
        rsite_place = rsite_places[i]
        full_sequences = full_sequences_per_rsite[i]
        start_seq, end_seq = start_end_sequences[i]
        print(i, rsite0, ename0, 'side:',rsite_place)
            
        # check if all full sequences have only one rsite0  
        if any([(full_sequence.count(rsite0) != 1) for full_sequence in full_sequences]):
            print("""Enzyme {} cannot satisfy the conditions.
            rsite {} is not unique in the full sequence""".format(ename0, rsite0))
            continue
            

        # 7.
        # Search for the most outer CS from MCS that are not in 4. 
        rsite1, ename1, rsite2, ename2 = find_compatible_MCS_rsites(MCS, rsitelist, enamelist, full_sequences, backbone_no_MCS_5, backbone_no_MCS_3)
        if not rsite1 or not rsite2:
            print("""Enzyme {} cannot satisfy the conditions.
            good rsites not found in MCS""".format(ename0, rsite1, rsite2))
            continue

        # 8. 
        full_plasmid = assemble_plasmid(backbone_no_MCS_5, backbone_no_MCS_3, full_sequences[0])

        if full_plasmid.count(rsite0) != 1:
            print("""Enzyme {} cannot satisfy the conditions.
            bad backbone for rsite0""".format(ename0, rsite0))
        else:
            if len(compatible_restriction_sites) == 0:
                optimal_plasmid = full_plasmid
                MCS_rsites = (rsite1, rsite2)
            compatible_restriction_sites.append((ename0, rsite_place, ename1, ename2))
            print("""Enzyme {} can be used for cuttin the {}.
            Enzymes that should be used for cloning inside the plasmid are {}, and {}.""".format(
                ename0, rsite_place, ename1, ename2
            ))
            start_name = '5UTR' if args.modality in [5,0] else 'Gene end'
            end_name = '3UTR' if args.modality in [3,0] else 'Gene start'
            print("Pieces that should be used for cloning are -- {}:{}, {}:{}".format(
                start_name, start_seq, end_name, end_seq
            ))

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
    parser.add_argument("--FPG_paths", nargs="*", type=str, default=[
        # './Test_examples/CLB2_3p_labeling/mCherry_FPG.fsa',
        # './Test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa',
        # './Test_examples/CLB2_3p_labeling/ymTq2_FPG.fsa'
    ], help="")

    args = parser.parse_args()
    print(args)

    main(args)

    # python.exe ./main.py --backbone_path './Test_examples/pETUL_backbone.fsa' --MCS_start_ind 1 --MCS_end_ind 108 --min_homology 70 --alpha 1.4 --Gene_path './Test_examples/GAL4_deletion/GAL4_pm_1000_gene.fsa' --modality 0 --enzyme_path './Chris_code/raw_enzyme_list.txt'

    # python.exe ./main.py --backbone_path './Test_examples/pETUL_backbone.fsa' --MCS_start_ind 1 --MCS_end_ind 108 --min_homology 70 --alpha 1.4 --Gene_path './Test_examples/CLB2_3p_labeling/CLB2_pm_1000.fsa' --linker_path './Test_examples/long_linker.fsa' --modality 5 --enzyme_path './Chris_code/raw_enzyme_list.txt' --FPG_paths './Test_examples/CLB2_3p_labeling/mCherry_FPG.fsa' './Test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa' './Test_examples/CLB2_3p_labeling/ymTq2_FPG.fsa'