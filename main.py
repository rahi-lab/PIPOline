import argparse
import numpy as np
from utils import *


class RSiteInfo:
    """ Contains all the information that describes one restriction site
    """

    def __init__(
        self, Gene, modality, R, min_homology, left_of_Gene, right_of_Gene, FPGs, FPGs_names, rsite0,
        ename0, rsite_place, full_sequences, gene_rsite_position, start_seq, end_seq, real_R
    ):
        self.Gene = Gene
        self.modality = modality
        self.R = R
        self.min_homology = min_homology
        self.left_of_Gene = left_of_Gene
        self.right_of_Gene = right_of_Gene
        self.FPGs = FPGs
        self.FPGs_names = FPGs_names
        self.rsite0 = rsite0
        self.ename0 = ename0
        self.rsite_place = rsite_place
        self.full_sequences = full_sequences
        self.gene_rsite_position = gene_rsite_position
        self.start_seq = start_seq
        self.end_seq = end_seq
        self.real_R = real_R


class PIPOline:
    """ Main class of the project, initialized with set of constraints (backbone, MCS, linker and enzymes).
    It finds and validates restriction sites for given Gene and modality (3'/5' tagging or deletion) 
    """

    def __init__(self, backbone_path, MCS_start_ind, MCS_end_ind, linker_path, enzyme_path):
        # Read the files (1. and 2.)
        print(
            '\n\n\nPIPOline by Stojkovic, Gligorovski, Rahi\n\n\n'
            '***************************\n'
            '********** Read input files\n\nVector backbone:'
        )
        _, self.backbone = read_from_fsa(backbone_path)
        
        print('\nMultiple cloning site extracted from vector backbone:')
        self.MCS = self.backbone[MCS_start_ind-1:MCS_end_ind]
        print(self.MCS,'\nLength:', len(self.MCS))
        self.backbone_no_MCS_5 = self.backbone[0:MCS_start_ind-1]
        self.backbone_no_MCS_3 = self.backbone[MCS_end_ind:]

        print('\nLinker:')
        _, self.linker = read_from_fsa(linker_path)

        self.rsitelist, self.enamelist = read_enzyme_list(enzyme_path)


    def find_rsite_locations(self, sequence, rsitelist, enamelist, starting_end, min_homology=0):
        """Searches for cutsites from "rsitelist" in a given sequence. 
        Starts by ommiting "min_homology" bp from the "starting_end" (either 5 or 3 as in 5' and 3', respectively)"""
        
        sequence_cut = sequence[min_homology:len(sequence)-min_homology]
        
        if starting_end == 5:
            rsite_position_list = np.array([sequence_cut.find(rsite) for rsite in rsitelist])
            sorted_inds = np.argsort(rsite_position_list)

        elif starting_end == 3:
            rsite_position_list = np.array([sequence_cut.rfind(rsite) for rsite in rsitelist]) 
            sorted_inds = np.argsort(-rsite_position_list)

        # remove not matched rsites
        rsite_position_list = rsite_position_list[sorted_inds]
        matched_inds = (rsite_position_list != -1)
        rsite_position_list = rsite_position_list[matched_inds]
        rsite_position_list += min_homology

        # apply sort and remove inds to rsitelist and enamelist
        rsitelist = rsitelist[sorted_inds][matched_inds]
        enamelist = enamelist[sorted_inds][matched_inds]

        return rsitelist, enamelist, rsite_position_list


    def generate_start_end_sequences(self, left_chunk, right_chunk, rsite_side, rsite_position, rsite, minhomology, R, stop_codon_offset=0):
        """Generates pop_in and pop_out sequences for the cases of tagging genes (at either 5' or 3' end)

        Splitted in two cases, depending on the position of the cutsite:
        rsite_side argument can have 'left' or 'right' values. 
        For 5' UTR tagging, 'left' means inside the gene starting from 5' end, 'right' means inside the right 1000 bp
        For 3' UTR tagging, 'left' means inside the left 1000 bp, 'right' means inside the gene starting from 3' end
        For gene deletion, 'left' means inside the left 1000 bp, 'right' means inside the right 1000 bp
        stop_codon_offset is used to remove the stop codon of the gene in case we are tagging the 3' end """
        

        if rsite_side not in ['left','right']:
            raise ValueError('Wrong riste_side value')

        if rsite_side == 'left':
            pop_in = left_chunk[rsite_position-minhomology:len(left_chunk)-stop_codon_offset]
            pop_out = right_chunk[:int(R*len(pop_in))]
            return pop_in, pop_out
        else:
            pop_in = right_chunk[:rsite_position+minhomology]
            pop_out = left_chunk[-int(R*len(pop_in)):len(left_chunk)-stop_codon_offset]
            return pop_out, pop_in


    def rsite_search_5_tag(self, Gene, R, min_homology, left_of_Gene, right_of_Gene, FPGs):
        """Finds rsites for 5` tag modality
        """

        Gene_and_right_of_Gene = Gene + right_of_Gene

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = self.find_rsite_locations(
            Gene_and_right_of_Gene, self.rsitelist, self.enamelist, 5, min_homology)
        rsitelist_5UTR, enamelist_5UTR, rsite_position_list_5TR = self.find_rsite_locations(
            left_of_Gene, self.rsitelist, self.enamelist, 3, min_homology)

        start_end_sequences_gene = [self.generate_start_end_sequences(
            left_of_Gene, Gene_and_right_of_Gene, 'right', rsite_pos, rsite, min_homology, R) for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)]
        start_end_sequences_5UTR = [self.generate_start_end_sequences(
            left_of_Gene, Gene_and_right_of_Gene, 'left', rsite_pos, rsite, min_homology, R) for rsite_pos, rsite in zip(rsite_position_list_5TR, rsitelist_5UTR)]

        full_sequences = [
            [start_seq+FPG[:-3]+self.linker+end_seq for FPG in FPGs] #We don't take the stop codon from FPG
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_5UTR
        ]
        real_Rs = [len(pop_out)/float(len(pop_in)) for pop_out,pop_in in start_end_sequences_gene] + \
            [len(pop_out)/float(len(pop_in)) for pop_in,pop_out in start_end_sequences_5UTR]

        return {
            'rsitelist': np.append(rsitelist_gene, rsitelist_5UTR),
            'enamelist': np.append(enamelist_gene, enamelist_5UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_5TR),
            'rsite_places': ['in gene+3`UTR']*len(rsitelist_gene)+['in 5`UTR']*len(rsitelist_5UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_5UTR,
            'real_Rs': real_Rs
        }


    def rsite_search_3_tag(self, Gene, R, min_homology, left_of_Gene, right_of_Gene, FPGs):
        """Finds rsites for 3` tag modality
        """

        left_of_Gene_and_Gene = left_of_Gene + Gene

        rsitelist_gene, enamelist_gene, rsite_position_list_gene = self.find_rsite_locations(
            left_of_Gene_and_Gene, self.rsitelist, self.enamelist, 3, min_homology)
        rsitelist_3UTR, enamelist_3UTR, rsite_position_list_3UTR = self.find_rsite_locations(
            right_of_Gene, self.rsitelist, self.enamelist, 5, min_homology)

        start_end_sequences_gene = [
            self.generate_start_end_sequences(left_of_Gene_and_Gene, right_of_Gene, 'left', rsite_pos, rsite, min_homology, R, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_gene, rsitelist_gene)
        ]
        start_end_sequences_3UTR = [
            self.generate_start_end_sequences(left_of_Gene_and_Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, R, stop_codon_offset=3)
            for rsite_pos, rsite in zip(rsite_position_list_3UTR, rsitelist_3UTR)
        ]

        full_sequences = [
            [start_seq+self.linker+FPG+end_seq for FPG in FPGs]
            for start_seq, end_seq in start_end_sequences_gene+start_end_sequences_3UTR
        ]
        real_Rs = [len(pop_out)/float(len(pop_in)) for pop_in,pop_out in start_end_sequences_gene] + \
            [len(pop_out)/float(len(pop_in)) for pop_out,pop_in in start_end_sequences_3UTR]

        return {
            'rsitelist': np.append(rsitelist_gene, rsitelist_3UTR),
            'enamelist': np.append(enamelist_gene, enamelist_3UTR),
            'rsite_position_list': np.append(rsite_position_list_gene, rsite_position_list_3UTR),
            'rsite_places': ['in 5`UTR+gene']*len(rsitelist_gene)+['in 3`UTR']*len(rsitelist_3UTR),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_gene+start_end_sequences_3UTR,
            'real_Rs': real_Rs
        }


    def rsite_search_delete(self, R, min_homology, left_of_Gene, right_of_Gene):
        """Finds rsites for delete modality
        """

        rsitelist_5, enamelist_5, rsite_position_list_5 = self.find_rsite_locations(left_of_Gene, self.rsitelist,
                                                                                self.enamelist, 3, min_homology)
        rsitelist_3, enamelist_3, rsite_position_list_3 = self.find_rsite_locations(right_of_Gene, self.rsitelist,
                                                                            self.enamelist, 5, min_homology)

        start_end_sequences_5UTR = [self.generate_start_end_sequences(
            left_of_Gene, right_of_Gene, 'left', rsite_pos, rsite, min_homology, R) for rsite_pos, rsite in zip(rsite_position_list_5, rsitelist_5)]
        start_end_sequences_3UTR = [self.generate_start_end_sequences(
            left_of_Gene, right_of_Gene, 'right', rsite_pos, rsite, min_homology, R) for rsite_pos, rsite in zip(rsite_position_list_3, rsitelist_3)]

        full_sequences = [
            [start_seq+end_seq]
            for start_seq, end_seq in start_end_sequences_5UTR+start_end_sequences_3UTR
        ]
        real_Rs = [len(pop_out)/len(pop_in) for pop_in,pop_out in start_end_sequences_5UTR] + \
            [len(pop_out)/len(pop_in) for pop_out,pop_in in start_end_sequences_3UTR]

        return {
            'rsitelist': np.append(rsitelist_5, rsitelist_3),
            'enamelist': np.append(enamelist_5, enamelist_3),
            'rsite_position_list': np.append(rsite_position_list_5, rsite_position_list_3),
            'rsite_places': ['in 5`']*len(rsitelist_5)+['in 3`']*len(rsitelist_3),
            'full_sequences': full_sequences,
            'start_end_sequences': start_end_sequences_5UTR+start_end_sequences_3UTR,
            'real_Rs': real_Rs
        }


    def rsite_search(self, Gene, modality, R, min_homology=0, left_of_Gene='', right_of_Gene='', FPGs=[], FPGs_names = []):
        """This is the central function that finds suitable cutsites for the given list of FPGs and the linker


        Modality is used to define the application of the program:
        5 for 5' tagging with given list of FPGs,
        3 for 3' taggging with given list of FPGs, and
        0 for deletion (FPGs neglected)

        Function returns a dictionary of parameters that describe the cloning """

        if modality not in [0, 3, 5]:
            raise ValueError('Wrong modality value')

        #5' tagging
        if modality == 5:
            return_dict = self.rsite_search_5_tag(Gene, R, min_homology, left_of_Gene, right_of_Gene, FPGs)

        # 3' tagging 
        elif modality == 3:
            return_dict = self.rsite_search_3_tag(Gene, R, min_homology, left_of_Gene, right_of_Gene, FPGs)

        # Deletion
        else:
            return_dict = self.rsite_search_delete(R, min_homology, left_of_Gene, right_of_Gene)

        full_sequences_per_rsite = return_dict['full_sequences']
        sort = np.argsort([len(x[0]) for x in full_sequences_per_rsite])
        for key, value in return_dict.items():
            return_dict[key] = list(np.array(value)[sort])

        # return return_dict
        return [
            RSiteInfo(
                Gene,
                modality,
                R,
                min_homology,
                left_of_Gene,
                right_of_Gene,
                FPGs,
                FPGs_names,
                return_dict['rsitelist'][i],
                return_dict['enamelist'][i],
                return_dict['rsite_places'][i],
                return_dict['full_sequences'][i],
                return_dict['rsite_position_list'][i],
                return_dict['start_end_sequences'][i][0],
                return_dict['start_end_sequences'][i][1],
                return_dict['real_Rs'][i]
            ) for i in range(len(full_sequences_per_rsite))
        ]

    def find_compatible_MCS_rsites(self, full_sequences):
        """Searches for the cutsites that can be used for opening up the bakcbone and putting in the insert. 
        \n Performs two searches: one starting from the 5' end and one starting from the 3' end. 
        We stop as soon we find the first good cutsite from each side. 
        The search strategy is thus optimal in a sense that it will remove as much cutsites from the MCS as it can.
        If rsite0 cuts the remaining MCS sequence it will also cut any other sequence found by further search. 
        For this reason we don't check for rsite0 presence within this function, just find the most outer cutsites in 
        the MCS that uniquely cut the backbone and don't cut the insert sequence
        """

        rsite1 = ''
        ename1 = ''
        rsite2 = ''
        ename2 = ''


        rsitelist_sorted_right, enamelist_sorted_right, _ = self.find_rsite_locations(self.MCS, self.rsitelist, self.enamelist, starting_end=5)
        for i, rsite in enumerate(rsitelist_sorted_right):
            if not any(rsite in full_sequence for full_sequence in full_sequences) and not rsite in (self.backbone_no_MCS_5 + self.MCS[:self.MCS.find(rsite)] + self.backbone_no_MCS_3):
                rsite1 = rsite
                ename1 = enamelist_sorted_right[i]
                break

        rsitelist_sorted_left, enamelist_sorted_left, _ = self.find_rsite_locations(self.MCS, self.rsitelist, self.enamelist, starting_end=3)
        for i, rsite in enumerate(rsitelist_sorted_left):
            if not any(rsite in full_sequence for full_sequence in full_sequences) and not rsite in (self.backbone_no_MCS_5 + self.MCS[:self.MCS.find(rsite1)] + self.MCS[self.MCS.find(rsite)+ len(rsite):] + self.backbone_no_MCS_3):
                rsite2 = rsite
                ename2 = enamelist_sorted_left[i]
                break
            
        return rsite1, ename1, rsite2, ename2


    def find_additional_cutsites(self, plasmids, rsitelist, enamelist):
        """Finds all appropropriate cutsites that could be added between different pieces of the insert in the plasmid 
        
        (e.g. between the gene and the linker or between the linker and the FPG etc).
        \n A good cutsite has to fullfill three criteria:
        1. Divisible by 3 (this could be circuimvented by accomidating the linker length)
        2. Does not cut the final plasmids (plural in case of several FPGs)
        3. Does not introduce a stop codon"
        4. When introduced, does not introduce a rsite0, rsite1 nor rsite2 around the joints of the new sequences"""

        good_rsite_list = []
        good_enzyme_list = []
        for rsite, enzyme in zip(rsitelist, enamelist):
            if len(rsite)%3 == 0 and dna_to_protein(rsite).count('*') == 0:
                if not any(rsite in plasmid for plasmid in plasmids):
                    good_rsite_list.append(rsite)
                    good_enzyme_list.append(enzyme)

        return good_enzyme_list, good_rsite_list
   
    def check_uniqueness_preservation(self, start_seqs, end_seqs, good_enzyme_list, good_rsite_list, rsite0, rsite1, rsite2): 
        #Checks whether the enzymes from the list can be introduced between start_seq and end_seq without introducing riste0, rsite1 nor riste2
        #if yes, it returns the chosen cutsite (first found in the list read from the top) and updates the lists by removing it
        #If not, it returns None, which means that the list of 'good' cutsites has to be expanded
        
        index = -1 #index of the found cutsite

        #either start_seq or end_seq can be a list in case they're FPGs
        for i, rsite in enumerate(good_rsite_list):
            satisfies_all = 1 #parameter that is 1 in case the ongoing cutsite satisfies all combinations of start_seq and end_seq; 0 otherwise
            for start_seq in start_seqs:
                for end_seq in end_seqs:
                    assembled_seq = start_seq + rsite + end_seq
                    if(rsite0 in assembled_seq or rsite1 in assembled_seq or rsite2 in assembled_seq):
                        satisfies_all = 0
          
            if(satisfies_all == 1):
                        index = i
                        good_site = rsite
                        break #breaks the biggest loop

        if(index == -1): #not found
            return None, good_enzyme_list, good_rsite_list
        else:
            #remove them by .pop() so that they're not used repeatedly for the next steps

            good_enzyme_list.pop(index)
            good_rsite_list.pop(index)
            return good_site, good_enzyme_list, good_rsite_list
    
    def verify_additional_cutsites(self, FPGs, linker, start_seq, end_seq, modality, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2):

        if(modality == 5):
            #the cutsites are added 1. between the start seq and the FPG 2. between the FPG and the linker 3. between the linker and the end_seq 

            good_site1, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(start_seq, FPGs, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
            good_site2, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(FPGs, linker, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
            good_site3, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(linker, end_seq, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
            
            if(any(good_site == None for good_site in [good_site1, good_site2, good_site3])):
                print("\nThere are not enough popular restriction enzymes in your list to assemble the insert.")
            else:
                print(
                    "\nAdded the cut sites 1. {first}, 2. {second} and 3. {third} between the gene-of-interest sequences,"
                    " the linker, and the fluorescent protein to create the final insert sequence".format(
                        first = good_site1, second=good_site2, third=good_site3)
                )
        
        if(modality == 3):
            #the cutsites are added 1. between the start seq and the linker 2. between the linker and the FPGs 3. between the FPGs and the end_seq 

            good_site1, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(start_seq, linker, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
            good_site2, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(linker, FPGs, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
            good_site3, good_enzyme_list, good_rsite_list = self.check_uniqueness_preservation(FPGs, end_seq, good_enzyme_list, good_rsite_list, rsite0, rsite1, riste2)
        


            if(any(good_site == None for good_site in [good_site1, good_site2, good_site3])):
                print("\nThere are not enough popular restriction enzymes in your list to assemble the insert.")
            else:
                print(
                    "\nAdded the cut sites 1. {first}, 2. {second} and 3. {third} between the gene-of-interest sequences,"
                    " the linker, and the fluorescent protein to create the final insert sequence".format(
                        first = good_site1, second=good_site2, third=good_site3)
                )

        #TODO organize this return of None to be recobnizible in the main if needed
        return good_site1, good_site2, good_site3



    def assemble_plasmid(self, backbone_no_MCS_5, backbone_no_MCS_3, sequence): 
        return backbone_no_MCS_5 + sequence + backbone_no_MCS_3

    def print_cutsite_location(self, cutsite, ename, left_part, right_part, FPGs_seq, FPGs_name, linker, modality):
        """Prints out the location of the cutsite that is present in the insert. 
        \n Goes first through all pieces of the insert and if it's not present in any concludes that it's in the joints"""

        in_isolated_parts = 0

        left_part_name = '5` UTR' if modality in [5,0] else 'Gene end'
        right_part_name = '3` UTR' if modality in [3,0] else 'Gene start'
        if(left_part.find(cutsite) > 0):
            print("Cutsite of {}, {} is present in the {}".format(ename, cutsite, left_part_name))
            in_isolated_parts = 1
        if(right_part.find(cutsite) > 0):
            print("Cutsite of {}, {} is present in the {}".format(ename, cutsite, right_part_name))
            in_isolated_parts = 1
        if(linker.find(cutsite) > 0):
            print("Cutsite of {}, {} is present in the linker".format(ename, cutsite, linker))
            in_isolated_parts = 1
        for FPG_name, FPG in zip(FPGs_name, FPGs_seq):
            if(FPG.find(cutsite) > 0):
                print("Cutsite of {}, {} is present in {}".format(ename, cutsite, FPG_name))
                in_isolated_parts = 1
        if(not in_isolated_parts):
            print("Cutsite of {}, {} does not cut any gene parts, linker or any of the FPGs but cuts the final sequence after joining them. Try introducing a buffer sequence between them".format(ename, cutsite))

    def validate_restriction_site(self, rsite_id, rsite_info):
        """ Validates given restriction site (RSiteInfo object)
        Validation checks following things:
        1. Uniqness of rsite in full_sequence
        2. If there are compatible MCS cutsites
        3. If rsite cuts final plasmid
        """

        print("""\n\n\n\n****** Testing restriction site {}, sequence: {}, enzyme: {}, location: {} {}""".format(str(
            rsite_id+1), rsite_info.rsite0, rsite_info.ename0, rsite_info.gene_rsite_position, rsite_info.rsite_place))
            
        # check if all full sequences have only one rsite0  
        if any([(full_sequence.count(rsite_info.rsite0) != 1) for full_sequence in rsite_info.full_sequences]):
            
            print(
                "\nEnzyme {} cannot be used for linearizing and integrating the"
                " plasmid because {} is not unique in the insert sequence.".format(
                    rsite_info.ename0, rsite_info.rsite0)
            )
            self.print_cutsite_location(rsite_info.rsite0, rsite_info.ename0, rsite_info.start_seq, rsite_info.end_seq, rsite_info.FPGs, rsite_info.FPGs_names, self.linker, rsite_info.modality)
            return (None, None, None, None, None)
            
        # Search for the most outer CS from MCS that are not in 4. 
        rsite1, ename1, rsite2, ename2 = self.find_compatible_MCS_rsites(rsite_info.full_sequences)
        if not rsite1 or not rsite2:
            print(
                "\n There are no good cutsites in the MCS to clone the insert"
                " which would use {} cutsite {}.".format(
                    rsite_info.ename0, rsite_info.rsite0)
            )
            return (None, None, None, None, None)

        #find the 5' part of the MCS that is left after cutting the backbone with rsite1, cutsite rsite1 itself not included
        cut_MCS_5 = self.MCS[:self.MCS.find(rsite1)]
        #find the 3' part of the MCS that is left after cutting the backbone with rsite2, cutsite rsite2 itself not included
        cut_MCS_3 = self.MCS[self.MCS.find(rsite2)+len(rsite2):]
        
        full_plasmid = self.assemble_plasmid(
            self.backbone_no_MCS_5 + cut_MCS_5, cut_MCS_3 + self.backbone_no_MCS_3,
            rsite1 + rsite_info.full_sequences[0] + rsite2
        )
        if full_plasmid.count(rsite_info.rsite0) > 1:
            print(
                "\nInsert cannot be used because {} (restriction site {}) can cut"
                " the backbone after insert integration.".format(
                    rsite_info.ename0, rsite_info.rsite0)
            )
            return (None, None, None, None, None)
        else:
            compatible_rsites = (rsite_info.ename0, rsite_info.rsite_place, ename1, ename2)
            print(
                "\nRestricion site {} {} cut by enzyme {} can be used for linearizing"
                " and integrating the plasmid.\n\nFor cloning the insert into the backbone,"
                " restriction sites {} (enzyme {}) and restriction site {} (enzyme {}) can be used.".format(
                    rsite_info.rsite0, rsite_info.rsite_place, rsite_info.ename0, rsite1, ename1, rsite2, ename2)
            )
            
            start_name = '5` UTR' if rsite_info.modality in [5,0] else 'Gene end'
            end_name = '3` UTR' if rsite_info.modality in [3,0] else 'Gene start'
            print(
                "\nInsert will be constructed with the following sequences from"
                " the gene-of-interest:\n{}:\n{}\n{}:\n{}".format(
                    start_name, rsite_info.start_seq, end_name, rsite_info.end_seq)
            )
            print(
                "\nDesired R: {}, Realized R: {:.2} (which can be less"
                " than the desired R if the gene-of-interest sequence is"
                " not long enough to find homology for the popout)".format(
                    rsite_info.R, rsite_info.real_R)
            )
            print(
                "\nTotal length of sequences from the gene-of-interest going into the insert is:",
                str((len(rsite_info.start_seq)+len(rsite_info.end_seq)))
            )

            return (rsite1, rsite2, cut_MCS_5, cut_MCS_3, compatible_rsites)


    def find_optimal_plasmid(self, modality, rsite_info, cut_MCS_5, cut_MCS_3, rsite1, rsite2, popular_enzyme_path, assembled_plasmid_name, Gene_path, FPG_names, FPG_seq, save_optimal_plasmid=False):
       
        """Find good additonal cutsites to add on the joints between the gene chunks, linker and FPG
        Go through the list and check if they are good for use with these sequences"""
        optimal_plasmid = None

        #to construct the optimal plasmid, around the FPG side of the two cutsites that surround it, add placeholder sequence to ensure high digestion efficency
        placeholder_code = 'GGTGGC' #two glycine codons (but will be cut out when cloning in FPG so in principle doesn't matter)
        full_plasmids = [
            self.assemble_plasmid(
                self.backbone_no_MCS_5 + cut_MCS_5, 
                cut_MCS_3 + self.backbone_no_MCS_3, 
                rsite1 + full_sequence + rsite2
            ) for full_sequence in rsite_info.full_sequences
        ]
        
        if(modality == 0):
            optimal_plasmid = self.assemble_plasmid(self.backbone_no_MCS_5 + cut_MCS_5, cut_MCS_3 + self.backbone_no_MCS_3, rsite1 + rsite_info.start_seq + rsite_info.end_seq + rsite2)
            
            if save_optimal_plasmid and assembled_plasmid_name != None:
                with open(assembled_plasmid_name, 'w') as f:
                    f.write(
                        '>Plasmid for deletion of {}. Can be integrated into'
                        ' the budding yeast genome using {}, {} cutsite.'.format(
                            Gene_path.split('/')[-1][0:4], rsite_info.ename0, rsite_info.rsite0)
                    )
                    f.write('\n')
                    f.write(optimal_plasmid)
                f.close()
                # print that we are saved the plasmid to the file
                print('\nOptimal plasmid found and saved to the file {}!'.format(assembled_plasmid_name))
        
        if(modality == 3 or modality == 5):

            popular_rsitelist, popular_enamelist = read_enzyme_list(popular_enzyme_path)
            good_pop_enzymes, good_pop_cutsites = self.find_additional_cutsites(full_plasmids, popular_rsitelist, popular_enamelist)
            
            #three additional cutsites that we introduce into the plasmid, to make the parts more modular
            rsite3, rsite4, rsite5 = self.verify_additional_cutsites(FPG_seq, self.linker, rsite_info.start_seq, rsite_info.end_seq, modality, good_pop_enzymes, good_pop_cutsites, rsite_info.rsite0, rsite1, rsite2)
            
            if any(rsite == None for rsite in [rsite3, rsite4, rsite5]):
                print("\nThere are not enough popular restriction enzymes in your list to assemble the insert.")
                return None
                
        if(modality == 5):
            optimal_plasmid = self.assemble_plasmid(
                self.backbone_no_MCS_5 + cut_MCS_5, 
                cut_MCS_3 + self.backbone_no_MCS_3,
                rsite1 + rsite_info.start_seq + rsite3 + 
                placeholder_code + rsite4 + self.linker + 
                rsite5 + rsite_info.end_seq + rsite2
            )
            print(
                '\nThe insert sequence with a placeholder sequence ({}) in place of the FPG sequence:\n{}'.format(
                    placeholder_code, rsite1 + rsite_info.start_seq + rsite3 +
                    placeholder_code + rsite4 + self.linker +
                    rsite5 + rsite_info.end_seq + rsite2)
            )
            if save_optimal_plasmid and assembled_plasmid_name != None:
                with open(assembled_plasmid_name, 'w') as f:
                    f.write(
                        '>Plasmid for N-terminal tagging of {} with fluorescent proteins. Can be integrated into'
                        ' the budding yeast genome using {}, {} cutsite. Intented to be labeled with {}. '
                        'To ensure high-efficency digestion during FPG cloning, there is'
                        ' a placeholder sequence between the gene piece and the linker'.format(
                            Gene_path.split('/')[-1][0:4], rsite_info.ename0, rsite_info.rsite0, FPG_names)
                    )
                    f.write('\n')
                    f.write(optimal_plasmid)
                f.close()
                # print that we are saved the plasmid to the file
                print('\nOptimal plasmid found and saved to the file {}!'.format(assembled_plasmid_name))

        if(modality == 3):
            optimal_plasmid = self.assemble_plasmid(
                self.backbone_no_MCS_5 + cut_MCS_5, 
                cut_MCS_3 + self.backbone_no_MCS_3, 
                rsite1 + rsite_info.start_seq + rsite3 + 
                self.linker + rsite4 + placeholder_code + 
                rsite5 + rsite_info.end_seq + rsite2
            )
            print(
                '\nThe insert sequence with a placeholder sequence ({}) in place of the FPG sequence:\n{}'.format(
                    placeholder_code, rsite1 + rsite_info.start_seq +
                    rsite3 + self.linker + rsite4 +
                    placeholder_code + rsite5 + rsite_info.end_seq + rsite2)
            )

            if save_optimal_plasmid and assembled_plasmid_name != None:
                with open(assembled_plasmid_name, 'w') as f:
                    f.write(
                        '>Plasmid for C-terminal tagging of {} with fluorescent proteins. Can be integrated into'
                        ' the budding yeast genome using {}, {} cutsite. Intented to be labeled with {}. '
                        'To ensure high-efficency digestion during FPG cloning, there is'
                        ' a placeholder sequence between the gene piece and the linker'.format(
                            Gene_path.split('/')[-1][0:4], rsite_info.ename0, rsite_info.rsite0, FPG_names)
                    )
                    f.write('\n')
                    f.write(optimal_plasmid)
                f.close()
                # print that we are saved the plasmid to the file
                print('\nOptimal plasmid found and saved to the file {}!'.format(assembled_plasmid_name))

            print("\nOther popular enzymes that can be used in place of the three listed above are {}".format(good_pop_enzymes[3:]))
        return optimal_plasmid


    def run(self, Gene_path, modality, R, min_homology, FPG_paths, popular_enzyme_path=None, assembled_plasmid_name=None):
        """Assembles the plasmid by using subfunctions"""
        
        print('\nGene-of-interest sequence (ORF +/- 1000 bps):')
        _, Gene_plus = read_from_fsa(Gene_path)
        left_of_Gene, Gene, right_of_Gene = read_gene_plus_string(Gene_plus)

        print('\nFluorescent protein genes:')
        FPGs = [read_from_fsa(path)[1] for path in FPG_paths] if len(FPG_paths) else []
        FPGs_names = [path.split('/')[-1] for path in FPG_paths]
        rsite_info_list = self.rsite_search(Gene, modality, R, min_homology, left_of_Gene, right_of_Gene, FPGs, FPGs_names)

        compatible_restriction_sites = []
        MCS_rsites = None
        shortest_optimal_plasmid = None
        # shortest_optimal_plasmid_length = float('inf')
        shortest_optimal_plasmid_R = float('inf')

        print(
            "\n\n\n*************************************************************************************************************\n"\
            "********** Loop over all restriction sites that can be used for linearizing the final plasmid for integration"
        )

        # iterate over rsites and test them
        for i, rsite_info in enumerate(rsite_info_list):
            
            # test restriction site, and find compatible MCS rsites
            rsite1, rsite2, cut_MCS_5, cut_MCS_3, compatible_rsites = self.validate_restriction_site(i, rsite_info)

            # restriction site is valid if there are compatible_rsites (ename0, rsite_place, ename1, ename2)
            if compatible_rsites:

                # save compatible rsites and MCS rsites
                compatible_restriction_sites.append(compatible_rsites)
                if len(compatible_restriction_sites) == 1:
                    MCS_rsites = (rsite1, rsite2)

                # find optimal plasmid
                if modality in [3,5,0]:
                    optimal_plasmid = self.find_optimal_plasmid(
                        modality, rsite_info, cut_MCS_5, cut_MCS_3, rsite1, rsite2,
                        popular_enzyme_path, assembled_plasmid_name, Gene_path, FPGs_names, FPGs,
                        save_optimal_plasmid=not shortest_optimal_plasmid
                    )
                    # save shortest optimal plasmid found
                    if optimal_plasmid and not shortest_optimal_plasmid:
                        shortest_optimal_plasmid = optimal_plasmid
                        shortest_optimal_plasmid_R = rsite_info.real_R
        
        return shortest_optimal_plasmid, compatible_restriction_sites, MCS_rsites, shortest_optimal_plasmid_R


def main(args):
    
    pipoline = PIPOline(args.backbone_path, args.MCS_start_ind, args.MCS_end_ind, args.linker_path, args.enzyme_path)

    pipoline.run(
        args.Gene_path, 
        args.modality, 
        args.R, 
        args.min_homology,
        args.FPG_paths, 
        args.popular_enzyme_path, 
        args.assembled_plasmid_name
    )



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PIPOline arguments")
    parser.add_argument("--backbone_path", type=str, default='./Test_examples/pETUL_backbone.fsa', help="Plasmid backbone used for PIPO")
    parser.add_argument("--MCS_start_ind", type=int, default=1, help="Starting nucleotide of multiple cloning sequence (1-based counting)")
    parser.add_argument("--MCS_end_ind", type=int, default=108, help="Ending nucleotide of multiple cloning sequence (1-based counting)")
    parser.add_argument("--min_homology", type=int, default=70, help="Minimal length of DNA sequence used for homologous recombination during pop-in")
    parser.add_argument("--R", type=float, default=2, help="Ratio between DNA sequence length used for pop-out vs. pop-in")
    parser.add_argument('--Gene_path', type=str, default='./Test_examples/CLB2_3p_labeling/CLB2_pm_1000.fsa', help='ORF of the gene with 1000 bp upstream and downstream (in FASTA format)')
    parser.add_argument("--linker_path", type=str, default='./Test_examples/long_linker.fsa', help="Linker between the gene and the FPG")
    parser.add_argument("--modality", type=int, default=0, help="Modality of the program that can assume one of the following values: 0 for gene deletion, 5 for gene tagging at the 5’ terminus or3 for gene tagging at the 3’ terminus")
    parser.add_argument("--enzyme_path", type=str, default='./Test_examples/raw_enzyme_list.txt', help='List of enzymes that can be used for cloning and plasmid linearization before yeast transformation')
    parser.add_argument("--popular_enzyme_path", type=str, help='List of enzymes to be considered for adding around linker and FPG')
    parser.add_argument("--FPG_paths", nargs="*", type=str, default=['./Test_examples/CLB2_3p_labeling/mCherry_FPG.fsa','./Test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa'], help="List of fluorescent protein coding genes that should be considered during plasmid design ")
    parser.add_argument("--assembled_plasmid_name", type = str, default='Assembled_plasmid.fsa', help="Name of the output FASTA file")

    args = parser.parse_args()

    main(args)
