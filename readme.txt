The purpose of the code is to generate a plasmid for pop-out-based tagging or deletion of genes in budding yeast.

Input:
--backbone_path (path towards the )
--MCS_start_ind (index of the first nucleotide of the backbone present in the MCS, counted from 1)
--MCS_end_ind (index of the first nucleotide of the backbone present in the MCS, counted from 1)
--min_homology (minimal number of basepairs around the cusite used for integration)
--alpha (ratio between the homology piece used for integration and the homology piece used for pop-out)
--Gene_path (path towards the fasta file with gene sequence. gene sequence should be in (ORF +- 1000 bp) format)
--linker_path (path towards the linker sequence used for tagging)
--modality (0 for deletion, 5 or 3 for tagging, depending on the end of the gene that FPGs are being concatenated to)
--enzyme_path (path towards the enzyme list)
--popular_enzymes (path towards the popular enzyme list)
--FPG_paths (paths towards the fluorescent protein genes. Used only in case of tagging. Can be multiple FPGs; if this is the case several pathways should be forwarded. They should be ORFs from start to stop, both included)

Note: plasmid backbone, gene, linker, and FPGs should be supplied as fasta (.fsa) files.

Output: 
For each cutsite present in the gene sequence and surounding 2*1000 bp:
    Print whether it can be used for integrating the final plasmid into the yeast genome.
    If not, say the reason (possible reasons: 1. cuts the insert, with any of the FPGs, more than once, 2. cuts the backbone 3. there are no good cutsites for cloning within the MCS - in certain cases there might be multiple reasons but only one is printed).
    If yes, generate the left and right piece of DNA sequence that should be used for pop-out based tagging/deletion (e.g. in the case of 5' tagging these are sequences left from the gene (terminator) and the beginning of the gene).
    In case of gene tagging: Generate an example of insert that contains the first loaded FPG and additional cutsites around each part of the sequence.

Note: Sequences in the output are sorted according to the length of gene sequences that should be used, from smallest to largest.

Example of a call:
    Gene deletion:
    python.exe ./main.py --backbone_path './Test_examples/pETUL_backbone.fsa' --MCS_start_ind 1 --MCS_end_ind 108 --min_homology 70 --alpha 1.4 --Gene_path './Test_examples/GAL4_deletion/GAL4_pm_1000_gene.fsa' --modality 0 --enzyme_path './Test_examples/raw_enzyme_list.txt'
    
    Gene tagging:
    python.exe ./main.py --backbone_path ./Test_examples/pETUL_backbone.fsa --MCS_start_ind 1 --MCS_end_ind 108 --min_homology 70 --alpha 1.4 --Gene_path ./Test_examples/CLB2_3p_labeling/CLB2_pm_1000.fsa --linker_path ./Test_examples/long_linker.fsa --modality 5 --enzyme_path ./Test_examples/raw_enzyme_list.txt --popular_enzyme_path ./Test_examples/Popular_enzymes.txt --FPG_paths ./Test_examples/CLB2_3p_labeling/mCherry_FPG.fsa ./Test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa ./Test_examples/CLB2_3p_labeling/ymTq2_FPG.fsa

    Note: for Lazar, the '' around sequence location in gene deletion examples are OK (didn't test whether necessary), for Vojislav they're not (the program cannot find the directory).

Pipeline:
    1. Import database of clean RS (clean means longer than 5, palindromes, and containing only A/T/G/C)
    2. Import plasmid, gene, FPGs, MCS position, min_homology, and alpha
    3. Starting from the desired end of the gene and direction, excluding the first min_homology basepairs (bp), search for cutsites that can be used for integration into the yeast genome ("rsite0")
    4. Assemble the full sequence that is going to be used as an insert. e.g. for the 3' tagging and cutting inside the end of the gene, this means merging X (min_homology+len(RS) without STOP codon) with the linker, FPG, and alphaX 
    5. Check for the uniqueness of CS from 3. in sequence from 4.
    6. If not 5. repeat 3. starting from the last CS
    7. Search for the most outer CS from MCS that are not in 4. 
        (two independent searches, one from left, one from right)
    8. Assemble the final plasmid without the MCS part that is replaced with sequence from 4.
    9. Check for the uniqueness of CS from 3. in sequence the final plasmid
    10. If not 9. repeat 3. starting from the last CS
    11. If found a valid CS, return X, alphaX and assembled plasmid, and CS and MCS cut sites
        11.a if the cutsite is present in the sequence but not valid, output the reason 
    12. Only in case of tagging: Find popular restriction enzymes that do not cut the assembled plasmid and add them between the pieces. Output the whole sequence and additional popular enzymes that might be used.
