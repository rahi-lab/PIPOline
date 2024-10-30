# PIPOline

PIPOline is a tool designed to generate plasmids for pop-out-based tagging or deletion of genes in *Saccharomyces cerevisiae* (budding yeast). It automates the design of DNA constructs for gene tagging at the 5' or 3' ends or for gene deletion, facilitating genetic manipulation experiments.

## Table of Contents

- [Introduction](#introduction)
- [Input Parameters](#input-parameters)
- [Output](#output)
- [Example Usage](#example-usage)
- [PIPOline Workflow](#pipoline-workflow)
- [Notes](#notes)

## Introduction

PIPOline automates the process of designing plasmids for gene manipulation in yeast. It identifies suitable restriction sites for integration and pop-out recombination, assembles the necessary DNA sequences, and ensures compatibility with fluorescent protein genes (FPGs) and linkers when tagging genes.

## Input Parameters

When running PIPOline, the following input parameters should be provided:

- `--backbone_path`: Path to the plasmid backbone file (FASTA format).
- `--MCS_start_ind`: Index of the first nucleotide of the backbone present in the multiple cloning site (MCS) (1-based counting).
- `--MCS_end_ind`: Index of the last nucleotide of the backbone present in the MCS (1-based counting).
- `--min_homology`: Minimum number of base pairs around the cut site used for integration.
- `--R`: Ratio between the homology piece used for integration and the homology piece used for pop-out.
- `--Gene_path`: Path to the gene sequence file (FASTA format). The gene sequence should include the open reading frame (ORF) plus 1000 bp upstream and downstream.
- `--linker_path`: Path to the linker sequence used for tagging (FASTA format).
- `--modality`: Operation mode of the program:
  - `0` for gene deletion.
  - `5` for gene tagging at the 5' terminus.
  - `3` for gene tagging at the 3' terminus.
- `--enzyme_path`: Path to the file containing a list of restriction enzymes for cloning and plasmid linearization before yeast transformation.
- `--popular_enzyme_path`: Path to the file containing a list of popular restriction enzymes to be considered for adding around the linker and FPGs.
- `--FPG_paths`: Paths to the fluorescent protein genes (FASTA format). Used only in the case of tagging. Multiple FPGs can be provided by specifying multiple paths. They should be ORFs including both start and stop codons.

**Note:** Plasmid backbone, gene, linker, and FPGs should be supplied as FASTA (`.fsa`) files.

## Output

For each restriction site present in the gene sequence and surrounding ±1000 bp:

- The program prints whether the site can be used for integrating the final plasmid into the yeast genome.
- If not usable, it provides the reason, which may include:
  - The site cuts the insert (with any of the FPGs) more than once.
  - The site cuts the backbone.
  - There are no suitable cut sites for cloning within the MCS.
- If usable, it generates the left and right DNA sequences to be used for pop-out-based tagging or deletion.
  - For 5' tagging, these are sequences left of the gene (terminator) and the beginning of the gene.
- In the case of gene tagging, it generates an example of the insert that contains the first loaded FPG and additional cut sites around each part of the sequence.

**Note:** Sequences in the output are sorted according to the length of gene sequences used, from smallest to largest.

## Example Usage

### Gene Deletion

```bash
python.exe ./main.py \
  --backbone_path './resources/test_examples/pETUL_backbone.fsa' \
  --MCS_start_ind 1 \
  --MCS_end_ind 108 \
  --min_homology 70 \
  --R 1.4 \
  --Gene_path './resources/test_examples/GAL4_deletion/GAL4_pm_1000_gene.fsa' \
  --modality 0 \
  --enzyme_path './resources/test_examples/raw_enzyme_list.txt'
```

### Gene Tagging
```bash
python.exe ./main.py \
  --backbone_path ./resources/test_examples/pETUL_backbone.fsa \
  --MCS_start_ind 1 \
  --MCS_end_ind 108 \
  --min_homology 70 \
  --R 1.4 \
  --Gene_path ./resources/test_examples/CLB2_3p_labeling/CLB2_pm_1000.fsa \
  --linker_path ./resources/test_examples/long_linker.fsa \
  --modality 5 \
  --enzyme_path ./resources/test_examples/raw_enzyme_list.txt \
  --popular_enzyme_path ./resources/test_examples/Popular_enzymes.txt \
  --FPG_paths ./resources/test_examples/CLB2_3p_labeling/mCherry_FPG.fsa \
               ./resources/test_examples/CLB2_3p_labeling/ymNeonGreen_FPG.fsa \
               ./resources/test_examples/CLB2_3p_labeling/ymTq2_FPG.fsa
```

## PIPOline Workflow

Below is a brief description of the steps performed by PIPOline:

1. **Import Clean Restriction Sites**:  
   Loads a database of clean restriction sites (longer than 5 bp, palindromic, and containing only A/T/G/C).

2. **Load Inputs**:  
   Imports the plasmid backbone, gene, FPGs, MCS positions, minimum homology length, and ratio `R`.

3. **Search for Integration Cut Sites**:  
   Starting from the desired end of the gene, excluding the first `min_homology` base pairs, searches for cut sites that can be used for integration into the yeast genome (`rsite0`).

4. **Assemble Insert Sequence**:  
   Assembles the full sequence for the insert.  
   - For 3' tagging and cutting inside the end of the gene, merges the pop-in sequence (excluding the stop codon) with the linker, FPG, and pop-out sequence.

5. **Check Uniqueness of Cut Site**:  
   Verifies that the cut site from step 3 is unique in the sequence from step 4.

6. **Iterate if Necessary**:  
   If the uniqueness check fails, repeats steps 3–5 starting from the last cut site.

7. **Find Suitable MCS Cut Sites**:  
   Searches for the outermost cut sites in the MCS that are not present in the assembled insert sequence.  
   - Performs two independent searches, one from the left and one from the right.

8. **Assemble Final Plasmid**:  
   Constructs the final plasmid by replacing the MCS with the insert sequence from step 4.

9. **Verify Cut Site in Final Plasmid**:  
   Checks that the cut site from step 3 remains unique in the final plasmid.

10. **Iterate if Necessary**:  
    If the check in step 9 fails, repeats steps 3–9 starting from the last cut site.

11. **Output Results**:  
    - If a valid cut site is found, outputs the pop-in sequence, pop-out sequence, assembled plasmid, and relevant cut sites.
    - If not valid, provides the reason.

12. **Add Additional Cut Sites (Tagging Only)**:  
    For gene tagging, identifies popular restriction enzymes that do not cut the assembled plasmid and adds them between sequence parts.  
    - Outputs the complete sequence and additional enzymes that may be used.

## Notes

- **Sequence Files**: All sequences (plasmid backbone, gene, linker, FPGs) must be provided in FASTA (`.fsa`) format.
- **Modality Parameter**:
  - `0`: Gene deletion.
  - `5`: Gene tagging at the 5' end.
  - `3`: Gene tagging at the 3' end.
- **Multiple FPGs**: When tagging, multiple FPGs can be specified to generate constructs with different fluorescent tags.
- **Output Sorting**: The output sequences are sorted by the length of the gene sequences used, from shortest to longest.
