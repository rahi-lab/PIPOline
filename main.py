import argparse
import os
from pipoline import PIPOline
from utils import search_gene_file, validate_gene_data

RESOURCES_FOLDER = os.path.join(os.path.dirname(__file__), 'resources')
INDIVIDUAL_GENES_FOLDER = os.path.join(RESOURCES_FOLDER, 'individual_genes')
OVERLAPPING_GENES_FILE = os.path.join(RESOURCES_FOLDER, 'overlapping_genes.txt')

def main(args):
    # Check if both args.gene_name and args.Gene_path are provided
    if bool(args.gene_name) == bool(args.gene_path):
        print("Please provide either gene_name or Gene_path, but not both.")
        return

    # Find gene path based on gene name or use the provided gene path
    if args.gene_name:
        gene_path = search_gene_file(INDIVIDUAL_GENES_FOLDER, args.gene_name)
        if gene_path is None:
            print(f"No gene file found matching the name '{args.gene_name}'.")
            return
    else:
        gene_path = args.Gene_path

    # Validate gene data
    validate_gene_data(args.gene_name, gene_path, INDIVIDUAL_GENES_FOLDER, OVERLAPPING_GENES_FILE)

    # Run PIPOline
    pipoline = PIPOline(args.backbone_path, args.MCS_start_ind, args.MCS_end_ind, args.linker_path, args.enzyme_path)
    pipoline.run(
        gene_path, 
        args.modality, 
        args.R, 
        args.min_homology,
        args.FPG_paths, 
        args.popular_enzyme_path, 
        args.assembled_plasmid_name
    )



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="PIPOline arguments")
    parser.add_argument("--backbone_path", type=str, default='./Test_examples/pETUL_backbone.fsa',
                        help="Plasmid backbone used for PIPO")
    parser.add_argument("--MCS_start_ind", type=int, default=1,
                        help="Starting nucleotide of multiple cloning sequence (1-based counting)")
    parser.add_argument("--MCS_end_ind", type=int, default=108,
                        help="Ending nucleotide of multiple cloning sequence (1-based counting)")
    parser.add_argument("--min_homology", type=int, default=70,
                        help="Minimal length of DNA sequence used for homologous recombination during pop-in")
    parser.add_argument("--R", type=float, default=2,
                        help="Ratio between DNA sequence length used for pop-out vs. pop-in")
    parser.add_argument('--Gene_path', type=str, default=None,
                        help='ORF of the gene with 1000 bp upstream and downstream (in FASTA format)')
    parser.add_argument("--gene_name", type=str, default=None,
                        help="Name of the gene that is being cloned")
    parser.add_argument("--linker_path", type=str, default=None,
                        help="Linker between the gene and the FPG")
    parser.add_argument("--modality", type=int, default=0,
                        help="Modality of the program that can assume one of the following values: 0 for gene deletion, 5 for gene tagging at the 5’ terminus or3 for gene tagging at the 3’ terminus")
    parser.add_argument("--enzyme_path", type=str, default=None,
                        help='List of enzymes that can be used for cloning and plasmid linearization before yeast transformation')
    parser.add_argument("--popular_enzyme_path", type=str,
                        help='List of enzymes to be considered for adding around linker and FPG')
    parser.add_argument("--FPG_paths", nargs="*", type=str, default=[],
                        help="List of fluorescent protein coding genes that should be considered during plasmid design ")
    parser.add_argument("--assembled_plasmid_name", type = str, default='Assembled_plasmid.fsa',
                        help="Name of the output FASTA file")

    args = parser.parse_args()

    main(args)
