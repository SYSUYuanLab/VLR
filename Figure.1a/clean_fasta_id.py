#start
import argparse
from Bio import SeqIO

def clean_fasta_ids(input_file, output_file):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.description = record.id  # "Removing the description, keeping only the ID.
            SeqIO.write(record, out_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Clean FASTA headers by keeping only the ID.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with cleaned headers")
    args = parser.parse_args()

    clean_fasta_ids(args.input, args.output)

if __name__ == "__main__":
    main()