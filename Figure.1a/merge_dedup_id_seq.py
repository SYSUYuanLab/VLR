import argparse
from Bio import SeqIO

def merge_and_deduplicate_by_id(input_files, output_file):
    seen_ids = set()
    unique_records = []

    for file in input_files:
        for record in SeqIO.parse(file, "fasta"):
            record_id = record.id
            if record_id not in seen_ids:
                seen_ids.add(record_id)
                unique_records.append(record)

    with open(output_file, "w") as out_handle:
        SeqIO.write(unique_records, out_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Merge FASTA files and remove records with duplicate IDs.")
    parser.add_argument("-i", "--inputs", nargs='+', required=True, help="Input FASTA files")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file with unique IDs")
    args = parser.parse_args()

    merge_and_deduplicate_by_id(args.inputs, args.output)

if __name__ == "__main__":
    main()