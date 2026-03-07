import re
import argparse
from Bio import SeqIO

def extract_last_transmembrane(tmhmm_file):
    TMhelix_dict = {}
    current_id = None
    transmem_list = []

    with open(tmhmm_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith("ID"):
                current_id=line.split()[1]
            elif line.startswith("FT") and "TRANSMEM" in line:
                match = re.search(r'TRANSMEM\s+(\d+)\s+(\d+)',line)
                start,end = int(match.group(1)),int(match.group(2))
                TMhelix_dict.setdefault(current_id, []).append((start, end))
     
    return TMhelix_dict

def check_motifs_in_inside(fasta_file, tmhmm_file):
    TMhelix_regions = extract_last_transmembrane(tmhmm_file)
    motif_pattern = re.compile(r'Y\w\w[IL]\w{6,12}Y\w\w[IL]') 

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id
        sequence = str(record.seq)
        sequence_length = len(sequence) 

        matches = [(m.start() + 1, m.group()) for m in motif_pattern.finditer(sequence)] 
        if not matches:
            continue

        print(f"\nSequence: {seq_id} | Sequence Length: {sequence_length}") 
        regions = TMhelix_regions.get(seq_id, [])  

        for pos, motif in matches:
            # Skipping if the current sequence has no transmembrane regions.
            if not regions:
                print(f"  Motif at position {pos}: {motif} --> No TMhelix region available <U+274C>")
                continue

            # Retrieving the transmembrane region containing the ITAM motif.
            relevant_region = regions[-1]  
            start_of_TMhelix, end_of_TMhelix = relevant_region 

            # Extracting the sequence of the transmembrane region.
            tmhelix_sequence = sequence[start_of_TMhelix - 1:end_of_TMhelix] 

            # Checking whether the start position of the ITAM motif is greater than the end position of the transmembrane
            if pos > end_of_TMhelix:
                print(f"  Motif at position {pos}: {motif} --> After TMhelix <U+2705>")
            else:
                print(f"  Motif at position {pos}: {motif} --> Before TMhelix <U+274C>")

            # Outputing the start position, end position, and sequence of the transmembrane region.
 print(f"    TMhelix region: Start={start_of_TMhelix}, End={end_of_TMhelix}")
            print(f"    TMhelix region sequence: {tmhelix_sequence}")

            # Outputing the start position and sequence of the motif.
            motif_sequence = sequence[pos - 4:pos + len(motif) - 1] 
            print(f"    complete Motif sequence: {motif_sequence}")
def main():
    parser = argparse.ArgumentParser(description="Check motif positions and whether they fall after TMHMM-predicted 'TMhelix' regions.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-t", "--tmhmm", required=True, help="TMHMM result file")
    args = parser.parse_args()

    check_motifs_in_inside(args.fasta, args.tmhmm)

if __name__ == "__main__":
    main()