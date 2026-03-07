import re
import argparse
from Bio import SeqIO

def load_TMhelix_regions(tmhmm_file)
    TMhelix_dict = {}  # Dictionary of transmembrane regions
    with open(tmhmm_file) as f
        for line in f
            parts = line.strip().split()
            if len(parts) = 5 and parts[2] == 'TMhelix'
                seq_id = parts[0]
                start, end = int(parts[3]), int(parts[4])
                TMhelix_dict.setdefault(seq_id, []).append((start, end))  # Storing transmembrane regions in a dictionary.
    return TMhelix_dict

def check_motifs_in_inside(fasta_file, tmhmm_file)
    TMhelix_regions = load_TMhelix_regions(tmhmm_file)
    motif_pattern = re.compile(r'Yww[IL]w{6,12}Yww[IL]')  # Defining a regular expression pattern for ITAM motifs

    for record in SeqIO.parse(fasta_file, fasta)
        seq_id = record.id
        sequence = str(record.seq)
        sequence_length = len(sequence)

        matches = [(m.start() + 1, m.group()) for m in motif_pattern.finditer(sequence)]  # Finding all matches of the ITAM motif.
        if not matches
            continue

        print(fnSequence {seq_id}  Sequence Length {sequence_length})  
        regions = TMhelix_regions.get(seq_id, [])  # Retrieving the transmembrane regions of the current sequence

        for pos, motif in matches
            # Skipping if the current sequence has no transmembrane regions.
            if not regions
                print(f  Motif at position {pos} {motif} -- No TMhelix region available U+274C)
                continue

            # Retrieving the transmembrane region containing the ITAM motif.
            relevant_region = regions[-1]  # Retrieving the first transmembrane region
            start_of_TMhelix, end_of_TMhelix = relevant_region

            # Extracting the sequence of the transmembrane region.
            tmhelix_sequence = sequence[start_of_TMhelix - 1end_of_TMhelix]  

            # Checking whether the start position of the ITAM motif is greater than the end position of the transmembrane region.
            if pos  end_of_TMhelix
                print(f  Motif at position {pos} {motif} -- After TMhelix U+2705)
            else
                print(f  Motif at position {pos} {motif} -- Before TMhelix U+274C)

            # Outputing the start position, end position, and sequence of the transmembrane region.
            print(f    TMhelix region Start={start_of_TMhelix}, End={end_of_TMhelix})
            print(f    TMhelix region sequence {tmhelix_sequence})

            # Outputing the start position and sequence of the motif.
            motif_sequence = sequence[pos - 4pos + len(motif) - 1] 
            print(f    Motif sequence {motif_sequence})
def main()
    parser = argparse.ArgumentParser(description=Check motif positions and whether they fall after TMHMM-predicted 'TMhelix' regions.)
    parser.add_argument(-f, --fasta, required=True, help=Input FASTA file)
    parser.add_argument(-t, --tmhmm, required=True, help=TMHMM result file)
    args = parser.parse_args()

    check_motifs_in_inside(args.fasta, args.tmhmm)

if __name__ == __main__
    main()