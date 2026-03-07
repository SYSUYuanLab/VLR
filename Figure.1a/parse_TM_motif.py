import argparse
import csv
import re

def parse_tmhmm_file(input_file):
    with open(input_file, "r", encoding="utf-8") as f:
        content = f.read()

    entries = content.strip().split("Sequence: ")[1:]
    results = []

    for entry in entries:
        lines = entry.strip().splitlines()
        header = lines[0]
        gene_match = re.match(r"(.+?)\s*\|\s*Sequence Length: (\d+)", header)
        if not gene_match:
            continue

        gene_id, length = gene_match.group(1), gene_match.group(2)
        motif_blocks = []
        tm_regions = []

        current_block = []
        for line in lines[1:]:
            if line.strip().startswith("Motif at position"):
                if current_block:
                    motif_blocks.append(current_block)
                current_block = [line]
            elif line.strip().startswith("TMhelix region") or line.strip().startswith("complete Motif sequence"):
                current_block.append(line)
        if current_block:
            motif_blocks.append(current_block)

        has_after = any("After TMhelix" in block[0] for block in motif_blocks)
        has_before = any("Before TMhelix" in block[0] for block in motif_blocks)

        if not has_after:
            continue  # skip entire entry if only Before exists

        # extract TMhelix regions (keep last one only)
        tm_positions = []
        for block in motif_blocks:
            for line in block:
                if line.strip().startswith("TMhelix region:"):
                    match = re.search(r"Start=(\d+), End=(\d+)", line)
                    if match:
                        tm_positions.append((int(match.group(1)), int(match.group(2))))
        last_tm_start, last_tm_end = tm_positions[-1] if tm_positions else ("", "")

        for block in motif_blocks:
            location_line = block[0]
            if "After TMhelix" in location_line or "Before TMhelix" in location_line:
                motif_region_match = re.search(r"Motif at position (\d+):", location_line)
                motif_region = motif_region_match.group(1) if motif_region_match else ""

                complete_seq = ""
                tm_seq = ""
                for line in block:
                    if line.strip().startswith("complete Motif sequence:"):
                        complete_seq = line.strip().split(":")[1].strip()
                    elif line.strip().startswith("TMhelix region sequence:"):
                        tm_seq = line.strip().split(":")[1].strip()

                results.append([
                    gene_id,
                    length,
                    complete_seq,
                    motif_region,
                    tm_seq,
                    f'"{last_tm_start}-{last_tm_end}"' if last_tm_start != "" else ""
                ])

    return results

def write_csv(results, output_file):
    headers = ["Proteins", "Length(AA)", "Motif", "Motif_region", "TM", "TM_region"]
    with open(output_file, "w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        writer.writerows(results)

def main():
    parser = argparse.ArgumentParser(description="Parse TMHMM-style motif prediction results.")
    parser.add_argument("-i", "--input", required=True, help="Input text file")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    args = parser.parse_args()

    results = parse_tmhmm_file(args.input)
    write_csv(results, args.output)

if __name__ == "__main__":
    main()