import argparse
import pandas as pd

def get_validation_args():
    parser = argparse.ArgumentParser(description="Validate microDNA circles by start position.")
    parser.add_argument('--circles',
                        type=str,
                        required=True,
                        help='Path to circles TSV file')
    parser.add_argument('--reads',
                        type=str,
                        required=True,
                        help='Path to soft-clipped reads TSV file')
    parser.add_argument('--start_positions',
                        type=int,
                        nargs='+',
                        required=True,
                        help='One or more start positions to validate (space-separated integers)')
    return parser.parse_args()

def validate_circles_by_start(circles_path, reads_path, start_positions):
    circles = pd.read_csv(circles_path, sep="\t")
    reads_df = pd.read_csv(reads_path, sep="\t")

    validated = []

    for start in start_positions:
        subset = circles[circles['start_pos'] == start]
        if subset.empty:
            print(f"[!] No circle found for start_pos: {start}")
            continue

        circle = subset.iloc[0]
        chrom = circle['chrom']
        end_pos = circle['end_pos']
        start_seq = circle['start_seq']
        end_seq = circle['end_seq']

        start_reads = reads_df[
            (reads_df['chrom'] == chrom) &
            (reads_df['pos'] == start) &
            (reads_df['clip_seq'] == start_seq) &
            (reads_df['clip_type'] == 'start')
        ]

        end_reads = reads_df[
            (reads_df['chrom'] == chrom) &
            (reads_df['pos'] == end_pos) &
            (reads_df['clip_seq'] == end_seq) &
            (reads_df['clip_type'] == 'end')
        ]

        validated.append({
            'start_pos': start,
            'end_pos': end_pos,
            'chrom': chrom,
            'start_seq': start_seq,
            'end_seq': end_seq,
            'start_support_reads': len(start_reads),
            'end_support_reads': len(end_reads),
            'total_support': len(start_reads) + len(end_reads)
        })

    return pd.DataFrame(validated)

def main():
    args = get_validation_args()
    df = validate_circles_by_start(args.circles, args.reads, args.start_positions)
    if df.empty:
        print("No circles validated.")
    else:
        print("\nValidated Circles:")
        print(df.to_string(index=False))
        df.to_csv("validated_circles.tsv", sep='\t', index=False)
        print("\nSaved to validated_circles.tsv")

if __name__ == "__main__":
    main()

