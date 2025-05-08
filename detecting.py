import pysam
from collections import defaultdict
import argparse
import pandas as pd

BIN_SIZE = 1000
MAX_DISTANCE = 1000
MIN_HOMOLOGY_SCORE = 0.8

def get_args():
    parser = argparse.ArgumentParser(description="Detect soft-clipped reads as microDNA candidates.")
    parser.add_argument('--bam',
                        type=str, 
                        required=True, 
                        help='Path to input BAM file')
    parser.add_argument('--min_clip', 
                        type=int, 
                        default=10, 
                        help='Minimum length of soft-clipped sequence')
    parser.add_argument('--out', 
                        type=str, 
                        default='soft_clipped_reads.tsv', 
                        help='Output TSV filename')
    return parser.parse_args()

def extract_soft_clipped_reads(bam_path, min_clip):
    soft_clipped_reads = []

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue

            cigar = read.cigartuples
            if not cigar:
                continue

            if cigar[0][0] == 4 and cigar[0][1] >= min_clip:
                clipped_seq = read.query_sequence[:cigar[0][1]]
                position = read.reference_start
                soft_clipped_reads.append({
                    "read_name": read.query_name,
                    "chrom": read.reference_name,
                    "pos": position,
                    "clip_type": "start",
                    "clip_len": cigar[0][1],
                    "clip_seq": clipped_seq,
                    "cigar": read.cigarstring
                })

            elif cigar[-1][0] == 4 and cigar[-1][1] >= min_clip:
                clipped_seq = read.query_sequence[-cigar[-1][1]:]
                position = read.reference_end
                soft_clipped_reads.append({
                    "read_name": read.query_name,
                    "chrom": read.reference_name,
                    "pos": position,
                    "clip_type": "end",
                    "clip_len": cigar[-1][1],
                    "clip_seq": clipped_seq,
                    "cigar": read.cigarstring
                })

    return soft_clipped_reads


def write_to_tsv(reads, out_file):
    with open(out_file, "w") as f:
        f.write("read_name\tchrom\tpos\tclip_type\tclip_len\tclip_seq\tcigar\n")
        for r in reads:
            f.write(f"{r['read_name']}\t{r['chrom']}\t{r['pos']}\t{r['clip_type']}\t{r['clip_len']}\t{r['clip_seq']}\t{r['cigar']}\n")

def sw_fill_matrix(A, B, gap, miss, match):
    H = [[0 for j in range(len(B) + 1)] for i in range(len(A) + 1)]

    for i in range(1, len(A) + 1):
        for j in range(1, len(B) + 1):
            H[i][j] = max( H[i-1][j-1] + (match if A[i-1] == B[j-1] else miss),
                           H[i-1][j] + gap,
                           H[i][j-1] + gap,
                           0)
    return H

def sw_traceback(H, A, B, gap, miss, match):
    i = None
    j = None
    max_score = 0
    for _i in range(len(A) + 1):
        for _j in range(len(B) + 1):
            if H[_i][_j] > max_score:
                max_score = H[_i][_j]
                i = _i
                j = _j

    if max_score == 0 or i is None or j is None:
        return "", "", 0

    score = H[i][j]
    align_A = []
    align_B = []
    while H[i][j] > 0:
        if H[i][j] == H[i-1][j-1] + (match if A[i-1] == B[j-1] else miss):
            align_A.append(A[i-1])
            align_B.append(B[j-1])
            i -= 1
            j -= 1
        elif H[i][j] == H[i-1][j] + gap:
            align_A.append(A[i-1])
            align_B.append('-')
            i -= 1
        elif H[i][j] == H[i][j-1] + gap:
            align_A.append('-')
            align_B.append(B[j-1])
            j -= 1
        else:
            break
    return ''.join(align_A[::-1]), ''.join(align_B[::-1]), score

def sw(seq1, seq2, gap=-2, miss=-1, match=1):
    H = sw_fill_matrix(seq1, seq2, gap, miss, match)
    align1, align2, score = sw_traceback(H, seq1, seq2, gap, miss, match)
    match = ['|' if align1[i] == align2[i] else ' ' for i in range(len(align1))]
    return score, align1, align2, ''.join(match)

def compute_microhomology(seq1, seq2):
    score, aligned_seq1, aligned_seq2, alignment_str = sw(seq1, seq2, gap=-2, miss=-1, match=1)
    max_possible = min(len(seq1), len(seq2))
    return score / max_possible if max_possible > 0 else 0.0

def cluster_and_score(start_df, end_df):
    start_df = start_df.copy()
    end_df = end_df.copy()
    start_df['bin'] = (start_df['pos'] // BIN_SIZE).astype(int)
    end_df['bin'] = (end_df['pos'] // BIN_SIZE).astype(int)

    candidate_pairs = []

    for offset in [-1, 0, 1]:
        temp = end_df.copy()
        temp['bin'] += offset

        merged = pd.merge(
            start_df, temp,
            on=['chrom', 'bin'],
            suffixes=('_start', '_end')
        )

        merged['distance'] = merged['pos_end'] - merged['pos_start']
        merged = merged[(merged['distance'] > 0) & (merged['distance'] <= MAX_DISTANCE)]

        microhomology_scores = [
            compute_microhomology(row['clip_seq_start'], row['clip_seq_end'])
            for _, row in merged.iterrows()
        ]
        merged['microhomology'] = microhomology_scores

        merged = merged[merged['microhomology'] >= MIN_HOMOLOGY_SCORE]

        merged['confidence'] = merged['microhomology'] * (1 - (merged['distance'] / MAX_DISTANCE))

        candidate_pairs.append(merged)

    if not candidate_pairs:
        return pd.DataFrame()  

    all_candidates = pd.concat(candidate_pairs, ignore_index=True)

    circles = all_candidates[[
        'chrom', 'pos_start', 'pos_end', 'clip_seq_start', 'clip_seq_end',
        'distance', 'microhomology', 'confidence'
    ]].rename(columns={
        'pos_start': 'start_pos',
        'pos_end': 'end_pos',
        'clip_seq_start': 'start_seq',
        'clip_seq_end': 'end_seq'
    }).sort_values(by='confidence', ascending=False).reset_index(drop=True)

    return circles

if __name__ == "__main__":
    args = get_args()
    reads = extract_soft_clipped_reads(args.bam, args.min_clip)
    print(f"Found {len(reads)} soft-clipped reads.")

    write_to_tsv(reads, args.out)
    print(f"Wrote output to: {args.out}")

    df = pd.DataFrame(reads)
    start_clips = df[df['clip_type'] == 'start']
    end_clips = df[df['clip_type'] == 'end']

    circles = cluster_and_score(start_clips, end_clips)
    print(f"Identified {len(circles)} candidate microDNA circles.")

    circles.to_csv("microdna_circles.tsv", sep='\t', index=False)
    print("Candidate circles written to microdna_circles.tsv")