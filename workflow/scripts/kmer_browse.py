import argparse
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Load k-mer count table into a dictionary
def load_kmer_count_table(kmer_count_file):
    kmer_dict = {}
    with open(kmer_count_file, 'r') as f:
        for line in f:
            kmer, count = line.strip().split()
            kmer_dict[kmer] = int(count)
    return kmer_dict

# Process the genome sequence using k-mer counts from the table
def get_kmer_max_counts_in_sequence(sequence, kmer_size, kmer_dict):
    seq_length = len(sequence)
    counts = [0] * seq_length  # Initialize a list to hold the max counts for each nucleotide
    
    # Iterate over the sequence with a sliding window of kmer_size
    for i in range(seq_length - kmer_size + 1):
        kmer = sequence[i:i+kmer_size]
        count = kmer_dict.get(kmer, 0)  # Get count from the table, default to 0 if not found
        
        # Propagate the maximum count to all nucleotides covered by the k-mer
        for j in range(i, i + kmer_size):
            counts[j] = max(counts[j], count)  # Keep the maximum count
    
    return counts

# Convert counts to ASCII-compatible quality scores for FASTQ
def convert_counts_to_quality_scores(counts):
    quality_scores = []
    for count in counts:
        # Cap the count between 0 and 93 (since quality scores typically range from 33 to 126)
        score = min(93, max(0, count))
        # Convert count to ASCII (Phred quality score + 33)
        quality_scores.append(chr(score + 33))
    return ''.join(quality_scores)

def write_bed_file(bed_file_handle, record_id, max_counts):
    """
    Writes a BED file with maximum k-mer counts as score for each nucleotide position.
    """
    for i, count in enumerate(max_counts):
        if count > 0:  # Only include positions with non-zero k-mer counts
            bed_file_handle.write(f"{record_id}\t{i}\t{i+1}\t{count}\n")

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate a FASTQ file with max k-mer counts in the quality field and a BED file with max k-mer count scores.")
    parser.add_argument("--kmer_count_file", required=True, help="Path to the dumped k-mer count file")
    parser.add_argument("--fasta_file", required=True, help="Path to the genome FASTA file")
    parser.add_argument("--kmer_size", type=int, required=True, help="Size of the k-mers")
    parser.add_argument("--output_fastq", required=True, help="Path to the output FASTQ file")
    parser.add_argument("--output_bed", required=True, help="Path to the output BED file")
    
    args = parser.parse_args()

    # Load the k-mer count table
    kmer_dict = load_kmer_count_table(args.kmer_count_file)

    # Create the output FASTQ and BED files
    with open(args.output_fastq, "w") as fastq_out, open(args.output_bed, "w") as bed_out:
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            sequence = str(record.seq)
            max_kmer_counts = get_kmer_max_counts_in_sequence(sequence, args.kmer_size, kmer_dict)
            
            # Write BED file for each contig (FASTA record)
            write_bed_file(bed_out, record.id, max_kmer_counts)
            
            # Create a FASTQ entry for each contig
            quality_scores = convert_counts_to_quality_scores(max_kmer_counts)
            seq_record = SeqRecord(
                Seq(sequence),  # Full sequence
                id=record.id,
                description=record.description,
                letter_annotations={"phred_quality": [ord(q) - 33 for q in quality_scores]}
            )
            SeqIO.write(seq_record, fastq_out, "fastq")
    
    print(f"Final FASTQ file saved to {args.output_fastq}")
    print(f"Final BED file saved to {args.output_bed}")

if __name__ == "__main__":
    main()
