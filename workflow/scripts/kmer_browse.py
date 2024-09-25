import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import matplotlib.pyplot as plt

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

# Write summarized BED format (consecutive positions with the same k-mer count in one row)
def write_summarized_bed_file(bed_file_handle, record_id, max_counts):
    start = 0
    current_count = max_counts[0]

    for i in range(1, len(max_counts)):
        if max_counts[i] != current_count:
            # Write the previous block
            bed_file_handle.write(f"{record_id}\t{start}\t{i}\t{current_count}\n")
            # Start a new block
            start = i
            current_count = max_counts[i]

    # Write the final block
    bed_file_handle.write(f"{record_id}\t{start}\t{len(max_counts)}\t{current_count}\n")

# Convert counts to FASTQ-compatible quality scores
def convert_counts_to_quality_scores(counts):
    quality_scores = []
    for count in counts:
        # Cap the count between 0 and 93 (since quality scores typically range from 33 to 126)
        score = min(93, max(0, count))
        # Convert count to ASCII (Phred quality score + 33)
        quality_scores.append(chr(score + 33))
    return ''.join(quality_scores)

# Generate and save the histogram of k-mer counts
def generate_histogram(quality_scores, output_path):
    plt.figure()
    plt.hist(quality_scores, bins=20, color='blue', edgecolor='black')
    plt.title('Distribution of Quality Scores')
    plt.xlabel('Phred Quality Score')
    plt.ylabel('Frequency')
    plt.savefig(output_path)
    plt.close()

# Generate and save the line plot of k-mer counts
def generate_line_plot(quality_scores, boundaries, output_path):
    plt.figure()
    plt.plot(quality_scores, color='green')
    plt.title('Quality Scores Across All Positions (with contig separators)')
    plt.xlabel('Position in Sequences')
    plt.ylabel('Phred Quality Score')

    # Add vertical lines to mark the end of each contig
    for boundary in boundaries[:-1]:  # Exclude the final boundary (end of the last contig)
        plt.axvline(x=boundary, color='red', linestyle='--', label='Contig Separator' if boundary == boundaries[0] else "")
    
    plt.savefig(output_path)
    plt.close()

def main():
    # Argument parsing
    parser = argparse.ArgumentParser(description="Generate a FASTQ file with max k-mer counts in the quality field and a summarized BED file with consecutive k-mer count blocks.")
    parser.add_argument("--kmer_count_file", required=True, help="Path to the dumped k-mer count file (e.g., output of 'jellyfish dump')")
    parser.add_argument("--fasta_file", required=True, help="Path to the genome FASTA file")
    parser.add_argument("--kmer_size", type=int, required=True, help="Size of the k-mers")
    parser.add_argument("--output_fastq", required=True, help="Path to the output FASTQ file")
    parser.add_argument("--output_bed", required=True, help="Path to the summarized output BED file")
    parser.add_argument("--output_histogram", required=True, help="Path to save the histogram plot")
    parser.add_argument("--output_line_plot", required=True, help="Path to save the line plot")
    
    args = parser.parse_args()

    # Load the k-mer count table
    kmer_dict = load_kmer_count_table(args.kmer_count_file)

    # Variables to store all quality scores and boundaries for the line plot
    all_quality_scores = []
    contig_boundaries = []

    # Create the output FASTQ and summarized BED files
    with open(args.output_fastq, "w") as fastq_out, open(args.output_bed, "w") as bed_out:
        current_pos = 0
        for record in SeqIO.parse(args.fasta_file, "fasta"):
            sequence = str(record.seq)
            max_kmer_counts = get_kmer_max_counts_in_sequence(sequence, args.kmer_size, kmer_dict)

            # Append quality scores and track the boundary for this contig
            all_quality_scores.extend(max_kmer_counts)
            current_pos += len(max_kmer_counts)
            contig_boundaries.append(current_pos)

            # Write summarized BED file for each contig (FASTA record)
            write_summarized_bed_file(bed_out, record.id, max_kmer_counts)
            
            # Create a FASTQ entry for each contig
            quality_scores = convert_counts_to_quality_scores(max_kmer_counts)
            fastq_entry = f"@{record.id}\n{sequence}\n+\n{''.join(quality_scores)}\n"
            fastq_out.write(fastq_entry)
    
    # Generate and save the plots
    generate_histogram(all_quality_scores, args.output_histogram)
    generate_line_plot(all_quality_scores, contig_boundaries, args.output_line_plot)
    
    print(f"Final FASTQ file saved to {args.output_fastq}")
    print(f"Final summarized BED file saved to {args.output_bed}")
    print(f"Histogram plot saved to {args.output_histogram}")
    print(f"Line plot saved to {args.output_line_plot}")

if __name__ == "__main__":
    main()
