import pysam
import argparse
import os
import sys
import multiprocessing
import time

COMPLEMENT_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def parse_args():
    parser = argparse.ArgumentParser(description="Split BAM files based on specific nucleotide conversions (e.g., T>C, G>A).")
    
    parser.add_argument("-i", "--input_bam", required=True, help="Path to input BAM file")
    parser.add_argument("-o", "--output_dir", required=True, help="Directory for output files")
    parser.add_argument("-g", "--genome", required=True, help="Path to indexed reference genome FASTA")
    parser.add_argument("-c", "--conversion", type=str, default="T>C", 
                        help="Conversion to look for (e.g., 'T>C', 'G>A'). Format must be Ref>Read. Default: T>C")
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads to use")
    parser.add_argument("--min_conversions", type=int, default=1, help="Min conversions required to label a read")
    
    return parser.parse_args()

def get_conversion_logic(conv_str):
    clean_str = conv_str.upper().replace(">", "").replace("2", "").replace("-", "")
    
    if len(clean_str) != 2:
        sys.exit(f"Error: Conversion format '{conv_str}' not recognized. Use format like 'T>C'.")
        
    fwd_ref = clean_str[0]
    fwd_alt = clean_str[1]
    
    # check reverse strand
    rev_ref = fwd_ref.translate(COMPLEMENT_MAP)
    rev_alt = fwd_alt.translate(COMPLEMENT_MAP)
    
    return (fwd_ref, fwd_alt), (rev_ref, rev_alt)

def check_conversions(read, fasta_handle, min_conv, fwd_pair, rev_pair):
    try:
        ref_seq = fasta_handle.fetch(read.reference_name, read.reference_start, read.reference_end).upper()
    except (ValueError, IndexError):
        return False

    pairs = read.get_aligned_pairs(matches_only=True)
    conversions = 0
    read_seq = read.query_sequence.upper()
    is_reverse = read.is_reverse
    
    
    if is_reverse:
        target_ref, target_read = rev_pair
    else:
        target_ref, target_read = fwd_pair

    for q_idx, r_idx in pairs:
        rel_r_idx = r_idx - read.reference_start
        if rel_r_idx >= len(ref_seq):
            continue
            
        ref_base = ref_seq[rel_r_idx]
        read_base = read_seq[q_idx]
        
        #check specific conversion
        if ref_base == target_ref and read_base == target_read:
            conversions += 1
            
        if conversions >= min_conv:
            return True
            
    return False

def worker_fn(args_tuple):
    (refs_to_process, input_bam_path, genome_path, temp_prefix, min_conv, fwd_pair, rev_pair) = args_tuple
    pysam.set_verbosity(0)

    bam_in = pysam.AlignmentFile(input_bam_path, "rb")
    fasta = pysam.FastaFile(genome_path)
    
    out_labeled_path = f"{temp_prefix}_labeled.bam"
    out_unlabeled_path = f"{temp_prefix}_unlabeled.bam"
    
    bam_labeled = pysam.AlignmentFile(out_labeled_path, "wb", template=bam_in)
    bam_unlabeled = pysam.AlignmentFile(out_unlabeled_path, "wb", template=bam_in)
    
    local_labeled_count = 0
    local_unlabeled_count = 0
    
    read_buffer = {}

    try:
        for chrom in refs_to_process:
            for read in bam_in.fetch(chrom):
                
                if read.is_secondary or read.is_supplementary:
                    continue

                qname = read.query_name

                if qname in read_buffer:
                    mate = read_buffer.pop(qname)
                    
                    is_read_labeled = check_conversions(read, fasta, min_conv, fwd_pair, rev_pair)
                    is_mate_labeled = check_conversions(mate, fasta, min_conv, fwd_pair, rev_pair)
                    
                    if is_read_labeled or is_mate_labeled:
                        bam_labeled.write(mate)
                        bam_labeled.write(read)
                        local_labeled_count += 2
                    else:
                        bam_unlabeled.write(mate)
                        bam_unlabeled.write(read)
                        local_unlabeled_count += 2
                else:
                    read_buffer[qname] = read

        
        for qname, read in read_buffer.items():
            if check_conversions(read, fasta, min_conv, fwd_pair, rev_pair):
                bam_labeled.write(read)
                local_labeled_count += 1
            else:
                bam_unlabeled.write(read)
                local_unlabeled_count += 1

    finally:
        bam_labeled.close()
        bam_unlabeled.close()
        bam_in.close()
        fasta.close()
        
    return (out_labeled_path, out_unlabeled_path, local_labeled_count, local_unlabeled_count)

def main():
    args = parse_args()
    
    # checks
    if not os.path.exists(args.input_bam):
        sys.exit(f"Error: Input file {args.input_bam} not found.")
    if not os.path.exists(args.genome):
        sys.exit(f"Error: Genome file {args.genome} not found.")
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    fwd_pair, rev_pair = get_conversion_logic(args.conversion)

    with pysam.AlignmentFile(args.input_bam, "rb") as f:
        all_refs = list(f.header.references)

    num_chunks = min(args.threads, len(all_refs))
    chunks = [[] for _ in range(num_chunks)]
    for i, ref in enumerate(all_refs):
        chunks[i % num_chunks].append(ref)

    worker_args = []
    base_name = os.path.splitext(os.path.basename(args.input_bam))[0]
    
    for i in range(num_chunks):
        temp_prefix = os.path.join(args.output_dir, f"temp_worker_{i}_{base_name}")
        worker_args.append((chunks[i], args.input_bam, args.genome, temp_prefix, args.min_conversions, fwd_pair, rev_pair))

    with multiprocessing.Pool(processes=args.threads) as pool:
        results = pool.map(worker_fn, worker_args)

    #aggregate
    temp_labeled_files = []
    temp_unlabeled_files = []
    total_labeled = 0
    total_unlabeled = 0

    for res in results:
        temp_labeled_files.append(res[0])
        temp_unlabeled_files.append(res[1])
        total_labeled += res[2]
        total_unlabeled += res[3]

    
    unsorted_labeled = os.path.join(args.output_dir, f"{base_name}_labeled_unsorted.bam")
    unsorted_unlabeled = os.path.join(args.output_dir, f"{base_name}_unlabeled_unsorted.bam")
    
    final_labeled = os.path.join(args.output_dir, f"{base_name}_labeled.bam")
    final_unlabeled = os.path.join(args.output_dir, f"{base_name}_unlabeled.bam")

    pysam.cat("-o", unsorted_labeled, *temp_labeled_files)
    pysam.cat("-o", unsorted_unlabeled, *temp_unlabeled_files)

    print("Sorting output BAMs...", file=sys.stderr)
    pysam.sort("-o", final_labeled, unsorted_labeled)
    pysam.sort("-o", final_unlabeled, unsorted_unlabeled)

    
    for f in temp_labeled_files + temp_unlabeled_files:
        if os.path.exists(f): os.remove(f)
    if os.path.exists(unsorted_labeled): os.remove(unsorted_labeled)
    if os.path.exists(unsorted_unlabeled): os.remove(unsorted_unlabeled)

    print("Indexing...", file=sys.stderr)
    pysam.index(final_labeled)
    pysam.index(final_unlabeled)
    
    # Stats
    total_reads = total_labeled + total_unlabeled
    if total_reads > 0:
        prop_labeled = total_labeled / total_reads
        prop_unlabeled = total_unlabeled / total_reads
    else:
        prop_labeled = 0
        prop_unlabeled = 0

    summary_path = os.path.join(args.output_dir, "subs_summary.tsv")
    write_header = not os.path.exists(summary_path)
    
    with open(summary_path, "a") as f:
        if write_header:
            f.write("Sample\tConversion\tLabeled\tUnlabeled\tProp_Labeled\tProp_Unlabeled\n")
        f.write(f"{base_name}\t{args.conversion}\t{total_labeled}\t{total_unlabeled}\t{prop_labeled:.5f}\t{prop_unlabeled:.5f}\n")
    
    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    main()
    
