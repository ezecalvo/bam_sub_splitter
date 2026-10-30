# split_bam_by_conversion

Split an aligned BAM file into **labeled** and **unlabeled** reads based on a specific nucleotide conversion. It's built for metabolic RNA labeling methods that cause base substitutions, such as **SLAM-seq** and **TimeLapse-seq / TT-TimeLapse-seq** (T>C), or other chemistries (for example G>A).

A read is **labeled** if it has at least `--min_conversions` of the chosen conversion compared with the reference genome. For paired-end data, if either mate is labeled, both mates go to the labeled BAM.

## Requirements

- Python ≥ 3.7
- [pysam](https://github.com/pysam-developers/pysam)


## Input

- A **coordinate-sorted, indexed BAM** (`.bam` + `.bai`)
- An **indexed reference genome FASTA** (`.fa` + `.fai`). Use the same one the reads were aligned to.

## Usage

```bash
python split_bam_by_conversion.py \
    -i sample.bam \
    -g genome.fa \
    -o output_dir/ \
    -c T>C \
    -t 8 \
    --min_conversions 2
```

### Arguments

| Argument | Description | Default |
|---|---|---|
| `-i`, `--input_bam` | Input BAM file (sorted and indexed) | required |
| `-o`, `--output_dir` | Output directory (created if missing) | required |
| `-g`, `--genome` | Indexed reference FASTA | required |
| `-c`, `--conversion` | Conversion to detect, as `Ref>Read` (e.g. `T>C`, `G>A`). `TC`, `T2C` and `T-C` also work | `T>C` |
| `-t`, `--threads` | Number of worker processes (work is split by chromosome) | `1` |
| `--min_conversions` | Minimum conversions needed to call a read labeled | `1` |

## Output

For an input `sample.bam`:

| File | Description |
|---|---|
| `sample_labeled.bam` (+ `.bai`) | Reads/pairs with ≥ N conversions, sorted and indexed |
| `sample_unlabeled.bam` (+ `.bai`) | All other reads, sorted and indexed |
| `subs_summary.tsv` | Per-sample counts (each run adds one row) |

`subs_summary.tsv` columns:

```
Sample  Conversion  Labeled  Unlabeled  Prop_Labeled  Prop_Unlabeled
```

Counts are **reads**, not pairs. Several samples written to the same output directory share one summary file, so you can run it in a loop:

```bash
for bam in bams/*.bam; do
    python split_bam_by_conversion.py -i "$bam" -g genome.fa -o split/ -t 8
done
```

## How strand is handled

Conversions are checked against the strand each read aligns to:

- **Forward-strand reads:** looks for the conversion as given (e.g. `T>C`)
- **Reverse-strand reads:** looks for the complement (e.g. `A>G`)

This assumes the **read orientation matches the RNA strand**, as in forward-stranded libraries like QuantSeq FWD, which is common for SLAM-seq. For **reverse-stranded** libraries (e.g. dUTP/TruSeq stranded, where read 1 is antisense), pass the complementary conversion instead (`-c A>G` to detect T>C on the RNA).

## Notes and limitations

- Secondary and supplementary alignments are skipped and don't appear in either output.
- Unmapped reads aren't included in the outputs.
- Mates mapped to different chromosomes are classified separately.
- Base quality isn't considered, and known SNPs aren't masked. Real T>C SNPs in your sample will create false labeling, so consider filtering SNP positions (or raising `--min_conversions`) if this matters for your analysis.
- Only mismatches at aligned positions are counted; insertions, deletions and soft-clipped bases are ignored.
