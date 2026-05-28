import argparse, re, os, glob, csv
import Plot_RAWs

def gc_content(sequence):
    """Calculate the GC content of a DNA sequence."""
    gc_count = sequence.count('G') + sequence.count('C')
    print(sequence)
    print(gc_count)
    return (gc_count / len(sequence)) * 100 if len(sequence) > 0 else 0

def process_fasta(input_file, output_file):
    """Read the .eg file, calculate GC content for each sequence, and write to output file."""
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            parts = line.strip().split('\t')
            if len(parts) != 3:
                continue  # Skip malformed lines
            sequence_index, position_sequence, sequence = parts
            gc = gc_content(sequence)
            kmer = len(sequence)
            outfile.write(f"{sequence_index}\t{position_sequence}\t{sequence}\t{gc:.2f}\t{kmer}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate GC content from .eg file.")
    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input .eg file(s). Globs like k*.eg are allowed (shell-expanded).")
    parser.add_argument(
    "-s", "--save",action="store_true",help="Keep the generated *-gc.eg files (default: delete after merge).")
    parser.add_argument("-p", "--plot", default=None, help="Generate plots. Provide output folder.")
    args = parser.parse_args()
    
    # Expand globs manually (for safety on some systems)
    input_files = []
    for pattern in args.input:
        input_files.extend(glob.glob(pattern))
    if not input_files:
        print("No input files found.")
        exit(1)

    # Process each file
    for input_file in input_files:
        base_name = os.path.basename(input_file)
        output_file = re.sub(r"\.eg$", "-gc.eg", base_name)
        output_file = os.path.join(os.path.dirname(input_file), output_file)

        process_fasta(input_file, output_file)
        print(f"Wrote {output_file}")

    # -----------------------------
    # Merge all -gc.eg into one CSV
    # -----------------------------
    merged_dir = os.path.dirname(os.path.abspath(input_files[0]))
    merged_csv = os.path.join(merged_dir, "merged_gc.csv")

    # Collect the generated -gc.eg files
    generated_gc_files = []
    for input_file in input_files:
        base_name = os.path.basename(input_file)
        gc_name = re.sub(r"\.eg$", "-gc.eg", base_name)
        gc_path = os.path.join(os.path.dirname(input_file), gc_name)
        if os.path.isfile(gc_path):
            generated_gc_files.append(gc_path)
        else:
            print(f"Warning: expected output not found, skipping merge: {gc_path}")

    if not generated_gc_files:
        print("No -gc.eg files found to merge. Skipping CSV merge.")
        return

    # Write merged CSV
    with open(merged_csv, "w", newline="") as csvfile:
        writer = csv.writer(csvfile)
        # writer.writerow(["Sequence_Index", "Position_Sequence", "Sequence", "GC_Content(%)", "Kmer_Length"])

        total_rows = 0
        for gc_file in generated_gc_files:
            with open(gc_file, "r") as fh:
                for line in fh:
                    parts = line.strip().split("\t")
                    if len(parts) != 5:
                        continue
                    sequence_index, position_sequence, sequence, gc, kmer = parts
                    writer.writerow([sequence_index, position_sequence, sequence, gc, kmer])
                    total_rows += 1

    print(f"Merged {len(generated_gc_files)} file(s) → {merged_csv}")

    if not args.save:
        deleted = 0
        for path in generated_gc_files:
            try:
                os.remove(path)
                deleted += 1
            except OSError as e:
                print(f"Warning: couldn't delete {path}: {e}")
    
    if args.plot:
        # Generate plots using Plot_RAWs module
        print("[>] Generating plots...")
        # Infer variant name from folder of the input files
        variant_name = os.path.basename(os.path.dirname(os.path.abspath(input_files[0])))
        print(f"[i] Variant detected: {variant_name}")

        outdir = (args.plot).rstrip("/\\") + "/"
        Plot_RAWs.run_plots(merged_csv, outdir, variant_name)
    
if __name__ == "__main__":
    main()

