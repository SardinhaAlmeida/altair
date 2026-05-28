import os, sys
import matplotlib.pyplot as plt

def get_variants(base_dir): # Iterate through the variants and their genome counts
    variants_counts = []

    for variant in sorted(os.listdir(base_dir)):
        variant_path = os.path.join(base_dir, variant)
        if not os.path.isdir(variant_path):
            continue 

        genome_count=0
        for file in os.listdir(variant_path):
            if file.lower().endswith((".fasta", ".fa", ".fna")):
                with open(os.path.join(variant_path, file), "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        if line.startswith(">"):
                            genome_count += 1

        variants_counts.append((variant, genome_count))
    
    sorted_variants_counts = sorted(variants_counts, key=lambda x: x[1], reverse=True)
    return sorted_variants_counts

def plot_histogram(variants, counts, virus_name, base_dir):

    plt.figure(figsize=(max(8, 0.8 * len(variants)), 6))
    bars = plt.bar(variants, counts, color="steelblue")
    plt.title(f"{virus_name}")
    plt.xlabel("Variant")
    plt.ylabel("Genome Quantity")
    plt.xticks(rotation=50, fontsize=6)

    # Add numbers above bars
    for bar, val in zip(bars, counts):
        plt.text(bar.get_x() + bar.get_width()/2, val + 0.2, str(val),
                ha='center', va='bottom', fontsize=8)

    plt.tight_layout()

    # --- 4. Save outputs ---
    out_png = os.path.join(base_dir, "Quantity.png")
    # out_pdf = os.path.join(base_dir, "genome_counts.pdf")
    plt.savefig(out_png, dpi=200)
    # plt.savefig(out_pdf)
    plt.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 Histogram_Quantity.py /path/to/pathogen_folder")
        sys.exit(1)

    base_dir = os.path.abspath(sys.argv[1])
    if not os.path.isdir(base_dir):
        print(f"Error: '{base_dir}' is not a valid directory.")
        sys.exit(1)

    virus_name = os.path.basename(base_dir.rstrip('/'))

    variants_counts = get_variants(base_dir)
    variants, counts = zip(*variants_counts) if variants_counts else ([], [])
    plot_histogram(variants, counts, virus_name, base_dir)

if __name__ == "__main__":
    main()
