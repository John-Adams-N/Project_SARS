from Bio import SeqIO
import pandas as pd

# File paths (adjusted to absolute paths)
fasta_files = {
    "Alpha": "C:/Users/Admin.DESKTOP-8TK90VT/Project_SARS/Data/Former VOC Alpha.fasta",
    "Omicron": "C:/Users/Admin.DESKTOP-8TK90VT/Project_SARS/Data/Former VOC Omicron.fasta"
}
metadata_files = {
    "Alpha": "C:/Users/Admin.DESKTOP-8TK90VT/Project_SARS/Data/Alpha patient metadata.tsv",
    "Omicron": "C:/Users/Admin.DESKTOP-8TK90VT/Project_SARS/Data/Omicron patient metadata.tsv"
}

# Read metadata files
def load_metadata(file_path):
    try:
        return pd.read_csv(file_path, sep="\t")
    except FileNotFoundError:
        print(f"\u274C ERROR: File not found - {file_path}")
        return None

# Read FASTA sequences and store in dictionary
def load_fasta(file_path):
    sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq)
    except FileNotFoundError:
        print(f"\u274C ERROR: File not found - {file_path}")
    return sequences

# Process each variant
merged_data = []
for variant, fasta_file in fasta_files.items():
    metadata = load_metadata(metadata_files[variant])
    if metadata is None:
        continue
    sequences = load_fasta(fasta_file)
    
    if not sequences:
        continue
    
    # Merge metadata with sequences
    metadata["Sequence"] = metadata["Accession ID"].map(sequences)
    
    # Remove entries with no matching sequence
    metadata = metadata.dropna(subset=["Sequence"])
    metadata["Variant"] = variant  # Add variant column
    merged_data.append(metadata)

# Combine both variants
if merged_data:
    final_df = pd.concat(merged_data)
    final_df.to_csv("C:/Users/Admin.DESKTOP-8TK90VT/Project_SARS/merged_sarscov2_data.csv", index=False)
    print("\u2705 Merging complete! Data saved as 'merged_sarscov2_data.csv'.")
else:
    print("\u274C No matching data found. Check metadata and FASTA files.")
