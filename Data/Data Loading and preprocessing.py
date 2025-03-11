# Data Loading and preprocessing

from Bio import SeqIO
import pandas as pd

# File paths (adjust if needed)
fasta_files = {
    "Alpha": "Former VOC Alpha.fasta",
    "Omicron": "Former VOC Omicron.fasta"
}
metadata_files = {
    "Alpha": "Alpha patient metadata.tsv",
    "Omicron": "Omicron patient metadata.tsv"
}

# Read metadata files
def load_metadata(file_path):
    return pd.read_csv(file_path, sep="\t")

# Read FASTA sequences and store in dictionary
def load_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Process each variant
merged_data = []
for variant, fasta_file in fasta_files.items():
    metadata = load_metadata(metadata_files[variant])
    sequences = load_fasta(fasta_file)
    
    # Merge metadata with sequences
    metadata["Sequence"] = metadata["Accession ID"].map(sequences)
    
    # Remove entries with no matching sequence
    metadata = metadata.dropna(subset=["Sequence"])
    metadata["Variant"] = variant  # Add variant column
    merged_data.append(metadata)

# Combine both variants
final_df = pd.concat(merged_data)

# Save to CSV
final_df.to_csv("merged_sarscov2_data.csv", index=False)

print("✅ Merging complete! Data saved as 'merged_sarscov2_data.csv'.")