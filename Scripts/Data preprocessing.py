from pathlib import Path
import pandas as pd
from Bio import SeqIO

# Define project root directory (parent of "Scripts")
BASE_DIR = Path(__file__).resolve().parent.parent  

# Define the data directory
DATA_DIR = BASE_DIR / "data"

# File paths
fasta_files = {
    "Alpha": DATA_DIR / "Former VOC Alpha.fasta",
    "Omicron": DATA_DIR / "Former VOC Omicron.fasta"
}
metadata_files = {
    "Alpha": DATA_DIR / "Alpha patient metadata.tsv",
    "Omicron": DATA_DIR / "Omicron patient metadata.tsv"
}

# Read metadata files
def load_metadata(file_path):
    if not file_path.exists():
        print(f"❌ ERROR: File not found - {file_path}")
        return None
    return pd.read_csv(file_path, sep="\t")

# Read FASTA sequences
def load_fasta(file_path):
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

# Process variants
merged_data = []
for variant, fasta_file in fasta_files.items():
    metadata_path = metadata_files[variant]

    metadata = load_metadata(metadata_path)
    if metadata is None:
        continue  # Skip if metadata file is missing

    sequences = load_fasta(fasta_file)
    
    # Merge metadata with sequences
    metadata["Sequence"] = metadata["Accession ID"].map(sequences)
    
    # Remove entries with no matching sequence
    metadata = metadata.dropna(subset=["Sequence"])
    metadata["Variant"] = variant  # Add variant column
    merged_data.append(metadata)

# Combine variants
if merged_data:
    final_df = pd.concat(merged_data)
    output_file = DATA_DIR / "merged_sarscov2_data.csv"
    final_df.to_csv(output_file, index=False)
    print(f"✅ Merging complete! Data saved as '{output_file}'")
else:
    print("⚠️ No data was merged. Check file paths.")
