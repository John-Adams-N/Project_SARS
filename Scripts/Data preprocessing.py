# ==============================================
# SARS-CoV-2 Data Preprocessing (Separated Merging)
# ==============================================
import os
import pandas as pd
from Bio import SeqIO

# Get the current working directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define file paths dynamically
data_dir = os.path.join(script_dir, "..", "data")

fasta_files = {
    "Alpha": os.path.join(data_dir, "Former VOC Alpha.fasta"),
    "Omicron": os.path.join(data_dir, "Former VOC Omicron.fasta")
}
metadata_files = {
    "Alpha": os.path.join(data_dir, "Alpha patient metadata.tsv"),
    "Omicron": os.path.join(data_dir, "Omicron patient metadata.tsv")
}

# Output files
merged_fasta_path = os.path.join(data_dir, "merged_sarscov2_sequences.fasta")
merged_metadata_path = os.path.join(data_dir, "merged_sarscov2_metadata.tsv")

# ==============================================
# Step 1: Merge FASTA Sequences
# ==============================================
def merge_fasta_files(output_path, fasta_dict):
    with open(output_path, "w") as outfile:
        for variant, file_path in fasta_dict.items():
            if not os.path.exists(file_path):
                print(f"❌ ERROR: File not found - {file_path}")
                continue
            
            # Read and write sequences to merged file
            for record in SeqIO.parse(file_path, "fasta"):
                record.description = f"{record.id} | Variant: {variant}"  # Add variant info
                SeqIO.write(record, outfile, "fasta")
    
    print(f"✅ Merged FASTA saved to: {output_path}")

# ==============================================
# Step 2: Merge Metadata
# ==============================================
def merge_metadata_files(output_path, metadata_dict):
    merged_data = []
    
    for variant, file_path in metadata_dict.items():
        if not os.path.exists(file_path):
            print(f"❌ ERROR: File not found - {file_path}")
            continue
        
        # Load metadata
        df = pd.read_csv(file_path, sep="\t")
        df["Variant"] = variant  # Add variant column
        merged_data.append(df)

    if merged_data:
        final_df = pd.concat(merged_data, ignore_index=True)
        final_df.to_csv(output_path, sep="\t", index=False)
        print(f"✅ Merged metadata saved to: {output_path}")
    else:
        print("❌ No metadata files found or loaded.")

# ==============================================
# Run Merging Functions
# ==============================================
merge_fasta_files(merged_fasta_path, fasta_files)
merge_metadata_files(merged_metadata_path, metadata_files)
