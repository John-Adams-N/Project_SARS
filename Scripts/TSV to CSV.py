import pandas as pd
import os

# Define file paths

# Get the current working directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define file paths dynamically
data_dir = os.path.join(script_dir, "..", "data")

# Input files 
data_dir = "data"
tsv_file = os.path.join(data_dir, "metadata.tsv")
csv_file = os.path.join(data_dir, "metadata.csv")

# Ensure data directory exists
os.makedirs(data_dir, exist_ok=True)

try:
    # Read the TSV file
    df = pd.read_csv(tsv_file, sep="\t")

    # Save as CSV
    df.to_csv(csv_file, index=False)

    print(f"✅ Successfully converted '{tsv_file}' to '{csv_file}'.")

except FileNotFoundError:
    print(f"❌ Error: The file '{tsv_file}' was not found. Please check the path.")
except Exception as e:
    print(f"❌ An unexpected error occurred: {e}")
