import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2

# Get the current working directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define file paths dynamically
data_dir = os.path.join(script_dir, "..", "data")

# CSV file for merged SARS-CoV-2 metadata
metadata_file = os.path.join(data_dir, "merged_sarscov2_metadata.csv")

# Ensure the file exists before loading
if not os.path.exists(metadata_file):
    raise FileNotFoundError(f"❌ Metadata file not found: {metadata_file}")

# Load the metadata
df = pd.read_csv(metadata_file, parse_dates=["Collection date"])

# Print a summary to confirm successful loading
print(f"✅ Loaded metadata file: {metadata_file}")
print(df.head())

# Filter Alpha and Omicron variants
alpha_df = df[df["Variant"] == "Alpha"]
omicron_df = df[df["Variant"] == "Omicron"]

# Extract and count mutations for each variant
alpha_mutations = alpha_df["AA Substitutions"].str.split(",").explode().value_counts()
omicron_mutations = omicron_df["AA Substitutions"].str.split(",").explode().value_counts()

# Identify shared and unique mutations
alpha_set = set(alpha_mutations.index)
omicron_set = set(omicron_mutations.index)
shared_mutations = alpha_set & omicron_set
alpha_unique = alpha_set - omicron_set
omicron_unique = omicron_set - alpha_set

# 🏆 Top 20 Mutations in Alpha and Omicron
plt.figure(figsize=(12,5))
sns.barplot(x=alpha_mutations.head(20).index, y=alpha_mutations.head(20).values, color="blue", label="Alpha")
sns.barplot(x=omicron_mutations.head(20).index, y=omicron_mutations.head(20).values, color="red", label="Omicron")
plt.xticks(rotation=90)
plt.title("Top 20 Frequent Mutations in Alpha vs. Omicron")
plt.xlabel("Mutation")
plt.ylabel("Count")
plt.legend()
plt.show()

# 🔵 Venn Diagram: Unique vs. Shared Mutations
plt.figure(figsize=(6,6))
venn2([alpha_set, omicron_set], set_labels=("Alpha", "Omicron"))
plt.title("Mutation Overlap Between Alpha and Omicron")
plt.show()

# Save results
alpha_mutations.to_csv(os.path.join(data_dir, "alpha_mutation_counts.csv"))
omicron_mutations.to_csv(os.path.join(data_dir, "omicron_mutation_counts.csv"))
