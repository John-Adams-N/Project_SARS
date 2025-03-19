import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Get script directory & define data path dynamically
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, "..", "data")
file_path = os.path.join(data_dir, "merged_sarscov2_metadata.csv")

# Load cleaned metadata
df = pd.read_csv(file_path)

# 🔹 Explicitly ensure "Collection date" is datetime
df["Collection date"] = pd.to_datetime(df["Collection date"], errors="coerce")

# 📊 Variant Distribution
plt.figure(figsize=(8, 5))
sns.countplot(data=df, x="Variant", hue="Variant", legend=False)  # Fixes warning
plt.title("Variant Distribution")
plt.xlabel("SARS-CoV-2 Variant")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.show()

# 📅 Collection Date Trends (Grouped by Month)
df["Collection Period"] = df["Collection date"].dt.to_period("M")  # Now this works!

plt.figure(figsize=(12, 6))
sns.countplot(data=df, x="Collection Period", order=sorted(df["Collection Period"].astype(str).unique()))
plt.title("Distribution of Collection Dates (Grouped by Month)")
plt.xlabel("Collection Period")
plt.ylabel("Frequency")
plt.xticks(rotation=45)
plt.show()

# 🧬 Mutation Patterns
plt.figure(figsize=(12, 5))
df["AA Substitutions"].str.split(",").explode().value_counts().head(20).plot(kind="bar", color="coral")
plt.title("Top 20 Most Frequent Mutations")
plt.xlabel("Mutation")
plt.ylabel("Count")
plt.xticks(rotation=90)
plt.show()