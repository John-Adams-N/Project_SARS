import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Load cleaned metadata
file_path = "c:/Users/Cyberplus/Project_SARS/data/merged_sarscov2_metadata.csv"
df = pd.read_csv(file_path, parse_dates=["Collection date"])

# 📊 Variant Distribution
plt.figure(figsize=(8,5))
sns.countplot(data=df, x="Variant", palette="coolwarm")
plt.title("Variant Distribution")
plt.xlabel("SARS-CoV-2 Variant")
plt.ylabel("Count")
plt.xticks(rotation=45)
plt.show()

# 📅 Collection Date Trends
plt.figure(figsize=(10,5))
df["Collection date"].hist(bins=20, color="skyblue", edgecolor="black")
plt.title("Distribution of Collection Dates")
plt.xlabel("Collection Date")
plt.ylabel("Frequency")
plt.xticks(rotation=45)
plt.show()

# 🧬 Mutation Patterns
plt.figure(figsize=(12,5))
df["AA Substitutions"].str.split(",").explode().value_counts().head(20).plot(kind="bar", color="coral")
plt.title("Top 20 Most Frequent Mutations")
plt.xlabel("Mutation")
plt.ylabel("Count")
plt.xticks(rotation=90)
plt.show()
