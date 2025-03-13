import pandas as pd

# Load dataset
file_path = "c:/Users/Cyberplus/Project_SARS/data/merged_sarscov2_metadata.csv"
df = pd.read_csv(file_path)

# Drop columns with too many missing values
df.drop(columns=["Last vaccinated", "Additional host information"], inplace=True)

# Convert collection date to datetime
df["Collection date"] = pd.to_datetime(df["Collection date"], errors="coerce")

# Fill missing values with placeholders (or choose another strategy)
df.fillna({"Additional location information": "Unknown", 
           "Sampling strategy": "Unknown", 
           "Specimen": "Unknown"}, inplace=True)

# Print cleaned dataset summary
print("\n✅ Cleaned Metadata Overview:")
print(df.info())
print(f"\n📊 Cleaned Dataset Shape: {df.shape[0]} rows × {df.shape[1]} columns")
