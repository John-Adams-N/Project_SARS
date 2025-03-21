import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import LabelEncoder
from matplotlib_venn import venn2
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import shap
import joblib
from sklearn.ensemble import RandomForestClassifier

# Get the current working directory of the script
script_dir = os.path.dirname(os.path.abspath(__file__))

# Define file paths dynamically
data_dir = os.path.join(script_dir, "..", "data")
model_dir = os.path.join(script_dir, "..", "models")

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

# Load data
data_path = os.path.join(data_dir, "merged_sarscov2_metadata.csv")  # Adjust path as needed
df = pd.read_csv(data_path)

# 🔹 Step 1: Extract Unique Mutations
df["Mutations"] = df["AA Substitutions"].str.split(",")  # Convert mutation lists
all_mutations = set(df["Mutations"].explode())  # Get unique mutations

# 🔹 Step 2: Create a Binary Mutation Matrix
mutation_df = pd.DataFrame(0, index=df.index, columns=list(all_mutations))

for i, mutations in enumerate(df["Mutations"]):
    if isinstance(mutations, list):
        mutation_df.loc[i, mutations] = 1  # Set 1 if mutation is present

# 🔹 Step 3: Encode Variant Labels
label_encoder = LabelEncoder()
df["Variant_Label"] = label_encoder.fit_transform(df["Variant"])

# 🔹 Step 4: Combine Features and Labels
final_df = pd.concat([mutation_df, df["Variant_Label"]], axis=1)

# Save processed dataset
final_df.to_csv(os.path.join(data_dir, "mutation_features.csv"), index=False)

print("Feature Engineering Complete! 🚀")

# Additional code for Logistic Regression Model

# Set paths
feature_file = os.path.join(data_dir, "mutation_features.csv")

# Load the feature dataset
df = pd.read_csv(feature_file)

# Split into features (X) and labels (y)
X = df.drop(columns=["Variant_Label"])
y = df["Variant_Label"]

# Train-test split (80% training, 20% testing)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# Train Logistic Regression Model
model = LogisticRegression(max_iter=500)
model.fit(X_train, y_train)

# Predictions
y_pred = model.predict(X_test)

# Evaluation Metrics
accuracy = accuracy_score(y_test, y_pred)
print(f"✅ Logistic Regression Accuracy: {accuracy:.4f}")
print("\nClassification Report:\n", classification_report(y_test, y_pred))

# Confusion Matrix
plt.figure(figsize=(6,4))
sns.heatmap(confusion_matrix(y_test, y_pred), annot=True, fmt="d", cmap="Blues", xticklabels=["Alpha", "Omicron"], yticklabels=["Alpha", "Omicron"])
plt.xlabel("Predicted")
plt.ylabel("Actual")
plt.title("Confusion Matrix")
plt.show()

# Train and save Random Forest model
model_dir = os.path.join(script_dir, "..", "models")
os.makedirs(model_dir, exist_ok=True)  # Ensure directory exists
model_path = os.path.join(model_dir, "random_forest_model.pkl")

rf_model = RandomForestClassifier(n_estimators=100, random_state=42)
rf_model.fit(X_train, y_train)

joblib.dump(rf_model, model_path)
print(f"✅ Random Forest model saved at: {model_path}")

# Additional code for SHAP analysis

# Load trained Random Forest model
rf_model_path = os.path.join(model_dir, "random_forest_model.pkl")

# Ensure the model file exists before loading
if not os.path.exists(rf_model_path):
    raise FileNotFoundError(f"❌ Random Forest model file not found: {rf_model_path}")

rf_model = joblib.load(rf_model_path)

# 📊 Bar Plot - Feature Importance from Random Forest
feature_importance = rf_model.feature_importances_
sorted_idx = np.argsort(feature_importance)[::-1][:10]  # Top 10 mutations

plt.figure(figsize=(10, 5))
sns.barplot(x=X.columns[sorted_idx], y=feature_importance[sorted_idx], palette="coolwarm")
plt.xticks(rotation=90)
plt.xlabel("Mutation")
plt.ylabel("Feature Importance Score")
plt.title("Top 10 Mutations Contributing to Variant Classification (Random Forest)")
plt.show()

# 📊 SHAP Summary Plot
explainer = shap.TreeExplainer(rf_model)
shap_values = explainer.shap_values(X)

# Ensure the shape of shap_values matches the shape of X
assert shap_values[1].shape == X.shape, "The shape of the shap_values matrix does not match the shape of the provided data matrix."

plt.figure(figsize=(10, 6))
shap.summary_plot(shap_values[1], X, show=False)  # Class 1 (Omicron)
plt.title("SHAP Summary Plot: Mutation Impact on Variant Classification")
plt.show()
