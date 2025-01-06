# Fetch metadata for the sequences
print("Fetching metadata...")
handle = Entrez.esummary(db="nucleotide", id=",".join(id_list), retmode="xml")
summary = Entrez.read(handle)
handle.close()

# Save metadata to a CSV file
import csv
output_metadata = "kenya_spike_metadata.csv"
with open(output_metadata, "w", newline="") as csvfile:
    fieldnames = ["Accession", "Title", "Length", "CreateDate"]
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

    for record in summary:
        writer.writerow({
            "Accession": record["AccessionVersion"],
            "Title": record["Title"],
            "Length": record["Length"],
            "CreateDate": record["CreateDate"],
        })

print(f"Metadata saved to {output_metadata}.")
