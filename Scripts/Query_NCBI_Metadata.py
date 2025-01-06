from Bio import Entrez

# Configure Entrez
Entrez.email = "johnadams9644@gmail.com"  # Always include your email when using Entrez
Entrez.api_key = "538c00a7151da6a3ff3f8afbc2ff987a4609" # Optional but recommended to avoid rate limits

# Define search query
query = "SARS-CoV-2[Organism] AND S[Gene Name] AND Kenya[Location]"

# Search NCBI for the query
print("Searching NCBI for sequences...")
handle = Entrez.esearch(db="nucleotide", term=query, retmax=100000)  # Adjust retmax to fetch all
record = Entrez.read(handle)
handle.close()

# Get the list of sequence IDs
id_list = record["IdList"]
print(f"Found {len(id_list)} sequences.")

# Fetch the sequences
if id_list:
    print("Fetching sequences...")
    batch_size = 100  # NCBI recommends fetching in batches to avoid timeouts
    output_file = "kenya_spike_sequences.fasta"

    with open(output_file, "w") as file:
        for start in range(0, len(id_list), batch_size):
            end = min(start + batch_size, len(id_list))
            print(f"Fetching sequences {start + 1} to {end}...")
            handle = Entrez.efetch(
                db="nucleotide", id=id_list[start:end], rettype="fasta", retmode="text"
            )
            sequences = handle.read()
            handle.close()
            file.write(sequences)

    print(f"Sequences saved to {output_file}.")
else:
    print("No sequences found for the given query.")
