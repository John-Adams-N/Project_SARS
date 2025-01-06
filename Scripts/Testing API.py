from Bio import Entrez
from pprint import pprint

Entrez.email = "johnadams9644@gmail.com"  # Always include your email when using Entrez
Entrez.api_key = "538c00a7151da6a3ff3f8afbc2ff987a4609" # Optional but recommended to avoid rate limits

def query_ncbi(database, term):
    # Perform the search
    handle = Entrez.esearch(db=database, term=term, retmax=10)  # Adjust retmax as needed
    record = Entrez.read(handle)
    return record

# Example usage
result = query_ncbi("nucleotide", "SARS-CoV-2 spike protein")

pprint(result)
