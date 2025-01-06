from Bio import Entrez

Entrez.email = "johnadams9644@gmail.com"

handle = Entrez.esearch(db="nucleotide", term="SARS-CoV-2[Organism] AND S[Gene Name] AND Kenya[Location]")
record = Entrez.read(handle)
handle.close()

id_list = record["IdList"]
handle = Entrez.efetch(db="nucleotide", id=",".join(id_list), rettype="fasta", retmode="text")
sequences = handle.read()
handle.close()

with open("kenya_spike.fasta", "w") as file:
    file.write(sequences)
