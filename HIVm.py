from Bio import Entrez, SeqIO
import os
import matplotlib.pyplot as plt
import pandas as pd

Entrez.email = "your_email@example.com"

#Create script to retrieve sequences

def fetch_sequences(query, retmax=10):
    handle = Entrez.esearch(db="nucleotide", term=query, retmax=retmax)
    record = Entrez.read(handle)
    handle.close()
    id_list = record["IdList"]
    sequences = []
    for seq_id in id_list:
        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="gb", retmode="text")
        seq_record = SeqIO.read(handle, "genbank")
        sequences.append(seq_record)
        handle.close()
    return sequences

hiv_sequences = fetch_sequences("HIV-1[organism]", retmax=5)
SeqIO.write(hiv_sequences, "hiv_sequences.fasta", "fasta")

#align script 

def align_sequences(input_file, output_file):
    os.system(f"clustalo -i {input_file} -o {output_file} --force --outfmt=clu")

align_sequences("hiv_sequences.fasta", "aligned_hiv_sequences.clustal")

#visualize script

def plot_mutations(mutations):
    df = pd.DataFrame(mutations, columns=["Sequence", "Position", "Reference Base", "Mutant Base"])
    plt.figure(figsize=(10, 6))
    plt.hist(df["Position"], bins=50, alpha=0.75)
    plt.xlabel("Position in Sequence")
    plt.ylabel("Number of Mutations")
    plt.title("Mutation Frequency Across HIV Genome")
    plt.show()

plot_mutations(mutations)

