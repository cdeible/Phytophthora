import pandas as pd
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import combinations
from Levenshtein import distance


def get_species_count(filename):
    df = pd.read_csv(filename)
    data = df["Species"].tolist()

    fams_cleaned = [re.sub(r'^.*Phytophthora\s*\.?\s*(sp\.|taxon|x)?\.?\s*', '', s, flags=re.IGNORECASE) for s in data]
    bio = SeqIO.parse("noRepeat.fasta", "fasta")
    data = {}

    for bio_data in bio:
        header = bio_data.description.split(' ')
        for word in header:
            if word in fams_cleaned:
                if word in data:
                    data[word] += 1
                else:
                    data[word] = 1
                break
    return data

#print(get_species_count("speciesCount.csv"))

def clean_species(species_name):
    gene_things = []
    bio = SeqIO.parse("noRepeat.fasta", "fasta")
    for bio_data in bio:
        header = bio_data.description.split(' ')
        if species_name in header:
            gene_things.append(bio_data)
    return gene_things

def calculate_pairwise_distances(sequences):
    """Calculates the Levenshtein distance between all pairs of sequences."""
    distances = {}
    for (id1, seq1), (id2, seq2) in combinations([(record.id, str(record.seq)) for record in sequences], 2):
        distances[(id1, id2)] = distance(seq1, seq2)
    return distances

def select_most_different(sequences, n): # Note: will only return an even number, so will add one to odd numbers
    """Selects the n most different sequences based on pairwise Levenshtein distance."""
    distances = calculate_pairwise_distances(sequences)
    
    # Sort distances in descending order
    sorted_distances = sorted(distances.items(), key=lambda item: item[1], reverse=True)
    
    selected_ids = set()
    selected_sequences = []
    
    for (id1, id2), dist in sorted_distances:
        if len(selected_sequences) >= n:
            break
        if id1 not in selected_ids and id2 not in selected_ids:
            seq1 = next(record for record in sequences if record.id == id1)
            seq2 = next(record for record in sequences if record.id == id2)
            selected_sequences.extend([seq1, seq2])
            selected_ids.update([id1, id2])
    
    #If less than n sequences are selected, add remaining sequences until n is reached
    if len(selected_sequences) < n:
        for record in sequences:
            if len(selected_sequences) >= n:
                break
            if record.id not in selected_ids:
                selected_sequences.append(record)
                selected_ids.add(record.id)
    
    return selected_sequences

masterCount = get_species_count("speciesCount.csv")
#names = pd.read_csv("names.txt")
names = ["ramorum"]
keeper = []
for names_data in names: #["Species"]:
    if masterCount[names_data] > 5:     # right now this is set at 5, but change as needed
        seqs = clean_species(names_data)
        keeper.extend(select_most_different(seqs, 5))       # also set at 5, but change as needed
    else:
        keeper.extend(clean_species(names_data))
SeqIO.write(keeper, "selectedSeqs.fasta", "fasta")
