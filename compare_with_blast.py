import os
import json

fname = "13.out"
if not os.path.exists("{fname}_1.json"):
    os.system(f"blastp -db exp_protein_sequences.fasta -query Q64331.fsa -out {fname} -outfmt 13")
with open("shortened_formatted.json") as w2v_file:
    for 
with open(f"{fname}_1.json", 'r') as f:
    data = json.load(f)
    hits = data["BlastOutput2"]["report"]["results"]["search"]["hits"]
    blast_matches = [hit["description"][0]["accession"] for hit in hits]

    
