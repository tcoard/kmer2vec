import json
num_scores_to_save = 6
ending_json = list()
with open("formatted.json") as f:
    data = json.load(f)
    for seq in data:
        filtered_seqs = sorted(seq["top_100_scores"], key=lambda k: k["score"], reverse=True)[0:num_scores_to_save+1]
        for filtered_seq in filtered_seqs:
            with open("exp_protein_sequences.fasta") as fasta:
                prev_line = "" 
                for line in fasta:
                    line = line.strip()
                    if line == filtered_seq["seq"]:
                        break
                    prev_line = line
                else:
                    raise ValueError(f"could not find {filtered_seq['seq']}")
                filtered_seq["id"] = prev_line.replace('>', '')

        ending_json.append(filtered_seqs)
    print(json.dumps(ending_json))

