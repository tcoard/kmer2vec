"""writes comparison of blast's results with our results"""
import os
import json


def main():
    num_to_save = 10
    final = dict()
    with open("our_matches_formatted.json") as our_match_file:
        data = json.load(our_match_file)
        for seq in data:
            seq_id = seq["original_id"]
            final[seq_id] = dict()
            fname = f"{seq_id}.out"
            if not os.path.exists(f"{fname}_1.json"):
                os.system(
                    f"blastp -db exp_protein_sequences.fasta -query Q64331.fsa -out blast_files/{fname} -outfmt 13"
                )
            with open(f"blast_files{fname}_1.json", "r") as f:
                data = json.load(f)
                hits = data["BlastOutput2"]["report"]["results"]["search"]["hits"]
                final[seq_id]["blast_matches"] = [
                    hit["description"][0]["accession"] for hit in hits
                ][0 : num_to_save + 1]
                final[seq_id]["our_matches"] = []
                final[seq_id]["in_both"] = []
                for i, match in enumerate(seq["matches"]):
                    # TODO there is a better way to do this
                    if i < num_to_save:
                        # break
                        our_id = match["id"]
                        final[seq_id]["our_matches"].append(our_id)
                    if our_id in final[seq_id]["blast_matches"]:
                        final[seq_id]["in_both"].append(our_id)

        print(json.dumps(final))


if __name__ == "__main__":
    main()
