import os
import json


def main():
    # TODO put in check(s) to make sure that I am actually evaluating the top scores from both
    num_to_save = 10
    final = dict()
    with open("blast_compare_ready.json") as our_match_file:
        data = json.load(our_match_file)
        for seq in data:
            seq_id = seq["original_id"]
            final[seq_id] = dict()
            fname = f"{seq_id}.out"
            if not os.path.exists(f"{fname}_1.json"):
                os.system(
                    f"blastp -db exp_protein_sequences.fasta -query Q64331.fsa -out {fname} -outfmt 13"
                )
            with open(f"{fname}_1.json", "r") as f:
                data = json.load(f)
                hits = data["BlastOutput2"]["report"]["results"]["search"]["hits"]
                final[seq_id]["blast_matches"] = [
                    hit["description"][0]["accession"] for hit in hits
                ][0 : num_to_save + 1]
                final[seq_id]["our_matches"] = []
                final[seq_id]["in_both"] = []
                for i, match in enumerate(seq["matches"]):
                    our_id = match["id"]
                    final[seq_id]["our_matches"].append(our_id)
                    if our_id in final[seq_id]["blast_matches"]:
                        final[seq_id]["in_both"].append(our_id)
                    # TODO there is a better way to do this
                    if i == 9:
                        break

        print(json.dumps(final))


if __name__ == "__main__":
    main()
