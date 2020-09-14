"""writes comparison of blast's results with our results"""
import os
import json
import subprocess


def main():
    num_to_save = 10
    final = dict()
    with open("matches_formatted.json") as match_file:
        match_data = json.load(match_file)
        for seq in match_data:
            orig_seq_id = seq["original_id"]
            original_seq = seq["original_seq"]
            final[orig_seq_id] = dict()
            fname = f"{orig_seq_id}.out"
            # TODO I am not sure which file is the most important, I am chose .pin randomly 
            if not os.path.exists("exp_protein_sequences.fasta.pin"):
                os.system('makeblastdb -in "exp_protein_sequences.fasta" -dbtype prot -parse_seqids')
            # if we have not run blast for this sequence before, do so
            if not os.path.exists(f"blast_files/{fname}_1.json"):
                # this is nasty and should be scripted differently
                # -query wants to take in a file, but I don't want to create a file
                # file substitution from a python subprocess is gross, but also so cool.
                subprocess.call(["bash", "-c", f"blastp -db exp_protein_sequences.fasta -query <(echo '> {orig_seq_id}\n{original_seq}') -out blast_files/{fname} -outfmt 13"])
            with open(f"blast_files/{fname}_1.json", "r") as f:
                blast_data = json.load(f)
                hits = blast_data["BlastOutput2"]["report"]["results"]["search"]["hits"]
                final[orig_seq_id]["blast_matches"] = [
                    hit["description"][0]["accession"] for hit in hits
                ][0 : num_to_save + 1]
                final[orig_seq_id]["matches"] = []
                final[orig_seq_id]["in_both"] = []
                for i, match in enumerate(seq["matches"]):
                    # TODO there is a better way to do this
                    if i < num_to_save:
                        # break
                        seq_id = match["id"]
                        final[orig_seq_id]["matches"].append(seq_id)
                    if seq_id in final[orig_seq_id]["blast_matches"]:
                        final[orig_seq_id]["in_both"].append(seq_id)

        print(json.dumps(final))


if __name__ == "__main__":
    main()
