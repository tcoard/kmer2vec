"""writes comparison of blast's results with our results"""
import os
import json
import subprocess
import re
from scipy.stats import spearmanr, kendalltau
from scipy.special import comb


def rank_data(their_data, our_data):
    rank = dict()
    for pos, val in enumerate(their_data):
        rank[val] = pos

    our_ranked_data = list()
    discordant = 0.0
    concordant = 0.0
    for pos, val in enumerate(our_data):
        our_ranked_data.append(rank.get(val, "NA"))
        for their_val in their_data[pos+1:]:
            if not rank.get(val) or rank[val] > rank[their_val]:
                discordant += 1
            else:
                concordant += 1
    # our data and their data are the same length
    kendall_t = (concordant - discordant) / comb(len(their_data), 2)
    return our_ranked_data, kendall_t


def main():
    max_num_to_save = 50
    final = dict()
    with open("matches.json", "r") as match_file, open("matches_compared.json", "w") as out_file:
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
                # TODO see if there is a better way to do this
                subprocess.call(
                    [
                        "bash",
                        "-c",
                        f"blastp -db exp_protein_sequences.fasta -query <(echo '> {orig_seq_id}\n{original_seq}') -out blast_files/{fname} -outfmt 13",
                    ]
                )
            with open(f"blast_files/{fname}_1.json", "r") as f:

                blast_data = json.load(f)
                hits = blast_data["BlastOutput2"]["report"]["results"]["search"]["hits"]
                num_to_save = min(len(hits), max_num_to_save)

                blast_ids = [hit["description"][0]["accession"] for hit in hits][0:num_to_save]
                our_ids = [seq["id"] for seq in seq["scores"]][0:num_to_save]

                our_ranked_data, kendall_t = rank_data(blast_ids, our_ids)
                blast_ranked_data = list(range(len(blast_ids)))
                chars_to_remove = ['"', "'", "[", "]"]
                blast_ranked_data = str(blast_ranked_data)
                our_ranked_data = str(our_ranked_data)
                for char in chars_to_remove:
                    blast_ranked_data = blast_ranked_data.replace(char, "")
                    our_ranked_data = our_ranked_data.replace(char, "")

                captured = subprocess.run(
                    [
                        "Rscript",
                        "-e",
                        f"library(irrNA); data <- data.frame(c({blast_ranked_data}), c({our_ranked_data})); kendallNA(data)",
                    ],
                    capture_output=True,
                )
                kendall_w = re.findall(r'\$([^\n]+)\n+\[1\] ([^\n]+)', captured.stdout.decode("utf-8"))
               # final[orig_seq_id]["kendall"] = dict()
                # for i, j in kendall_w:
                #     final[orig_seq_id]["kendall"][i] = j
                # final[orig_seq_id]["kendall"]["tau"] = kendall_t

                final[orig_seq_id]["kendall"]= kendall_t

                # final[orig_seq_id]["blast_matches"] = blast_ids
                # final[orig_seq_id]["our_matches"] = our_ids
                # turning it back to list because json does not like sets
                # final[orig_seq_id]["in_both"] = list(set(blast_ids).intersection(set(our_ids)))
                final[orig_seq_id]["in_both/total"] = [round(float(len(set(blast_ids).intersection(set(our_ids))))/float(num_to_save), 3), f"{len(set(blast_ids).intersection(set(our_ids)))}/{num_to_save}"]

                # rho, s_pval = spearmanr(blast_ids, our_ids)
                # tau, k_pval = kendalltau(blast_ids, our_ids)
                # final[orig_seq_id]["spearman"] = dict()
                # final[orig_seq_id]["spearman"]["rho"] = rho
                # final[orig_seq_id]["spearman"]["pval"] = s_pval
                # final[orig_seq_id]["kendall"]["tau"] = tau
                # final[orig_seq_id]["kendall"]["pval"] = k_pval


        print(json.dumps(final, indent=2), file=out_file)


if __name__ == "__main__":
    main()
