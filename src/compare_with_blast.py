"""writes comparison of blast's results with our results"""
import os
import json
import subprocess
import re

# from scipy.stats import spearmanr, kendalltau
from scipy.special import comb
from matplotlib import pyplot as plt


def get_lims(data):
    buffer_len = (max(data) - min(data)) * 0.1
    lower_lim = min(data) - buffer_len
    upper_lim = max(data) + buffer_len
    return lower_lim, upper_lim


def rank_data(their_data, our_data):
    rank = dict()
    for pos, val in enumerate(their_data):
        rank[val] = pos

    our_ranked_data = list()
    discordant = 0.0
    concordant = 0.0
    for val in our_data:
        ranked = rank.get(val)
        if ranked:
            our_ranked_data.append(ranked)

    their_ranked_data = sorted(our_ranked_data)
    for pos, rank in enumerate(our_ranked_data):
        for their_ranked_val in their_ranked_data[pos + 1 :]:
            if rank > their_ranked_val:
                discordant += 1
            else:
                concordant += 1
    # our data and their data are the same length
    kendall_t = 0
    if float(comb(len(our_ranked_data), 2)) != 0:
        kendall_t = float(concordant - discordant) / float(
            comb(len(our_ranked_data), 2)
        )
    return our_ranked_data, kendall_t


def main(run_variables):
    kendall_t_avg = 0.0
    kendall_w_avg = 0.0
    shared_percent_avg = 0.0

    max_num_to_save = 50
    final = dict()
    with open(run_variables["query_model_out"], "r") as match_file, open(
        run_variables["compare_with_blast"], "w"
    ) as out_file:
        match_data = json.load(match_file)
        our_scores = list()
        bit_scores = list()
        sorted_our_scores_bit = list()
        e_values = list()
        sorted_our_scores_e = list()
        for seq in match_data:
            orig_seq_id = seq["original_id"]
            original_seq = seq["original_seq"]
            final[orig_seq_id] = dict()
            fname = f"{orig_seq_id}.out"
            # TODO I am not sure which file is the most important, I am chose .pin randomly
            if not os.path.exists("exp_protein_sequences.fasta.pin"):
                os.system(
                    'makeblastdb -in "exp_protein_sequences.fasta" -dbtype prot -parse_seqids'
                )
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
                # hits - 1 because we will be removing the original ID from the search result
                num_to_save = min(len(hits) - 1, max_num_to_save)
                if num_to_save == 0:
                    continue
                bit_id_score = dict()
                e_id_score = dict()
                for hit in hits:
                    blast_id = hit["description"][0]["accession"]
                    if orig_seq_id != blast_id:
                        bit_id_score[blast_id] = hit["hsps"][0]["bit_score"]
                        e_id_score[blast_id] = hit["hsps"][0]["evalue"]
                        bit_scores.append(hit["hsps"][0]["bit_score"])
                        e_values.append(float(hit["hsps"][0]["evalue"]))

                our_id_score = dict()
                for sub_seq in seq["scores"]:
                    if sub_seq["id"] != orig_seq_id:
                        our_id_score[sub_seq["id"]] = sub_seq["score"]
                        our_scores.append(sub_seq["score"])

                for bit_id in bit_id_score:
                    sorted_our_scores_bit.append(our_id_score.get(bit_id, 0))

                for e_id in e_id_score:
                    sorted_our_scores_e.append(our_id_score.get(e_id, 0))

                # our_scores = [seq["id"] for seq in seq["scores"] if seq["id"] != orig_seq_id][0:num_to_save]

                blast_ids = [
                    hit["description"][0]["accession"]
                    for hit in hits
                    if hit["description"][0]["accession"] != orig_seq_id
                ][0:num_to_save]
                our_ids = [
                    sub_seq["id"]
                    for sub_seq in seq["scores"]
                    if sub_seq["id"] != orig_seq_id
                ][0:num_to_save]
                our_ranked_data, kendall_t = rank_data(blast_ids, our_ids)
                # blast_ranked_data = list(range(len(our_ranked_data)))
                blast_ranked_data = sorted(our_ranked_data)
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
                kendall_w = re.findall(
                    r"\$([^\n]+)\n+\[1\] ([^\n]+)", captured.stdout.decode("utf-8")
                )
                final[orig_seq_id]["Kendall's W statistics"] = dict()
                for i, j in kendall_w:
                    final[orig_seq_id]["Kendall's W statistics"][i.strip("`")] = j
                final[orig_seq_id]["Kendall's Tau"] = kendall_t

                # final[orig_seq_id]["kendall"]= kendall_t

                # final[orig_seq_id]["blast_matches"] = blast_ids
                # final[orig_seq_id]["our_matches"] = our_ids
                # turning it back to list because json does not like sets
                # final[orig_seq_id]["in_both"] = list(set(blast_ids).intersection(set(our_ids)))
                perc_shared_decimal = round(
                    float(len(set(blast_ids).intersection(set(our_ids))))
                    / float(num_to_save),
                    3,
                )
                final[orig_seq_id]["in_both/total"] = [
                    perc_shared_decimal,
                    f"{len(set(blast_ids).intersection(set(our_ids)))}/{num_to_save}",
                ]

                kendall_t_avg += float(kendall_t)
                kendall_w_avg += float(
                    final[orig_seq_id]["Kendall's W statistics"].get("Kendall's W", 0)
                )
                shared_percent_avg += float(perc_shared_decimal)

                # rho, s_pval = spearmanr(blast_ids, our_ids)
                # tau, k_pval = kendalltau(blast_ids, our_ids)
                # final[orig_seq_id]["spearman"] = dict()
                # final[orig_seq_id]["spearman"]["rho"] = rho
                # final[orig_seq_id]["spearman"]["pval"] = s_pval
                # final[orig_seq_id]["kendall"]["tau"] = tau
                # final[orig_seq_id]["kendall"]["pval"] = k_pval

        final["Averages"] = dict()
        final["Averages"]["Kendall's W"] = round(
            kendall_t_avg / float(len(match_data)), 3
        )
        final["Averages"]["Kendall's Tau"] = round(
            kendall_w_avg / float(len(match_data)), 3
        )
        final["Averages"]["Shared Percentage"] = round(
            shared_percent_avg / float(len(match_data)), 3
        )
        print(json.dumps(final, indent=2, sort_keys=True), file=out_file)

        # fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
        fig, axs = plt.subplots(2, 3, tight_layout=True)

        # We can set the number of bins with the `bins` kwarg
        bit_lims = get_lims(bit_scores)
        e_lims = get_lims(e_values)
        our_lims = get_lims(our_scores)

        n_bins = 25
        axs[0][0].hist(bit_scores, bins=n_bins)
        axs[0][0].set_xlim(bit_lims)
        axs[0][0].set_xlabel("Blast Bit Score")
        axs[0][0].set_ylabel("Frequency")
        axs[0][1].hist(our_scores, bins=n_bins)
        axs[0][1].set_xlim(our_lims)
        axs[0][1].set_xlabel("Kmer2Vec Cosine Similarity")
        axs[0][1].set_ylabel("Frequency")
        axs[0][2].set_xlim(bit_lims)
        axs[0][2].set_ylim(our_lims)
        axs[0][2].set_xlabel("Blast Bit Score")
        axs[0][2].set_ylabel("Kmer2Vec Cosine Similarity")
        axs[0][2].scatter(bit_scores, sorted_our_scores_bit, marker=".")

        axs[1][0].hist(e_values, bins=n_bins)
        axs[1][0].set_xlim(e_lims)
        axs[1][0].set_xlabel("Blast E Values")
        axs[1][0].set_ylabel("Frequency")
        axs[1][2].set_xlim(e_lims)
        axs[1][2].set_ylim(our_lims)
        axs[1][2].set_xlabel("Blast E Values")
        axs[1][2].set_ylabel("Kmer2Vec Cosine Similarity")
        axs[1][2].scatter(e_values, sorted_our_scores_e, marker=".")
        fig.savefig(run_variables["plotted_data"])


if __name__ == "__main__":
    main()
