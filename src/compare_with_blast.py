"""writes comparison of blast's results with our results"""
import os
import json
import subprocess
import re
import csv

# from scipy.stats import spearmanr, kendalltau
from scipy.special import comb
from scipy.stats import kendalltau, spearmanr
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

    final = dict()
    error_count = 0

    with open(run_variables["query_model_out"], "r") as match_file:
        match_data = json.load(match_file)

    for seq in match_data:
        final[seq] = dict()
        fname = f"{seq}.out"
        # TODO I am not sure which file is the most important, I am chose .pin randomly
        if not os.path.exists("exp_protein_sequences.fasta.pin"):
            os.system(
                'makeblastdb -in "exp_protein_sequences.fasta" -dbtype prot -parse_seqids'
            )
        # if we have not run blast for this sequence before, do so
        if not os.path.exists(f"blast_files/{fname}"):
            # this is nasty and should be scripted differently
            # -query wants to take in a file, but I don't want to create a file
            # file substitution from a python subprocess is gross, but also so cool.
            # TODO see if there is a better way to do this
            subprocess.call(
                [
                    "bash",
                    "-c",
                    f"blastp -db exp_protein_sequences.fasta -query <(echo '> {seq}\n{original_seq}') -out blast_files/{fname} -outfmt 7",
                ]
            )
        with open(f"blast_files/{fname}", "r") as f:
            blast_table = csv.reader(f, delimiter="\t")
            for row in blast_table:
                if row[0].startswith("#") or row[1] == seq:
                    continue
                if match_data[seq].get(row[1]):
                    final[seq].update(
                        {
                            str(row[1]): {
                                "seq_sim": row[2],
                                "cos_sim": match_data[seq].get(row[1], -1),
                                "bit_score": row[-1],
                                "e_value": row[-2],
                            }
                        }
                    )
                else:
                    error_count += 1

    seq_sims = list()
    cos_sims = list()
    bit_scores = list()
    e_values = list()
    seq_sims_ge80 = list()
    cos_sims_ge80 = list()
    seq_sims_le20 = list()
    cos_sims_le20 = list()

    for seq in final:
        for subseq in final[seq]:
            seq_sims.append(float(final[seq][subseq]["seq_sim"]))
            cos_sims.append(float(final[seq][subseq]["cos_sim"]))
            bit_scores.append(float(final[seq][subseq]["bit_score"]))
            e_values.append(float(final[seq][subseq]["e_value"]))
            if float(final[seq][subseq]["seq_sim"]) >= 80:
                seq_sims_ge80.append(float(final[seq][subseq]["seq_sim"]))
                cos_sims_ge80.append(float(final[seq][subseq]["cos_sim"]))
            elif float(final[seq][subseq]["seq_sim"]) <= 20:
                seq_sims_le20.append(float(final[seq][subseq]["seq_sim"]))
                cos_sims_le20.append(float(final[seq][subseq]["cos_sim"]))

    print(kendalltau(seq_sims, cos_sims))
    print(spearmanr(seq_sims, cos_sims))

    seq_lims = get_lims(seq_sims)
    cos_lims = get_lims(cos_sims)

    n_bins = 25
    fig, axs = plt.subplots(3, 3, tight_layout=True)
    axs[0][0].hist(seq_sims, bins=n_bins)
    axs[0][0].set_xlim(seq_lims)
    axs[0][0].set_xlabel("Seq Simularity")
    axs[0][0].set_ylabel("Frequency")
    axs[0][1].hist(cos_sims, bins=n_bins)
    axs[0][1].set_xlim(cos_lims)
    axs[0][1].set_xlabel("Cosine Similarity")
    axs[0][1].set_ylabel("Frequency")
    axs[0][2].set_xlim(seq_lims)
    axs[0][2].set_ylim(cos_lims)
    axs[0][2].set_xlabel("Seq Simularity")
    axs[0][2].set_ylabel("Cosine Similarity")
    axs[0][2].scatter(seq_sims, cos_sims, marker=".")

    seq_lims = get_lims(seq_sims_ge80)
    cos_lims = get_lims(cos_sims_ge80)
    axs[1][0].hist(seq_sims_ge80, bins=n_bins)
    axs[1][0].set_xlim(seq_lims)
    axs[1][0].set_xlabel("Seq Simularity >=80")
    axs[1][0].set_ylabel("Frequency")
    axs[1][1].hist(cos_sims_ge80, bins=n_bins)
    axs[1][1].set_xlim(cos_lims)
    axs[1][1].set_xlabel("Cos Sim (Seq >= 80)")
    axs[1][1].set_ylabel("Frequency")
    axs[1][2].set_xlim(seq_lims)
    axs[1][2].set_ylim(cos_lims)
    axs[1][2].set_xlabel("Seq Simularity >=80")
    axs[1][2].set_ylabel("Cos Sim (Seq >= 80)")
    axs[1][2].scatter(seq_sims_ge80, cos_sims_ge80, marker=".")

    seq_lims = get_lims(seq_sims_le20)
    cos_lims = get_lims(cos_sims_le20)
    axs[2][0].hist(seq_sims_le20, bins=n_bins)
    axs[2][0].set_xlim(seq_lims)
    axs[2][0].set_xlabel("Seq Simularity <= 20")
    axs[2][0].set_ylabel("Frequency")
    axs[2][1].hist(cos_sims_le20, bins=n_bins)
    axs[2][1].set_xlim(cos_lims)
    axs[2][1].set_xlabel("Cos Sim (Seq <= 20)")
    axs[2][1].set_ylabel("Frequency")
    axs[2][2].set_xlim(seq_lims)
    axs[2][2].set_ylim(cos_lims)
    axs[2][2].set_xlabel("Seq Simularity <= 20")
    axs[2][2].set_ylabel("Cos Sim (Seq <= 20)")
    axs[2][2].scatter(seq_sims_le20, cos_sims_le20, marker=".")

    # seq_lims = get_lims(seq_sims)
    # bit_lims = get_lims(bit_scores)
    # axs[3][0].hist(seq_sims, bins=n_bins)
    # axs[3][0].set_xlim(seq_lims)
    # axs[3][0].set_xlabel("Seq Simularity")
    # axs[3][0].set_ylabel("Frequency")
    # axs[3][1].hist(bit_scores, bins=n_bins)
    # axs[3][1].set_xlim(bit_lims)
    # axs[3][1].set_xlabel("Bit Scores")
    # axs[3][1].set_ylabel("Frequency")
    # axs[3][2].set_xlim(seq_lims)
    # axs[3][2].set_ylim(bit_lims)
    # axs[3][2].set_xlabel("Seq Simularity")
    # axs[3][2].set_ylabel("Bit Scores")
    # axs[3][2].scatter(seq_sims, bit_scores, marker=".")

    # seq_lims = get_lims(seq_sims)
    # e_lims = get_lims(e_values)
    # axs[4][0].hist(seq_sims, bins=n_bins)
    # axs[4][0].set_xlim(seq_lims)
    # axs[4][0].set_xlabel("Seq Simularity")
    # axs[4][0].set_ylabel("Frequency")
    # axs[4][1].hist(e_values, bins=n_bins)
    # axs[4][1].set_xlim(e_lims)
    # axs[4][1].set_xlabel("E Values")
    # axs[4][1].set_ylabel("Frequency")
    # axs[4][2].set_xlim(seq_lims)
    # axs[4][2].set_ylim(e_lims)
    # axs[4][2].set_xlabel("Seq Simularity")
    # axs[4][2].set_ylabel("E Values")
    # axs[4][2].scatter(seq_sims, e_values, marker=".")

    print("Errors:", error_count)

    fig.savefig(run_variables["plotted_data"])


if __name__ == "__main__":
    main()
