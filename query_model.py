"""Finds the cosine similarity between a random subset of sequences and all other sequences
"""

import random
import json
import warnings
import gensim
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# some of our imports have a lot of deprication warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def kmers_to_seq(kmers):
    """convers kmers to a sequence
    Example: ["ACDEF", "CDEFG", "DEFGH"] becomes "ACDEFGH"
    """
    seq = ""
    for kmer in kmers:
        seq += kmer[0]
    seq += kmers[-1][1:4]
    return seq


def avg_seq_vector(words, model, num_features):
    """averages all word (kmers) vectors in a given paragraph (sequence)"""
    feature_vec = np.zeros((num_features,), dtype="float32")

    for word in words:
        feature_vec = np.add(feature_vec, model[word])

    feature_vec = np.divide(feature_vec, len(words))
    return feature_vec


def main():
    model = gensim.models.Word2Vec.load("w2v_model_4_256_5_50_100_1e-06_10_model.pkl")
    kmer_line_nums = []
    kmer_seqs = []
    num_seq_comparing = 5
    random.seed(1)
    for i in range(num_seq_comparing):
        # random nubmers at seed 1: 17611 8271 33432 15455 64937
        # TODO make sure that I am not generating any duplicates
        kmer_line_nums.append(
            random.randrange(68940)
        )  # 68940 is the number of sequences in the original fasta

    # get all of the kmer-ed sequences
    with open("uniprot_sprot_4_kmers.csv", "r") as kmer_file:
        for i, line in enumerate(kmer_file):
            if i in kmer_line_nums:
                kmer_seqs.append(line)

    with open("matches.json", "w") as out_file:
        # printing the json structure of this document piecemeal so that if it stops early, we still have partial data
        print("[", file=out_file)
        for i, kmer_seq in enumerate(kmer_seqs):
            kmer_avg_vector = avg_seq_vector(kmer_seq.split(), model=model, num_features=256)
            with open("uniprot_sprot_4_kmers.csv", "r") as kmer_file, open("log.errors", "w") as error:
                num_scores_to_save = 100
                scores = [{"score": 0}] * num_scores_to_save
                for line in kmer_file:
                    try:
                        kmer_2_avg_vector = avg_seq_vector(line.split(), model=model, num_features=256)
                    # this will get key errors when it encounters kmers that it has not saved
                    except KeyError as err:
                        print(err, file=error)
                        continue

                    try:
                        # reshaping the vector here because cosine_similarity was
                        # complaining about our vector's dimensions
                        # I _think_ this if fine but,
                        # TODO I should double check the reshape
                        kmer_similarity = cosine_similarity(
                            kmer_avg_vector.reshape(1, -1),
                            kmer_2_avg_vector.reshape(1, -1),
                        )[0][0]
                    except ValueError as err:
                        print(err, file=error)

                    if kmer_similarity > scores[-1]["score"]:
                        scores.append(
                            {
                                "score": kmer_similarity,
                                "seq": kmers_to_seq(line.split()),
                            }
                        )
                        # TODO this is super inefficient, but it should work and developer time > compute time
                        # sort the list and then remove the smallest score
                        scores = sorted(scores, key=lambda k: k["score"], reverse=True)[0 : num_scores_to_save - 1]

                print(
                    json.dumps(
                        {
                            "original_seq": kmers_to_seq(kmer_seq.split()),
                            "scores": json.loads(str(scores).replace("'", '"')),
                        },
                        indent=2,
                    )
                    # adding the comma manually, but not to the last item
                    + ("," if (i + 1) != num_seq_comparing else ""),
                    file=out_file,
                )

        print("]", file=out_file)


if __name__ == "__main__":
    main()
