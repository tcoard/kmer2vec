import random
import gensim
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

# get random sentences.
# compare them to the rest of the sentences
# convert top 10 sentences to AA and find uniprot code
# convert chosen original sentence and see what it compares to on blast


def kmers_to_seq(kmers):
    # we want kmer 0, 4, 8, etc
    seq = ""
    for kmer in kmers:
        seq += kmer[0]
    seq += kmers[-1][1:4]
    return seq


def avg_seq_vector(words, model, num_features):
    # function to average all words vectors in a given paragraph
    feature_vec = np.zeros((num_features,), dtype="float32")
    nwords = 0

    for word in words:
        nwords = nwords + 1
        feature_vec = np.add(feature_vec, model[word])

    feature_vec = np.divide(feature_vec, nwords)
    return feature_vec


def main():
    model = gensim.models.Word2Vec.load("w2v_model_4_256_5_50_100_1e-06_10_model.pkl")
    kmer_line_nums = []
    kmer_seqs = []
    random.seed(1)
    for i in range(5):
        # random nubmers at seed 1: 17611 8271 33432 15455 64937
        kmer_line_nums.append(random.randrange(68940))

    prev_line = ""
    with open("uniprot_sprot_4_kmers.csv", "r") as f:
        for i, line in enumerate(f):
            if i in kmer_line_nums:
                kmer_seqs.append({"seq": line, "id": )
            prev_line = ""

    for kmer_seq in kmer_seqs:
        test_kmer_avg_vector = avg_seq_vector(kmer_seq.split(), model=model, num_features=256)
        with open("uniprot_sprot_4_kmers.csv", "r") as f, open("log.errors", 'w') as error:
            num_scores_to_save = 100
            scores = [{"score": 0}] * num_scores_to_save
            for i, line in enumerate(f):
                try:
                    kmer_2_avg_vector = avg_seq_vector(line.split(), model=model, num_features=256)
                # this will get key errors when it encounters kmers that it has not saved
                # TODO see which one's it is not saving and look into why
                except KeyError as e:
                    print(e, file=error)
                    continue

                try:
                    kmer_similarity = cosine_similarity(
                        test_kmer_avg_vector.reshape(1, -1), kmer_2_avg_vector.reshape(1, -1)
                    )[0][0]
                except ValueError as e:
                    print(e, file=error)

                if kmer_similarity > scores[0]["score"]:
                    scores.append({"score": kmer_similarity, "seq": kmers_to_seq(line.split())})
                    # this is super inefficient, but it should work and developer time > compute time
                    # sort the list and then remove the smallest score
                    scores = sorted(scores, key=lambda k: k["score"])[1:num_scores_to_save]

                # if i > 20:
                #     break

            print({"original_seq": kmers_to_seq(kmer_seq.split()), f"top_{num_scores_to_save}_scores": scores})


if __name__ == "__main__":
    main()
