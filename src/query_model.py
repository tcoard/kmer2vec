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


def main(run_variables):
    num_seq_comparing = 100
    random.seed(1)
    # num_scores_to_save = 100
    # flat[0] because np.save saves as an array
    id_norm_vec = np.load(run_variables["normalize_seq_out"], allow_pickle=True).flat[0]
    random_ids = random.sample(list(id_norm_vec), num_seq_comparing)
    final = dict()
    with open(run_variables["query_model_out"], "w") as out_file, open(
        run_variables["query_model_err"], "w"
    ) as error:
        # printing the json structure of this document piecemeal so that if it stops early, we still have partial data
        for i, random_id in enumerate(random_ids):
            # scores = [{"score": 0}] * num_scores_to_save
            scores = dict()
            random_vec = id_norm_vec[random_id]["vec"]
            random_seq = id_norm_vec[random_id]["seq"]
            # if i > 5:
            #     break
            for j, vec_id in enumerate(id_norm_vec):
                # if j > 5:
                #     break
                vec = id_norm_vec[vec_id]["vec"]
                seq = id_norm_vec[vec_id]["seq"]
                try:
                    kmer_similarity = cosine_similarity(random_vec, vec)[0][0]
                except ValueError as err:
                    print(err, file=error)

                scores[vec_id] = float(kmer_similarity)
                # if kmer_similarity > scores[-1]["score"]:
                #     scores.append({"score": kmer_similarity, "id": vec_id, "seq": seq})
                #     # TODO this is super inefficient, but it should work and developer time > compute time
                #     # sort the list and then remove the smallest score
                #     scores = sorted(scores, key=lambda k: k["score"], reverse=True)[
                #         0 : num_scores_to_save - 1
                #     ]
            final[random_id] = scores

        print(json.dumps(final, indent=2), file=out_file)

if __name__ == "__main__":
    main()
