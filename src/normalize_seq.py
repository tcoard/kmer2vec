"""
"""
import warnings
import gensim
import numpy as np

# some of our imports have a lot of deprication warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def find_id(seq):
    """Search through the original Fasta for the sequence's ID"""
    with open("exp_protein_sequences.fasta") as fasta:
        prev_line = ""
        for line in fasta:
            line = line.strip()
            if line == seq:
                break
            prev_line = line
        else:
            raise ValueError(f"could not find {seq}")

        return prev_line.replace(">", "")


def kmers_to_seq(kmers):
    """convers kmers to a sequence
    Example: ["ACDEF", "CDEFG", "DEFGH"] becomes "ACDEFGH"
    """
    seq = ""
    for kmer in kmers:
        seq += kmer[0]
    try:
        seq += kmers[-1][1:4]
    except IndexError as e:
        breakpoint()
    return seq


def avg_seq_vector(words, model, num_features):
    """averages all word (kmers) vectors in a given paragraph (sequence)"""
    feature_vec = np.zeros((num_features,), dtype="float32")

    for word in words:
        feature_vec = np.add(feature_vec, model[word])

    feature_vec = np.divide(feature_vec, len(words))
    return feature_vec


def main(run_variables):
    """normalize sequences and make a ID -> normalized vector dict"""
    model = gensim.models.Word2Vec.load("run_variables['create_model_out']")
    # TODO make a better name for this var
    id_norm_vec = dict()
    # get all of the kmer-ed sequences
    with open(run_variables["create_kmers_out"], "r") as kmer_file, open(
        run_variables["normalize_seq_out"], "w"
    ) as error:
        for i, kmer_seq in enumerate(kmer_file):
            kmer_seq = kmer_seq.split()
            if not kmer_seq:
                print(f"No seq on line: {i}", file=error)
                continue
            try:
                kmer_2_avg_vector = avg_seq_vector(
                    kmer_seq, model=model, num_features=256
                )
            # this will get key errors when it encounters kmers that it has not saved
            except KeyError as err:
                print(err, file=error)
                continue
            seq = kmers_to_seq(kmer_seq)
            id_norm_vec[find_id(seq)] = {
                "vec": kmer_2_avg_vector.reshape(1, -1),
                "seq": seq,
            }

    np.save(run_variables["normalize_seq_vec"], id_norm_vec)


if __name__ == "__main__":
    main()
