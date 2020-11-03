"""Turns sequences into a list of kmers
Example: "ACDEFGH" becomes ["ACDEF", "CDEFG", "DEFGH"]
"""


def make_kmers(seq, k):
    kmers = []
    for i in range(len(seq) - k + 1):
        kmers.append(seq[i : i + k])
    return kmers


def main(run_variables):
    with open(run_variables["fasta_file"], "r") as in_file, open(
        run_variables["create_kmers_out"], "w"
    ) as out_file:
        for i, line in enumerate(in_file):
            # skipping ID lines for now. They are readded later.
            if ">" in line:
                continue

            kmers = make_kmers(line.strip(), 4)
            print(" ".join(kmers), file=out_file)
            print(i)


if __name__ == "__main__":
    main()
