import json


def find_id(seq):
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


def main():
    num_scores_to_save = 6
    ending_json = list()
    with open("formatted.json") as f:
        data = json.load(f)
        for seq in data:
            filtered_seqs = sorted(
                seq["top_100_scores"], key=lambda k: k["score"], reverse=True
            )[0 : num_scores_to_save + 1]
            for filtered_seq in filtered_seqs:
                filtered_seq["id"] = find_id(filtered_seq["seq"])

            ending_json.append(
                {"original_id": find_id(seq["original_seq"]), "matches": filtered_seqs}
            )
        print(json.dumps(ending_json))


if __name__ == "__main__":
    main()
