"""Adds IDs to the output of our analysis output.

This could be done in the query_model, but because this takes time
and we want query model to be as quick as possible
(so that it is less likely to be interrupted), we are doing this after
"""
import json


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


def main():
    num_scores_to_save = 10
    ending_json = list()
    with open("our_matches.json", "r") as our_matches, open("our_matches_formatted.json", "w") as formatted:
        data = json.load(our_matches)
        # for each sequence, find the top num_scores_to_save result's and their IDs
        for seq in data:
            # TODO this could be done more efficiently, but pragmatically, it won't save much more time
            # keep the top 10 results
            filtered_seqs = sorted(seq["scores"], key=lambda k: k["score"], reverse=True)[0 : num_scores_to_save + 1]
            for filtered_seq in filtered_seqs:
                filtered_seq["id"] = find_id(filtered_seq["seq"])

            ending_json.append(
                {"original_id": find_id(seq["original_seq"]), "matches": filtered_seqs}
            )
        print(json.dumps(ending_json, indent=2), file=formatted)


if __name__ == "__main__":
    main()
