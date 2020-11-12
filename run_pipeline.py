from os import path, mkdir, stat
from src import (
    compare_with_blast,
    create_kmers,
    create_model,
    query_model,
    normalize_seq,
)

run_variables = {"kmer_len": 4, "extra_info": None}

run_var_formatted = "_".join([
    f"{var}-{run_variables[var]}"
    for var in run_variables
    if run_variables[var] is not None
])
output_dir = f"data_{run_var_formatted}/"
if not path.exists(output_dir):
    mkdir(output_dir)

# format data for model
run_variables["fasta_file"] = "exp_protein_sequences.fasta"
run_variables["create_kmers_out"] = f"{output_dir}kmers.csv"
if not path.exists(run_variables["create_kmers_out"]):
    print("Running: create_kmers.py")
    create_kmers.main(run_variables)
else:
    print("Skipping: create_kmers.py")

# create model
run_variables["create_model_out"] = f"{output_dir}w2v_model.csv"
if not path.exists(run_variables["create_model_out"]):
    print("Running: create_model.py")
    create_model.main(run_variables)
else:
    print("Skipping: create_model.py")

# create a normalized vector for each seq
run_variables["normalize_seq_out"] = f"{output_dir}id_norm_vec.npy"
run_variables["normalize_seq_err"] = f"{output_dir}normalize_seq.errors"
if not path.exists(run_variables["normalize_seq_out"]):
    print("Running: normalize_seq.py")
    normalize_seq.main(run_variables)
else:
    print("Skipping: normalize_seq.py")

# use model to create matches
run_variables["query_model_out"] = f"{output_dir}matches.json"
run_variables["query_model_err"] = f"{output_dir}query_model.errors"
if not path.exists(run_variables["query_model_out"]) or stat(run_variables["query_model_out"]).st_size == 0:
    print("Running: query_model.py")
    query_model.main(run_variables)
else:
    print("Skipping: query_model.py")

# compare with blast's results
# NOTE: blast must be installed locally
# TODO: Write checks for blast and R and and R packages that we are using
run_variables["compare_with_blast"] = f"{output_dir}matches_compared.json"
run_variables["plotted_data"] = f"{output_dir}plotted_data.png"
# if not path.exists(run_variables["compare_with_blast"]):
if not path.exists(run_variables["plotted_data"]) or stat(run_variables["plotted_data"]).st_size == 0:
    print("Running: compare_with_blast.py")
    compare_with_blast.main(run_variables)
else:
    print("Skipping: compare_with_blast.py")
