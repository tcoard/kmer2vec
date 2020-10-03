import add_id, compare_with_blast, create_kmers, example_w2v, query_model

# TODO don't have all of the data dump to the root dir
# TODO make data file names more descriptive

# format data for model
#create_kmers.main()

# create model
# if .pkl
#example_w2v.main()

# create a normalized vector for each seq
#normalize_seq.py

# use model to create matches
query_model.main()

# compare with blast's results
# NOTE: blast must be installed locally
compare_with_blast.main()
