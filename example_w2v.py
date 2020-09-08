#!/usr/bin/env python

import os
import sys
import logging

import gensim
from gensim.models import Word2Vec
from gensim.models.word2vec import LineSentence
import joblib

assert gensim.models.word2vec.FAST_VERSION > -1

logging.basicConfig(format="%(asctime)s : %(levelname)s : %(message)s", level=logging.INFO)

n_cores = 4
seed = 423
# seed = 564 #original seed

k = 4
d = 256
w = 50
neg_samps = 100
samp_freq = 1e-06
n_min = 10

epochs = 5

name = "uniprot_sprot"

# This is not called...?
# ids_fn = f"{name}_{k}_ids.pkl"
kmers_fn = f"{name}_{k}_kmers.csv"
model_fn = f"w2v_model_{k}_{d}_{epochs}_{w}_{neg_samps}_{str(samp_freq).replace('0.','')}_{n_min}_model.pkl"

print(kmers_fn, model_fn, flush=True)

kmers_init = LineSentence(kmers_fn, max_sentence_length=100000)

model = Word2Vec(
    kmers_init,
    sg=1,
    size=d,
    window=w,
    min_count=n_min,
    negative=neg_samps,
    sample=samp_freq,
    iter=epochs,
    workers=n_cores,
    seed=seed,
)

model.save(model_fn)

w2v_model = {}
for item in model.wv.vocab:
    w2v_model[item] = model[item]
joblib.dump(w2v_model, f"w2v_k{k}.joblib")
