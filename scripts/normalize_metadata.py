"""Put the AGP metadata in a fixed order so its md5 is stable.

redbiom can write metadata columns in a different order on each fetch.
Sorting columns and rows (values unchanged) makes the file byte-identical.
Usage: python scripts/normalize_metadata.py data/raw/agp/agp_metadata.tsv
"""
import sys

import pandas as pd

path = sys.argv[1]
df = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, index_col=0)
df = df.sort_index().sort_index(axis=1)
df.to_csv(path, sep="\t")
print(f"normalized {path}: {df.shape[0]} samples x {df.shape[1]} columns")
