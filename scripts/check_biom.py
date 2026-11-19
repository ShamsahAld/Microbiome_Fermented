"""md5 of the data in a BIOM table, ignoring its embedded creation date.

Hashes the sorted sample IDs, ASV IDs and counts.
  python scripts/check_biom.py TABLE.biom --write RECORD   # record the hash
  python scripts/check_biom.py TABLE.biom RECORD           # verify
"""
import hashlib
import sys

import biom

table_path = sys.argv[1]
t = biom.load_table(table_path).sort(axis="sample").sort(axis="observation")
m = t.matrix_data.tocsr()
h = hashlib.md5()
h.update("\n".join(t.ids("sample")).encode())
h.update("\n".join(t.ids("observation")).encode())
for arr in (m.data, m.indices, m.indptr):
    h.update(arr.tobytes())
digest = h.hexdigest()

if sys.argv[2] == "--write":
    open(sys.argv[3], "w").write(f"{digest}  {table_path} (data only)\n")
    print(f"wrote {sys.argv[3]}")
else:
    want = open(sys.argv[2]).read().split()[0]
    ok = digest == want
    print(f"{table_path}: {'OK' if ok else 'FAILED'} (data only)")
    sys.exit(0 if ok else 1)
