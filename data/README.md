# Construction data

`constructions.json` is the compact, human-readable record of all three even-party code constructions and the projection-derived odd-party states. The coordinate order for the two group-circulant matrices is explicitly fixed.

`huber-wyderka-table-rows_2026-08-01.json` is the archived snapshot of the Huber-Wyderka AME table used by Fig. 7 and by the historical claims in the comparison section.

Each `theta_*.npy` file is accompanied by a plain-text integer matrix. The stored 18-party arrays may retain diagonal quadratic phase terms produced by the local-Clifford row reduction; balanced cut blocks do not use the diagonal. The 17-party arrays are the recorded principal deletions used by the archived derivations. The verification scripts also check every one of the 18 possible deletions for q=11 and q=13.

Integrity is checked mathematically rather than by hashes: `verification/check_consistency.py` confirms that every copy of every matrix in this directory agrees with the copies embedded in the verifiers and re-establishes the balanced-cut certificates directly from the plain-text arrays.
