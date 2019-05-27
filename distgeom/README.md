# distgeom
## 2019/05/27
Conducted initial evaluation (9cf8b4).
Evaluating all molecules in platinum dataset took too long time, so evaluated only 1000.

```
$ time obabel platinum1000.smi -O platinum1000-ob-dg.sdf --gen3d dg
real	4m14.312s
user	4m13.651s
sys	0m0.601s
$ python evaluation.py platinum1000-ob-dg.sdf
same: 939/1000, rmsd: 2.1238465281154926
```
