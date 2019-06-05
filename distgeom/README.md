# distgeom
## 2019/06/05
Finished implementing BFGS, but not working.

1. initial coordinate generation is wrong. Calculated coordinates do not yield correct distance matrix. e.g.
```
$ obabel -:"C" -osdf --gen3d dg
distMat:
      0 1.06214   1.077 1.05546 1.07546
1.06214       0 1.76926 1.75836 1.72459
  1.077 1.76926       0 1.73532 1.73928
1.05546 1.75836 1.73532       0 1.75684
1.07546 1.72459 1.73928 1.75684       0
(omission)
coord:
 -0.00213807
    0.769541
   -0.600926
-5.06441e-09
  -0.0176648
   -0.397392
    0.654739
-6.69999e-09
  -0.0050872
   -0.622523
   -0.560298
 3.26094e-09
   -0.564044
   -0.551822
   -0.541658
-4.66617e-09
    -1.12794
    0.283521
    0.264277
 2.43322e-09
generated distance matrix
       0  1.71425  1.39266   1.4371  1.50074
 1.71425        0  1.23578  1.32429  1.35971
 1.39266  1.23578        0 0.563718  1.66181
  1.4371  1.32429 0.563718        0  1.29047
 1.50074  1.35971  1.66181  1.29047        0
 ```

2. BFGS did not converge. I believe this is partially because of wrong initial coordinates, so I should solve 1 first. 

## 2019/05/29
Try approach A (use distance geometry for out-of-list fragments).

Distance geometry for `O=C1OC(C(=O)N1)C` did not stop, so generation did not stop.

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
