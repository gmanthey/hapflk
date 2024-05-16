# Welcome to hapFLK documentation


hapflk is a software implementing the hapFLK [1] and FLK [2] tests for
the detection of selection signatures based on multiple population
genotyping data.

[1] [Fariello et al., 2013, Detecting Signatures of Selection Through
Haplotype Differentiation Among Hierarchically Structured
Populations. Genetics
193(3):929-941.](http://www.genetics.org/content/193/3/929.abstract)

[2] [Bonhomme et al., 2010, Detecting selection in population trees:
The Lewontin and Krakauer test extended. Genetics 186(1)
241-262](http://www.genetics.org/content/186/1/241.abstract)

# Installing and running hapflk

`hapflk` is distributed as a python package and is simply installed using pip:

```bash
pip install hapflk
```

Of course, it is always a good idea to do that within a python virtual environment to not mingle with other packages. 

Once this is done you should be able to run the software and get something like this:

```console
$ hapflk
          Start @ 2024-05-16 15:53:58
usage: hapflk [-h] (--bfile PREFIX | --sfile PREFIX) [-p PREFIX] [--ncpu N] [--eigen] [--reynolds] [--kinship FILE]
              [--reynolds-snps L] [--outgroup POP] [--keep-outgroup] [--covkin] [-K K] [--nfit NFIT] [--phased] [--kfrq] [--legacy]
              [--annot]

options:
  -h, --help            show this help message and exit
  --bfile PREFIX        PLINK bfile prefix (bim,fam,bed) (default: None)
  --sfile PREFIX        ShapeIT file prefix (haps,sample) (default: None)
  -p PREFIX, --prefix PREFIX
                        prefix for output files (default: hapflk)
  --ncpu N              Use N processors when possible (default: 40)
  --eigen               Perform eigen decomposition of tests (default: False)
  --reynolds            Force writing down Reynolds distances (default: False)
  --annot               Shortcut for --eigen --reynolds --kfrq (default: False)

Population kinship :
  Set parameters for getting the population kinship matrix

  --kinship FILE        Read population kinship from file (if None, kinship is estimated) (default: None)
  --reynolds-snps L     Number of SNPs to use to estimate Reynolds distances (default: 100000)
  --outgroup POP        Use population POP as outgroup for tree rooting (if None, optimize root location) (default: None)
  --keep-outgroup       Keep outgroup in population set (default: False)
  --covkin              Use covariance matrix as kinship (default: False)

hapFLK and LD model:
  Switch on hapFLK calculations and set parameters of the LD model

  -K K                  Set the number of clusters to K. hapFLK calculations switched off if K<0 (default: -1)
  --nfit NFIT           Set the number of model fit to use (default: 10)
  --phased, --inbred    Haplotype data provided (default: False)
  --kfrq                Write Cluster frequencies (Big files) (default: False)
  --legacy              Use Legacy fastPHASE (default: False)

```



The meaning of all this is explained in the documentation which you can navigate with the links on the left.
