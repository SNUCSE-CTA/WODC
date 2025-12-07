# Efficient Defective Clique Enumeration and Search with Worst-Case Optimal Search Space (SIGMOD 2026)
This repository implements algorithms for maximal k-defective clique enumeration and maximum k-defective clique search.


## How to compile
### Maximal k-defective clique enumeration
```bash
mkdir build
cd build
cmake .. -DMAXIMAL=ON
make
```
### Maximum k-defective clique search
```bash
mkdir build
cd build
cmake .. -DMAXIMUM=ON
make
```

An executable named `main.out` will be generated in the `build/main` directory.

## How to run 
In the `build` directory:
### Maximal k-defective clique enumeration
```bash
./main/main.out {graph-path} {k} {q}
```

### Maximum k-defective clique search
```bash
./main/main.out {graph-path} {k}
```

For example, to run the maximal k-defective clique enumeration algorithm:
```bash
./main/main.out ../datasets/soc-FourSquare/ 1 30
```

## Data format
Our data format follows the binary format described here:
https://lijunchang.github.io/Cohesive_subgraph_book/datasets



## Datasets
* Real-world graph datasets: https://lcs.ios.ac.cn/~caisw/Resource/realworld%20graphs.tar.gz
* Large-scale graph datasets: https://law.di.unimi.it/datasets.php
