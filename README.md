# weighted-diameter

(Author : Laurent Viennot, Inria 2023)

Efficient computation of the diameter, radius, eccentricities of a weighted graph in C++

## Compile

```
make
```

## Examples
```
./mdiam -weighted -diam test_digraph.txt   # diameter of largest strongly connected component
cat test_digraph.txt | cut -d' ' -f 1,2 | ./mdiam -symmetrize -simple -rad -   # radius of unweighted version symmetrized
cat test_digraph.txt | ./mdiam -weighted -closeness-all -no-cols -   # compute closeness out-centrality of all nodes
```

## TODO

Fix issues with locks for the multithreaded parts (use `-n-thread 8` without any guarantee).

## License

GNU LGPL.

