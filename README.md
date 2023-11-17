# weighted-diameter

(Author : Laurent Viennot, Inria 2023)

Efficient computation of the diameter, radius, and even all eccentricities of a weighted graph in C++.

## Compile

```
make
```

## Examples
```
./mdiam -weighted -largest-scc -diam test_digraph.txt   # diameter of largest strongly connected component
cat test_digraph.txt | cut -d' ' -f 1,2 | ./mdiam -symmetrize -simple -rad -   # radius of an unweighted version which is symmetrized (weak connectivity is assumed)
cat test_digraph.txt | ./mdiam -weighted -closeness-all - 2> /dev/null   # compute closeness out-centrality of all nodes
```

## Usage
```
Usage: ./mdiam [options] [filename]

Read or generate a graph and compute some of its parameters such as radius,
 diameter, or all eccentricities. The graph is considered directed (use
 -symmetrize to obtain a symmetric graph equivalent to an undirected graph) and
 can be weighted (use option -weighted). The graph is either generated from
 options (see, e.g. -grid), or read from a file (use filename - for reading from
 the standard input). The input format consists of a sequence of [u v] pairs
 (for arc u -> v) for unweighted graphs or [u v wgt] triples for weighted graphs
 (space and newline characters are ignored). Nodes can be identified by any
 string (excluding whitespaces).

Distances: this program computes various paramters related to distances in the
 input graph. The distance d(u,v) from u to v is defined as the length of a
 shortest path from u to v (non-negative weights are assumed). The length of a
 path is its number of arcs if the graph is unweighted, and the sum of weigts of
 its arcs if it is weighted. It is considered infinite if v is not reachable
 from u. Note that d(u,v) can be different from d(v,u) in a directed graph.
Efficiency: nodes are numbered as ints as they are encountered, this limits to
 graphs of 2^31 nodes at most. Some computations use multithreading (see
 -n-thread).

Possible options:

  -all-ecc        Compute all eccentricities. The eccentricity e(u) of a node u
 is max_v d(u,v) where d(u,v).

  -closeness-all  Compute the closeness centrality of all nodes. Outputs for
 each node u a quadruple [u r s h] where u is the label of u (as given in
 input), r is the size |R(u)| where R(u) denotes the set of nodes reachable from
 u, s is the sum of distances s = sum_{v in R(u)} d(u,v), and h is the harmonic
 sum h = sum_{v in R(u)} 1 / d(u,v). Closeness centrality is classically defined
 as 1/s, a normalized value can be obtained with (r-1)^2 / (s(n-1)).

  -diameter       Compute the diameter, that is the maximum eccentricity of a
 node (see -all-ecc). Consider -largest-scc for non (strongly) connected
 graphs.

  -grid l         Generate a l x l grid.

  -n-thread       Specify the number of threads to use.

  -power-law b    Generate a random graph according to the configuration model
 such that the degree sequence follows a power law with parameter b.

  -radius         Compute the radius, that is the minimum eccentricity of a node
 (see -all-ecc). Consider -largest-scc for non (strongly) connected graphs.

  -reverse        Reverse all arcs.

  -largest-scc    Restrict the graph to its largest strongly connected
 component.

  -simple         Ensure that the graph is simple by removing duplicate arcs:
 among all u -> v arcs, keep one with minimum weight.

  -symmetrize     Add arc v -> u for each arc u -> v in the original graph.

  -weighted       Read a weighted graph: the sequence read is interpreted as a
 list of triples [u v wgt] instead of pairs [u v].
```


## License

GNU LGPL.

