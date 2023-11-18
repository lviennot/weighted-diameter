# weighted-diameter

(Author : Laurent Viennot, Inria 2023)

Efficient computation of the diameter, radius, and even all eccentricities of a weighted graph in C++.

## Compile

```
make
```

## Examples
```
$ ./mdiam -weighted -largest-scc -diameter graph_example.txt   # diameter of largest strongly connected component
#,diam:17,n:7,m:11,dir:1,wgt:1
17 7 11 1 1 end
#_1_diam 2_n 3_m 4_dir 5_wgt 

$ cat graph_example.txt | cut -d' ' -f 1,2 | ./mdiam -symmetrize -simple -radius -   # radius of an unweighted version which is symmetrized (weak connectivity is assumed)
#,rad:3,n:10,m:29,dir:0,wgt:0
3 10 29 0 0 end
#_1_rad 2_n 3_m 4_dir 5_wgt 

$ cat graph_example.txt | ./mdiam -weighted -closeness-all -no-params -   # compute closeness out-centrality of all nodes
# v |Reach(v)| dist_sum harm_sum
0 10 82 2.43908
4 10 91 1.30466
8 3 5 0.833333
...
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
 graphs of 2^31 nodes at most (the code can easily be patched to use larger
 integers). Some computations use multithreading (see -n-thread).

Possible options:

  -eccentricity-all  Compute all eccentricities. The eccentricity e(u) of a node
 u is max_v d(u,v) where d(u,v).

  -closeness-all  Compute the closeness centrality of all nodes. Outputs for
 each node u a quadruple [u r s h] where u is the label of u (as given in
 input), r is the size |R(u)| where R(u) denotes the set of nodes reachable from
 u, s is the sum of distances s = sum_{v in R(u)} d(u,v), and h is the harmonic
 sum h = sum_{v in R(u)} 1 / d(u,v). Closeness centrality is classically defined
 as 1/s, a normalized value can be obtained with (r-1)^2 / (s(n-1)). Use
 -columns-verb 0 to output only those lines.

  -params-verb v   Computed parameters are output to stdout on one line,
 separated by spaces. The first four fields are [n m dir wgt]: the number n of
 nodes, the number m of edges, dir=1 for a directed graph (0 if -symmetrize is
 used) and wgt=1 for a weighted graph (use -weighted). Use -diameter, -radius
 for diameter, and/or radius. Nothing is output with v=0, the default is 1,
 while more information about the computation is output with 2.

  -diameter       Compute the diameter, that is the maximum eccentricity of a
 node (see -eccentricity-all). Consider -largest-scc for non (strongly)
 connected graphs.

  -grid l         Generate a l x l grid.

  -n-thread       Specify the number of threads to use. Options using
 multithreading are -closeness-all, 

  -power-law b    Generate a random graph according to the configuration model
 such that the degree sequence follows a power law with parameter b.

  -radius         Compute the radius, that is the minimum eccentricity of a node
 (see -eccentricity-all). Consider -largest-scc for non (strongly) connected
 graphs.

  -reverse        Reverse all arcs.

  -largest-scc    Restrict the graph to its largest strongly connected
 component.

  -simple         Ensure that the graph is simple by removing duplicate arcs:
 among all u -> v arcs, keep one with minimum weight.

  -symmetrize     Add arc v -> u for each arc u -> v in the original graph.

  -verbosity v    Print on stderr various information depending on v=0,1,2
 (defaults to 0 for nothing).

  -weighted       Read a weighted graph: the sequence read is interpreted as a
 list of triples [u v wgt] instead of pairs [u v].
```


## License

GNU LGPL.

