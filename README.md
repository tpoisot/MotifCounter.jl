# Motif counting in Julia

Step 1 - generate the motif dict

```sh
julia motiftable.jl
```

Step 2 - run the analysis for the IWDB dataset

```sh
julia test.jl
```

# Some notes

It's not that elegant. The first step generates a giant hash table (which really
is the transformation of the matrix into a string). Then motif counting is just
an aggregation on this.

It is currently limited to motifs with 2 and 3 nodes because of the sheer number
of combinations for 4 nodes and more. This will have to wait until I find a way
to distribute the operations properly.
