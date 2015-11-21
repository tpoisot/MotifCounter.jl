# Motif counting in Julia

Step 1 - generate the motif dict

``` sh
julia motiftable.jl
```

Step 2 - run the analysis for the IWDB dataset

``` sh
julia test.jl
```

# Some notes

It's more elegant than the first version! The initial unique vectors are
generated using a binary tree. Then they are folded into matrices, and the
matrices are compared, to identify the unique motifs. This is actually fast
up to 4 nodes, and still feasible for 5 and 6.

The motifs are called `vertices_edges_id`, so `3_3_2` is a motif with three
vertices, three edges, and the unique conformation 2. Look at the `hashes`
object to see which motif this is.

The motif table only needs to be generated once. But after this, the operations
become *fast*. Running a complete motif analysis on the `test.jl` file,
without any parallelism, takes about a minute.
