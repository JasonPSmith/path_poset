# path_poset

Computes the Euler characteristic and homology of the multipath complex studied [here](https://arxiv.org/abs/2208.04656).

Compile in unix using `make`.

Run with `./path_poset n t edges`
where:
* n is number of edges
* t is umber of threads
* edges is the address if a txt file containing the edges of the graph (where each line contains two numbers a b, indicating an edge a->b).

For example: 
`./path_poset 6 8 test.edges`

To compute the homology use:
`./path_poset_homology 6 8 test.edges test.homology`

Note that path_poset_homology requires the boost library.

Alternatively you can use path_poset_sage. Which creates a file that can be fed into [SageMath](https://www.sagemath.org/).

For example:
`./path_poset_sage 6 8 test.edges test.sage`
creates a file test.sage, which when run with
`./sage test.sage`
returns the homology using sagemath

The function path_poset_homology is faster, but only returns Betti's, whilst path_poset_sage is slower, but returns the full homology.
