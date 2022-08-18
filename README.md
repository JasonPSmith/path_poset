# path_poset

Computes the Mobiues function (or equivalently Euler Characteristic of the order complex) of the path poset studied [here](https://arxiv.org/abs/2110.11206).

Compile in unix using `make`.

Run with `./path_poset n t edges`
where:
* n is number of edges
* t is umber of threads
* edges is the address if a txt file containing the edges of the graph (where each line contains two numbers a b, indicating an edge a->b).

For example: 
`./path_poset 6 8 test.edges`

Also contained here is path_poset_sage and path_poset_deltser. These have an extra input "out_address", and will create a file that can be feed into sagemath or deltser, respectively, to compute homology.

For example:
`./path_poset_sage 6 8 test.edges test.sage`
creates a file test.sage, which when run with
`./sage test.sage`
returns the homology using sagemath

`./path_poset_deltser 6 8 test.edges test.delt`
creates a file test.delt, which when run with
`./deltser test.sage test.out`
returns the homology using deltser

The function path_poset_deltser is faster, but only returns Betti's, whilst path_poset_sage is slower, but returns the full homology.
