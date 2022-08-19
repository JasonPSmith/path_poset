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

Also contained here is path_poset_sage and path_poset_deltser. These compute homology using [SageMath](https://www.sagemath.org/) or [deltser](https://github.com/JasonPSmith/deltser), respectively. They have an extra input "out_address", and will create a file that can be fed into [SageMath](https://www.sagemath.org/) or [deltser](https://github.com/JasonPSmith/deltser).

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

To compile path_poset_deltser requires the [boost](https://www.boost.org/) library. If you do not wish to compile path_poset_deltser, then remove its name from the first line of Makefile.
