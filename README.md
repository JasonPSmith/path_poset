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
