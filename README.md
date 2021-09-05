## Compile the project:

g++ --std=c++14 main.cpp VirtualTour.cpp -o VirtualTour

## Run the project:
./VirtualTour < path_to_the_input_file

## Run project on sample input:
./VirtualTour < input.txt

## Input Format:(Look at the sample input in `input.txt`)
_(For details see the input description)_
***
_Input Assumption:_
The graph should be **connected** i.e. any node should be reachable from the start node.
***

1. **No of nodes**(int)
2. **Adjacency Matrix**(square matrix(float) of dimension nodes X nodes):
3. **Amplitude**(maximum velocity can be achieved)(float)
4. **No of components**(No of sine component is fourier series)(int)
5. **Introduce Hypothetical edge**(boolean)
6. **Start Node**(if the tour can start from any node specify -1)(int)
7. **Max allowed time**(if no time restriction specify -1)(float)

## Output Format:
1. **Optimal Time**: The least amount of time required to complete the tour.(Or to
            visit most no of nodes within the specified maximum allowed time.)
2. **Optimal Path**: The sequence of nodes to visit starting from the start node.

## Input Details:
1. No of nodes: No of nodes in the graph.
2. Adjacency Matrix: Adjacency matrix (nodes X nodes matrix) of the graph. If no edge exists
       between a pair of nodes (u, v) then (u, v)th entry of the adjacency matrix
       is denoted by -1 otherwise it is the distance between two nodes.
       adj_matrix[u][v] = -1 if no edge between nodes `u` and `v`
       adj_matrix[u][v] = distance(u, v) if there is an edge between `u` and `v`.
3. Amplitude: Maximum possible velocity (F_m * step_size * visual_gain_factor).
4. No of components: No of sines used in the fourier series.
5. Introduce Hypothetical edge:   If set, introduce the additional edges which
    are not present in original graph. If this parameter is set then input
    graph will be converted into a complete graph. If edge between two nodes u,
    v is not present then we introduce a new edge with length equal to the
    shortest path between u, v.
    E(u, v) = shortest_path(u, v).
6. Start Node: Start node of the tour (-1 if tour can start from any node).
7. Max allowed time:   Max allowed time. The tour should conclude within this
    time. Specify (-1) if there is no time restriction.









