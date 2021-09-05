#include "VirtualTour.h"

double PI = 3.1415;
using namespace std;

VirtualTour::VirtualTour(const VirtualTourParams& params)
    : nodes(params.nodes),
      amplitude(params.amplitude),
      no_of_component(params.no_of_component),
      introduce_hypothetical_edges(params.introduce_hypothetical_edges),
      start_node(params.start_node),
      max_allowed_time(params.max_allowed_time),
      adjacency_matrix(params.adjacency_matrix) {
  // Check the size of the adjacency matrix.
  if (adjacency_matrix.size() != nodes) {
    cout << "Adjacency matrix does not have " << nodes << "rows\n";
    exit(0);
  }
  if (adjacency_matrix[0].size() != nodes) {
    cout << "Adjacency matrix does not have " << nodes << "columns\n";
    exit(0);
  }

  this->adjacency_matrix_with_hypothetical_edges.resize(nodes);
  this->all_pair_shortest_paths.resize(nodes);
  this->time_matrix.resize(nodes);
  if (params.node_weights.size() == nodes) {
    this->node_weights = params.node_weights;
  } else {
    this->node_weights.clear();
    this->node_weights.resize(nodes, 1.0);
  }
  for (int i = 0; i < nodes; i++) {
    this->adjacency_matrix_with_hypothetical_edges[i].resize(nodes);
    this->all_pair_shortest_paths[i].resize(nodes);
    this->time_matrix[i].resize(nodes);
  }
  this->cost = NULL;
  this->parent = NULL;
}

void VirtualTour::InitialiseThePathMatrices() {
  for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      if (i == j) {
        adjacency_matrix_with_hypothetical_edges[i][j] = 0;
        all_pair_shortest_paths[i][j] = -1;
        continue;
      }
      adjacency_matrix_with_hypothetical_edges[i][j] =
          adjacency_matrix[i][j] < 0 ? MAX_DOUBLE : adjacency_matrix[i][j];
      all_pair_shortest_paths[i][j] = adjacency_matrix[i][j] < 0 ? -1 : i;
    }
  }
}

void VirtualTour::ComputeAllPairShortestPaths() {
  // Initialise the matrices
  InitialiseThePathMatrices();
  // Run Floyd Warshall algorithm to compute all pair shortest paths.
  for (int k = 0; k < nodes; k++) {
    for (int v = 0; v < nodes; v++) {
      for (int u = 0; u < nodes; u++) {
        // If vertex `k` is on the shortest path from `v` to `u`,
        // then update the value of --->
        // 1. adjacency_matrix_with_hypothetical_edges[v][u]
        // 2. all_pair_shortest_paths[v][u]

        if (adjacency_matrix_with_hypothetical_edges[v][k] != MAX_DOUBLE &&
            adjacency_matrix_with_hypothetical_edges[k][u] != MAX_DOUBLE &&
            (adjacency_matrix_with_hypothetical_edges[v][k] +
             adjacency_matrix_with_hypothetical_edges[k][u]) <
                adjacency_matrix_with_hypothetical_edges[v][u]) {
          adjacency_matrix_with_hypothetical_edges[v][u] =
              adjacency_matrix_with_hypothetical_edges[v][k] +
              adjacency_matrix_with_hypothetical_edges[k][u];
          all_pair_shortest_paths[v][u] = all_pair_shortest_paths[k][u];
        }
      }

      // if diagonal elements become negative, the
      // graph contains a negative-weight cycle
      if (adjacency_matrix_with_hypothetical_edges[v][v] < 0) {
        cout
            << "Negative-weight cycle found!!. Check the adjacency matrix input"
            << endl;
        exit(0);
      }
    }
  }
}

double VirtualTour::GetOptimumTime() {
 return optimum_cost;
}

vector<int> VirtualTour::GetOptimumPath() {
  if(introduce_hypothetical_edges) {
    return path_with_hypothtical_edges;
  }
  return original_path;
}

double VirtualTour::ComputeTimeBetweenNodes(double distance) {
  // Formula for calculating time between the nodes.
  // t(d) = 1/amplitude * distance * (2 - 8/pi^2(1 + 1/3^ + 1/5^2 ... ))
  // The no of terms of this infinite series to use is specified by the
  // parameter `no_of_component`.
  double squre_series_value = 0;
  for (int i = 0; i < no_of_component; i++) {
    squre_series_value += 1.0 / (double)((i + 1.0) * (i + 1.0));
  }
  squre_series_value *= 8.0 / (PI * PI);
  double time_between_nodes = distance / amplitude * (2 - squre_series_value);
  return time_between_nodes;
}

void VirtualTour::ComputeTimeMatrix() {
  for (int i = 0; i < nodes; i++) {
    for (int j = 0; j < nodes; j++) {
      if (adjacency_matrix_with_hypothetical_edges[i][j] < 0 ||
          adjacency_matrix_with_hypothetical_edges[i][j] >= MAX_DOUBLE) {
        cout << "Check input some node pairs have infinite distance after "
                "hypothetical edges were introduced\n";
        exit(0);
      }
      time_matrix[i][j] = ComputeTimeBetweenNodes(
          adjacency_matrix_with_hypothetical_edges[i][j]);
    }
  }
}

void VirtualTour::ComputeOptimalPath() {
  int bit_mask_len = (1 << nodes) - 1;
  if (nodes <= 1) {
    cout << "One or less nodes present. Tour not possible.\n";
    exit(0);
  }
  cost = (double**)malloc((bit_mask_len + 1) * sizeof(double*));
  parent = (int**)malloc((bit_mask_len + 1) * sizeof(int*));
  subset_weight = (double*)malloc((bit_mask_len + 1) * sizeof(double));
  for (int i = 0; i < bit_mask_len + 1; i++) {
    cost[i] = (double*)malloc(nodes * sizeof(double));
    parent[i] = (int*)malloc(nodes * sizeof(int));
    subset_weight[i] = -1.0;
    for (int j = 0; j < nodes; j++) {
      cost[i][j] = MAX_DOUBLE;
      parent[i][j] = -1;
    }
  }
  for (int i = 0; i < nodes; i++) {
    cost[1 << i][i] = 0.0;
  }
  subset_weight[0] = 0.0;
  for (int i = 0; i < bit_mask_len + 1; i++) {
    for (int j = 0; j < nodes; j++) {
      if (i & (1 << j)) {
        if (subset_weight[i] < 0) {
          subset_weight[i] = subset_weight[i ^ (1 << j)] + node_weights[j];
        }
        for (int k = 0; k < nodes; k++) {
          if (k == j) continue;
          if ((i & (1 << k)) && cost[i ^ (1 << j)][k] < MAX_DOUBLE) {
            double total_cost = cost[i ^ (1 << j)][k] + time_matrix[j][k];
            if (total_cost < cost[i][j]) {
              cost[i][j] = total_cost;
              parent[i][j] = k;
            }
          }
        }
      }
      if (DEBUG_LEVEL > 1) {
        cout << i << " " << j << " " << cost[i][j] << " " << parent[i][j]
             << endl;
      }
    }
  }
}

// Returns the path without the start node `u`.
vector<int> VirtualTour::FindShortestPath(int u, int v) {
  vector<int> path;
  if (u < 0 || v < 0 || u >= nodes || v >= nodes) return path;
  path.push_back(v);
  int current_node = all_pair_shortest_paths[u][v];
  while (current_node != u) {
    path.push_back(current_node);
    if (current_node < 0) break;
    current_node = all_pair_shortest_paths[u][current_node];
  }
  reverse(path.begin(), path.end());
  return path;
}

void VirtualTour::SetOptimalCost(int subset_mask, int tour_start_node) {
  optimum_cost = cost[subset_mask][tour_start_node];
  if (DEBUG_LEVEL > 0) {
    cout << "Optimal time: " << optimum_cost << endl;
  }
}

void VirtualTour::SetOriginalPath() {
  if (path_with_hypothtical_edges.empty()) return;
  original_path.push_back(path_with_hypothtical_edges[0]);
  for (int i = 1; i < path_with_hypothtical_edges.size(); i++) {
    vector<int> path_extension =
        FindShortestPath(original_path.back(), path_with_hypothtical_edges[i]);
    original_path.insert(original_path.end(), path_extension.begin(),
                         path_extension.end());
  }
  if (DEBUG_LEVEL > 0) {
    cout << "Original path\n";
    for (int node : original_path) {
      cout << node << " ";
    }
    cout << endl;
  }
}

void VirtualTour::SetPathWithHypothticalEdges(int subset_mask,
                                              int tour_start_node) {
  if (subset_mask < 0 || subset_mask >= (1 << nodes) || tour_start_node < 0 ||
      tour_start_node >= nodes) {
    cout << "Error in calculating the optimal path\n";
    exit(0);
  }
  // Set the optimal cost of the the virtual tour.
  SetOptimalCost(subset_mask, tour_start_node);
  while (tour_start_node != -1) {
    path_with_hypothtical_edges.push_back(tour_start_node);
    int current_node = tour_start_node;
    tour_start_node = parent[subset_mask][tour_start_node];
    subset_mask = subset_mask ^ (1 << current_node);
  }
  if (DEBUG_LEVEL > 0) {
    cout << "Optimal path node: " << endl;
    for (int x : path_with_hypothtical_edges) {
      cout << x << " ";
    }
    cout << endl;
  }
  SetOriginalPath();
}

void VirtualTour::FindOptimalPath() {
  int subset_mask;
  int tour_start_node = start_node < 0 ? -1 : start_node;
  if (max_allowed_time < 0) {
    subset_mask = (1 << nodes) - 1;
    if (start_node < 0) {
      double min_time = cost[subset_mask][0];
      tour_start_node = 0;
      for (int i = 0; i < nodes; i++) {
        if (cost[subset_mask][i] < min_time) {
          min_time = cost[subset_mask][i];
          tour_start_node = i;
        }
      }
    }
  } else {
    double total_node_weight = 0.0;
    double min_time = MAX_DOUBLE;
    subset_mask = 0;
    int total_size = (1 << nodes);
    for (int i = 0; i < total_size; i++) {
      if (subset_weight[i] < total_node_weight) continue;
      if (start_node < 0) {
        for (int j = 0; j < nodes; j++) {
          if (cost[i][j] > max_allowed_time) continue;
          if (subset_weight[i] > total_node_weight || (cost[i][j] < min_time)) {
            subset_mask = i;
            total_node_weight = subset_weight[i];
            tour_start_node = j;
            min_time = cost[i][j];
          }
        }
      } else {
        if (cost[i][start_node] <= max_allowed_time) {
          if (subset_weight[i] > total_node_weight ||
              (cost[i][start_node] < min_time)) {
            subset_mask = i;
            total_node_weight = subset_weight[i];
            min_time = cost[i][start_node];
          }
        }
      }
    }
  }
  if (DEBUG_LEVEL > 0) {
    cout << subset_mask << " " << tour_start_node << endl;
  }
  SetPathWithHypothticalEdges(subset_mask, tour_start_node);
}

void VirtualTour::ComputeVirtualTour() {
  // Compute the adjacency matrix with hypothetical edges.
  ComputeAllPairShortestPaths();
  // Once the complete adjacency matrix is computed, compute the time matrix
  // which will have the time taken to travel between any pair of node
  // directly or through other nodes.
  ComputeTimeMatrix();
  // For debugging print out the computed the time and modified distance
  // matrix.
  if (DEBUG_LEVEL > 0) {
    cout << "Printing the adjacency matrix " << endl;
    PrintTimeMatrixAndModifiedAdjacencyMatrix();
  }
  ComputeOptimalPath();
  FindOptimalPath();
}

void PrintMatrix(const string& matrix_name,
                 const vector<vector<double>>& matrix) {
  cout << matrix_name << ": " << endl;
  for (int i = 0; i < matrix.size(); i++) {
    for (int j = 0; j < matrix[0].size(); j++) {
      cout << matrix[i][j] << " ";
    }
    cout << endl;
  }
}

void VirtualTour::PrintTimeMatrixAndModifiedAdjacencyMatrix() {
  PrintMatrix("Modified Adjacency Matrix",
              adjacency_matrix_with_hypothetical_edges);
  PrintMatrix("Time Matrix", time_matrix);
}
