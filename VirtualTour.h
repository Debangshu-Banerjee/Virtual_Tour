#ifndef VIRTUALTOUR_H_
#define VIRTUALTOUR_H_

#include <bits/stdc++.h>
#define MAX_DOUBLE 1e9
#define DEBUG_LEVEL 0

using namespace std;
struct VirtualTourParams {
  // No of nodes in the graph.
  int nodes;
  // Maximum possible velocity (F_m * step_size * visual_gain_factor).
  double amplitude;
  // No of sines used in the fourier series.
  int no_of_component;
  // Start node of the tour (-1 if tour can start from any node).
  int start_node;
  // Max allowed time. The tour should conclude within this time. Specify (-1)
  // if there is no time restriction.
  double max_allowed_time;
  // If set, introduce the additional edges which are not present in original
  // graph. If this parameter is set then input graph will be converted into a
  // complete graph. If edge between two nodes u, v is not present then we
  // introduce a new edge with length equal to the shortest path between u, v.
  // E(u, v) = shortest_path(u, v).
  bool introduce_hypothetical_edges;
  // Adjacency matrix (nodes X nodes matrix) of the graph. If no edge exists
  // between a pair of nodes (u, v) then (u, v)th entry of the adjacency matrix
  // is denoted by -1 otherwise it is the distance between two nodes.
  // adj_matrix[u][v] = -1 if no edge between nodes `u` and `v`
  // adj_matrix[u][v] = distance(u, v) if there is an edge between `u` and `v`.
  vector<vector<double>> adjacency_matrix;
  // Weight of the ith node. Input is given as list whose length isequal to no
  // of nodes. (Currently not read as input).
  vector<double> node_weights;
};

class VirtualTour {
 public:
  // Constructor.
  VirtualTour(const VirtualTourParams& params);
  // Computes the virtual tour.
  void ComputeVirtualTour();
  // Get the optimal time required for the virtual tour under the provided
  // constraints. Should be called after ComputeVirtualTour.
  double GetOptimumTime();
  // Get the optimal path(sequence of nodes).
  vector<int> GetOptimumPath();

 private:
  void InitialiseThePathMatrices();
  void ComputeAllPairShortestPaths();
  double ComputeTimeBetweenNodes(double distance);
  void ComputeTimeMatrix();
  void PrintTimeMatrixAndModifiedAdjacencyMatrix();
  void ComputeOptimalPath();
  void FindOptimalPath();
  void SetOptimalCost(int subset_mask, int tour_start_node);
  void SetPathWithHypothticalEdges(int subset_mask, int tour_start_node);
  void SetOriginalPath();
  vector<int> FindShortestPath(int u, int v);

  const int nodes;
  const double amplitude;
  const int no_of_component;
  const bool introduce_hypothetical_edges;
  const int start_node;
  const double max_allowed_time;
  const vector<vector<double>> adjacency_matrix;
  vector<double> node_weights;
  vector<vector<double>> adjacency_matrix_with_hypothetical_edges;
  vector<vector<double>> all_pair_shortest_paths;
  vector<vector<double>> time_matrix;
  double** cost;
  int** parent;
  double* subset_weight;
  vector<int> path_with_hypothtical_edges;
  vector<int> original_path;
  double optimum_cost;
};

#endif  // VIRTUALTOUR_H_
