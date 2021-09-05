#include "VirtualTour.h"

using namespace std;

void ReadInput(VirtualTourParams& params) {
  cin >> params.nodes;
  if (params.nodes <= 0) {
    cout << "No of sites should be positive " << endl;
    exit(0);
  }
  params.adjacency_matrix.resize(params.nodes);
  for (int i = 0; i < params.nodes; i++) {
    params.adjacency_matrix[i].resize(params.nodes);
    for (int j = 0; j < params.nodes; j++) {
      cin >> params.adjacency_matrix[i][j];
    }
  }
  params.amplitude = 1.0;
  cin >> params.amplitude;
  cin >> params.no_of_component;
  cin >> params.introduce_hypothetical_edges;
  cin >> params.start_node;
  cin >> params.max_allowed_time;
}

void PrintInput(const VirtualTourParams& params) {
  cout << "Read Input: " << endl;
  cout << "No of nodes: " << params.nodes << endl;
  cout << "Adjacency Matrix: " << endl;
  for (int i = 0; i < params.nodes; i++) {
    for (int j = 0; j < params.nodes; j++) {
      cout << params.adjacency_matrix[i][j] << " ";
    }
    cout << endl;
  }
  cout << "Amplitude: " << params.amplitude << endl;
  cout << "No of component: " << params.no_of_component << endl;
  cout << "Hypothetical Edges: " << params.introduce_hypothetical_edges << endl;
  cout << "Start Node: " << params.start_node << endl;
  cout << "Maximum Allowed Time: " << params.max_allowed_time << endl;
  cout << endl << endl;
}

int main() {
  VirtualTourParams params;
  ReadInput(params);
  if (DEBUG_LEVEL > -1) {
    PrintInput(params);
  }

  VirtualTour virtual_tour(params);
  // Compute the virtual tour.
  virtual_tour.ComputeVirtualTour();
  // Print output to console.
  cout << "Output: " << endl;
  cout << "Optimal Time: " << virtual_tour.GetOptimumTime() << endl;
  cout << "Optimal Path: " << endl;
  vector<int> optimal_path = virtual_tour.GetOptimumPath();
  for (int i =0; i < optimal_path.size(); i++){
    cout << optimal_path[i] << " ";
  }
  cout << endl;
  return 0;
}
