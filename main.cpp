/*
 * AXEL DANIEL PADILLA REYES - A01642700
 * PABLO ESTEBAN REYES HERRERA - A01643307
 * SEBASTIAN GERRITSEN ORTIZ - A01643364
 * Network Infrastructure Planning and Analysis System
 * This program solves four distinct graph/computational geometry problems:
 * 1. Minimum Spanning Tree (Kruskal's Algorithm) for fiber network planning
 * 2. Traveling Salesman Problem (Nearest Neighbor heuristic) for route optimization
 * 3. Maximum Flow (Ford-Fulkerson Algorithm) for network capacity analysis
 * 4. Voronoi Diagram generation for geographical partitioning
 * test1.txt, Basic functionality
 * test2.txt, Scalability
 * test3.txt, Numerical stability
 * All test cases maintain matrix symmetry (required for undirected graphs),
 * have zeros on the diagonal (no self-loops) and
 * the coordinates for Voronoi diagrams are all positive in each case
 */

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <queue>
#include <limits>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

// Type definitions for CGAL library components
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Delaunay_triangulation_2<K> Delaunay;

/**
 * PART 4: Voronoi Diagram Generation
 * Time Complexity: O(n log n) where n is the number of points
 * Space Complexity: O(n)
 */
Delaunay createVoronoiDiagram(const std::vector<std::pair<int, int>>& points) {
		// Convert input points to CGAL point format
		std::vector<Point_2> cgal_points;
		for (const auto& p : points) {
				cgal_points.push_back(Point_2(p.first, p.second));
		}

		// Create Delaunay triangulation (dual of Voronoi diagram)
		Delaunay dt;
		dt.insert(cgal_points.begin(), cgal_points.end());
		return dt;
}

/**
 * Helper function to print Voronoi diagram vertices
 * Time Complexity: O(v + e) where v is number of vertices and e is number of edges
 */
void printVoronoi(const Delaunay& dt) {
		for (Delaunay::Finite_vertices_iterator v = dt.finite_vertices_begin(); 
				 v != dt.finite_vertices_end(); ++v) {
				std::vector<Point_2> cell_vertices;

				// Get circulator of faces around current vertex
				Delaunay::Face_circulator fc = dt.incident_faces(v);
				Delaunay::Face_circulator done = fc;

				if (fc != nullptr) {
						do {
								if (!dt.is_infinite(fc)) {
										cell_vertices.push_back(dt.dual(fc));
								}
								fc++;
						} while (fc != done);

						// Print vertices in counterclockwise order
						for (size_t i = 0; i < cell_vertices.size(); ++i) {
								std::cout << "(" << cell_vertices[i].x() << ", " << cell_vertices[i].y() << ")";
								if (i < cell_vertices.size() - 1) {
										std::cout << " ";
								}
						}
						std::cout << "\n";
				}
		}
}

using std::cin;
using std::cout;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::queue;
using std::numeric_limits;
using std::min;

// Edge structure for Kruskal's algorithm
struct Edge {
		int src, dest, weight;
};

/**
 * PART 1: Minimum Spanning Tree Helpers
 * Implementation of Union-Find data structure
 * Time Complexity: O(α(n)) where α is the inverse Ackermann function
 */
int find(int u, vector<int>& parent) {
		if (parent[u] != u)
				parent[u] = find(parent[u], parent);
		return parent[u];
}

void unionSet(int u, int v, vector<int>& parent, vector<int>& rank) {
		int rootU = find(u, parent);
		int rootV = find(v, parent);
		if (rank[rootU] < rank[rootV])
				parent[rootU] = rootV;
		else if (rank[rootU] > rank[rootV])
				parent[rootV] = rootU;
		else {
				parent[rootV] = rootU;
				rank[rootU]++;
		}
}

/**
 * PART 1: Kruskal's MST Algorithm
 * Time Complexity: O(E log E) where E is the number of edges
 * Space Complexity: O(V + E) where V is the number of vertices
 */
pair<int, vector<pair<int, int>>> kruskalMST(const vector<vector<int>>& graph) {
		int V = graph.size();
		vector<Edge> edges;

		// Convert adjacency matrix to edge list
		for (int i = 0; i < V; ++i) {
				for (int j = i + 1; j < V; ++j) {
						if (graph[i][j] != 0) {
								edges.push_back({i, j, graph[i][j]});
						}
				}
		}

		// Sort edges by weight
		sort(edges.begin(), edges.end(), [](Edge a, Edge b) {
				return a.weight < b.weight;
		});

		// Initialize Union-Find data structure
		vector<int> parent(V);
		vector<int> rank(V, 0);
		for (int i = 0; i < V; ++i) parent[i] = i;

		int mstWeight = 0;
		vector<pair<int, int>> mstEdges;

		// Build MST using Kruskal's algorithm
		for (const Edge& edge : edges) {
				int u = edge.src;
				int v = edge.dest;

				if (find(u, parent) != find(v, parent)) {
						unionSet(u, v, parent, rank);
						mstWeight += edge.weight;
						mstEdges.push_back({u, v});
						if (mstEdges.size() == V - 1) break;
				}
		}

		return { mstWeight, mstEdges };
}

/**
 * PART 2: Traveling Salesman Problem (Nearest Neighbor)
 * Time Complexity: O(V^2) where V is the number of vertices
 * Space Complexity: O(V)
 * Note: This is a heuristic approach and doesn't guarantee optimal solution
 */
vector<int> tspNearestNeighbor(const vector<vector<int>>& graph, int start) {
		int n = graph.size();
		vector<bool> visited(n, false);
		vector<int> route;

		int current = start;
		route.push_back(current);
		visited[current] = true;

		// Find nearest unvisited neighbor
		for (int i = 1; i < n; i++) {
				int nearest = -1;
				int minDistance = INT_MAX;

				for (int j = 0; j < n; j++) {
						if (!visited[j] && graph[current][j] < minDistance) {
								minDistance = graph[current][j];
								nearest = j;
						}
				}

				route.push_back(nearest);
				visited[nearest] = true;
				current = nearest;
		}

		// Complete the cycle
		route.push_back(start);

		return route;
}

// Helper function to print route using neighborhood labels
void printRoute(const vector<int>& route) {
		for (size_t i = 0; i < route.size(); i++) {
				char neighborhood = 'A' + route[i];
				cout << neighborhood;
				if (i < route.size() - 1) cout << " -> ";
		}
		cout << endl;
}

/**
 * PART 3: Maximum Flow Helpers
 * BFS implementation for Ford-Fulkerson
 * Time Complexity: O(V + E) where V is vertices and E is edges
 */
bool bfs(const vector<vector<int>>& residualGraph, int source, int sink, vector<int>& parent) {
		int n = residualGraph.size();
		vector<bool> visited(n, false);
		queue<int> q;

		q.push(source);
		visited[source] = true;
		parent[source] = -1;

		while (!q.empty()) {
				int u = q.front();
				q.pop();

				for (int v = 0; v < n; v++) {
						if (!visited[v] && residualGraph[u][v] > 0) {
								parent[v] = u;
								visited[v] = true;
								q.push(v);

								if (v == sink) {
										return true;
								}
						}
				}
		}
		return false;
}

/**
 * PART 3: Ford-Fulkerson Algorithm for Maximum Flow
 * Time Complexity: O(VE^2) where V is vertices and E is edges
 * Space Complexity: O(V^2)
 */
int fordFulkerson(const vector<vector<int>>& graph, int source, int sink) {
		int n = graph.size();
		vector<vector<int>> residualGraph = graph;
		vector<int> parent(n);
		int maxFlow = 0;

		// Find augmenting paths and update residual graph
		while (bfs(residualGraph, source, sink, parent)) {
				int pathFlow = numeric_limits<int>::max();

				// Find minimum residual capacity
				for (int v = sink; v != source; v = parent[v]) {
						int u = parent[v];
						pathFlow = min(pathFlow, (int)residualGraph[u][v]);
				}

				// Update residual capacities
				for (int v = sink; v != source; v = parent[v]) {
						int u = parent[v];
						residualGraph[u][v] -= pathFlow;
						residualGraph[v][u] += pathFlow;
				}

				maxFlow += pathFlow;
		}

		return maxFlow;
}

int main() {
		// Read number of nodes
		int n;
		cin >> n;

		// Read MST/TSP graph
		vector<vector<int>> mst_graph(n, vector<int>(n, 0));
		for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
						int v;
						cin >> v;
						mst_graph[i][j] = v;
				}
		}

		// Read max flow graph
		vector<vector<int>> maxflow_graph(n, vector<int>(n, 0));
		for (int i = 0; i < n; i++) {
				for (int j = 0; j < n; j++) {
						int v;
						cin >> v;
						maxflow_graph[i][j] = v;
				}
		}

		// Read coordinates for Voronoi diagram
		vector<pair<int, int>> pairs;
		for (int i = 0; i < n; i++) {
				char ch;
				int a, b;
				cin >> ch >> a >> ch >> b >> ch;
				pairs.emplace_back(a, b);
		}

		// PART 1: Find minimum spanning tree
		cout << "Way of wiring the neighborhoods with fiber:" << endl;
		pair<int, vector<pair<int, int>>> result = kruskalMST(mst_graph);
		vector<pair<int, int>> mstEdges = result.second;
		for (const pair<int, int>& edge : mstEdges) {
				cout << "(" << edge.first << ", " << edge.second << ")" << endl;
		}

		// PART 2: Solve TSP using nearest neighbor
		cout << "Route for mail delivery personnel: " << endl;
		vector<int> route = tspNearestNeighbor(mst_graph, 0);
		printRoute(route);

		// PART 3: Calculate maximum flow
		int maxFlow = fordFulkerson(maxflow_graph, 0, n - 1);
		cout << "Maximum Information Flow: " << endl << maxFlow << endl;

		// PART 4: Generate Voronoi diagram
		cout << "list of polygons:" << endl;
		Delaunay diagram = createVoronoiDiagram(pairs);
		printVoronoi(diagram);

		return 0;
}