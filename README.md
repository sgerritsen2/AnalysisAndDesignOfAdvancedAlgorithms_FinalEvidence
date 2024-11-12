# 📘 Network Infrastructure Planning and Analysis System

**Authors:**
- Axel Daniel Padilla Reyes - A01642700
- Pablo Esteban Reyes Herrera - A01643307
- Sebastián Gerritsen Ortiz - A01643364

**Creation Date:** November 12, 2024

## Project Description

This project implements a system for planning and analyzing network infrastructure using various graph and computational geometry algorithms. It aims to solve four specific problems often encountered in network design and geographical analysis:

1. **Minimum Spanning Tree (MST)**: Utilizes Kruskal's Algorithm for efficient fiber network planning.
2. **Traveling Salesman Problem (TSP)**: Employs the Nearest Neighbor heuristic for optimizing delivery routes.
3. **Maximum Flow Analysis**: Applies the Ford-Fulkerson Algorithm to analyze network capacity.
4. **Voronoi Diagram Generation**: Creates Voronoi partitions for geographical area segmentation.

## Functionalities

### 1. Minimum Spanning Tree (Fiber Network Planning)
- **Description**: Finds the optimal way to connect nodes (neighborhoods) in a network with minimum cost.
- **Algorithm**: Kruskal's Algorithm.
- **Time Complexity**: \(O(E \log E)\), where \(E\) is the number of edges.
- **Space Complexity**: \(O(V + E)\), where \(V\) is the number of vertices.

### 2. Traveling Salesman Problem (Route Optimization)
- **Description**: Optimizes routes for delivery personnel across all neighborhoods.
- **Algorithm**: Nearest Neighbor heuristic.
- **Time Complexity**: \(O(V^2)\), where \(V\) is the number of vertices.
- **Space Complexity**: \(O(V)\).
- **Note**: This is a heuristic approach and does not guarantee the optimal solution.

### 3. Maximum Flow Analysis (Network Capacity)
- **Description**: Determines the maximum possible data flow between nodes in the network.
- **Algorithm**: Ford-Fulkerson Algorithm with BFS.
- **Time Complexity**: \(O(VE^2)\), where \(V\) is the number of vertices and \(E\) is the number of edges.
- **Space Complexity**: \(O(V^2)\).

### 4. Voronoi Diagram Generation (Geographical Partitioning)
- **Description**: Divides a geographical area into distinct regions around given points, useful for service area delineation.
- **Algorithm**: Delaunay Triangulation (dual of Voronoi diagram).
- **Time Complexity**: \(O(n \log n)\), where \(n\) is the number of points.
- **Space Complexity**: \(O(n)\).

## Test Cases

- **test1.txt**: Validates basic functionality with small, symmetrical matrices.
- **test2.txt**: Assesses scalability of the algorithms on larger datasets.
- **test3.txt**: Tests numerical stability with a variety of input magnitudes.

## Main Functions and Their Complexities

### Kruskal’s Algorithm (Minimum Spanning Tree)
- **Complexity**: \(O(E \log E)\)
- **Description**: Builds the MST by connecting nodes with the least total weight.

### Nearest Neighbor Heuristic (Traveling Salesman Problem)
- **Complexity**: \(O(V^2)\)
- **Description**: Constructs a TSP route by selecting the nearest unvisited neighbor, starting from a chosen node.

### Ford-Fulkerson Algorithm (Maximum Flow)
- **Complexity**: \(O(VE^2)\)
- **Description**: Finds the maximum flow in a flow network by augmenting paths until no more are available.

### Voronoi Diagram (Delaunay Triangulation)
- **Complexity**: \(O(n \log n)\)
- **Description**: Partitions a 2D plane into regions based on proximity to a set of input points.

## Running the Program

1. Compile the program:
   ```bash
   make
Run the executable:
./main < input.txt
input.txt should contain data in the expected format for each section:
MST/TSP adjacency matrix.
Max flow capacity matrix.
Coordinates for Voronoi diagram.
References

GfG. (2023a, November 16). Substring in C. GeeksforGeeks. Link
GfG. (2023b, November 21). Bubble Sort data structure and algorithm tutorials. GeeksforGeeks. Link
GfG. (2024a, January 10). Binary Search Data Structure and algorithm Tutorials. GeeksforGeeks. Link
GfG. (2024b, March 22). QuickSort Data structure and algorithm tutorials. GeeksforGeeks. Link

---
