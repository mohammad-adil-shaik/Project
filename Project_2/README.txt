Structured Overview of the Project:

---

### Project Overview
This project implements three classic graph algorithms:
1. Dijkstra's Shortest Path: It solves the single-source shortest path problem, computing the minimum-distance path beginning with a given source vertex.
2. Kruskal's Minimum Spanning Tree (MST): It builds minimum spanning tree for undirected graphs.
3. DFS for Cycle Detection and Topological Sorting: detect cycle if it exists, allows topological sorting in directed graphs.


Users can load graph data from formatted text files, choose an algorithm to execute, and view the computed results.

---

### Prerequisites
Prior to the project execution, make sure:

1. Java JDK: If not available, it needs to be downloaded and installed.
2. Graph Input Files – Prepare or download graph input files formatted as described below.

---

### Input Structure

Each algorithm reads data from input files with a specific format:

- Dijkstra's Algorithm Input:
  - First line: `numofVertices numofEdges type`
    - `numofVertices` – Total number of vertices in the graph.
    - `numofEdges` – Total number of edges.
    - `type` - 'D' for directed graphs, 'U' for undirected graphs.
    - `Following lines` - Each line describes an edge in the format `<from> <to> <weight>`.

- Kruskal's Algorithm Input:
  - Structure: Similar to Dijkstra's input, but Kruskal’s algorithm processes only undirected graphs.

- DFS Input:
  - Structure: Follows a similar format and uses 'D' for directed graphs, with additional functionality for cycle detection.



### Running the Algorithm

1. Compile the Program:  
   Open a terminal or command prompt, navigate to the project directory, and run:
   ```bash
   javac Main.java
   ```

2. Run the Program:  
   ```bash
   java Main
   ```

3. Choose an Algorithm:  
   When prompted, select:
   - `1` for Dijkstra's Algorithm
   - `2` for Kruskal's Algorithm
   - `3` for DFS - Topological Sorting and Cycle Detection

4. Select an Input File:  
   After choosing an algorithm, select from one of four input file options to load the corresponding graph file for processing.

5. Review Results:  
   - Dijkstra's Algorithm – Displays the shortest path from the source vertex.
   - Kruskal's Algorithm – Shows the edges in the Minimum Spanning Tree and calculates the total cost.
   - DFS – If the graph is cyclic, displays detected cycles; if acyclic, provides the topological order.

6. Execution Time:  
   The program prints the execution time for each algorithm in nanoseconds.

---

### Configuring File Paths

To customize input file paths, adjust the paths in the methods `getInputFilePathDijkstra`, `getInputFilePathKruskal`, and `getInputFilePathDFS` within the code, allowing flexibility to update default paths as needed.

--- 

This overview provides a clear guide for users and developers on how to set up, run, and modify the project as needed.