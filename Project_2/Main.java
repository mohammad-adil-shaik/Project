import java.io.File;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Scanner;
import java.util.Set;
import java.io.File;
import java.io.FileReader;
import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;

public class Main {

    // Applying Dijkstra’s algorithm for pathfinding
    static class Dijkstra_Algorithm {
        FileReader reader = null;
        Scanner scanner = null;
        int numofVertices, numofEdges, originVertex;
        int[] distances;
        int[] previousVertices;
        String vertexLabels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        String currentLine = null;
        boolean isDirectedGraph = false;
        List<List<VertexNode>> adjacencyList;
        int edgeFrom;
        int edgeTo;
        int edgeWeight;

        public void fileReader(String filePath) throws Exception {
            reader = new FileReader(new File(filePath));
            if (reader != null) {
                scanner = new Scanner(reader);

                if ((currentLine = scanner.nextLine()) != null) {
                    String[] tokens = currentLine.split(" ");

                    if (tokens[2].equals("D")) {
                        isDirectedGraph = true;
                        System.out.println("Type of Graph: Directed");
                    } else {
                        System.out.println("Type of Graph: Undirected");
                    }

                    numofVertices = Integer.parseInt(tokens[0]);
                    System.out.println("Total Vertices Count: " + numofVertices);

                    numofEdges = Integer.parseInt(tokens[1]);
                    System.out.println("Total Edges Count: " + numofEdges);

                    // Neighboring vertices list (Adjacency list)
                    adjacencyList = new ArrayList<List<VertexNode>>();
                    // Setting up adjacency list for each vertex
                    for (int index = 0; index < numofVertices; index++) {
                        adjacencyList.add(new ArrayList<VertexNode>());
                    }

                    for (int index = 0; index < numofEdges; index++) {
                        currentLine = scanner.nextLine();
                        tokens = currentLine.split(" ");

                        edgeFrom = vertexLabels.indexOf(tokens[0].charAt(0));
                        edgeTo = vertexLabels.indexOf(tokens[1].charAt(0));
                        edgeWeight = Integer.parseInt(tokens[2]);

                        adjacencyList.get(edgeFrom).add(new VertexNode(edgeTo, edgeWeight));
                        // In case the graph is undirected
                        if (!isDirectedGraph) {
                            adjacencyList.get(edgeTo).add(new VertexNode(edgeFrom, edgeWeight));
                        }
                    }

                    // Condition for source presence in input data
                    if (scanner.hasNext()) {
                        currentLine = scanner.nextLine();
                        originVertex = vertexLabels.indexOf(currentLine.charAt(0));
                        System.out.println("\nProvided source vertex:  " + currentLine.charAt(0) + "\n");
                    } else {
                        originVertex = 0;
                        System.out.println(
                                "\nSource vertex not specified. Defaulting to \n" + vertexLabels.charAt(0) + " as source vertex.\n");
                    }
                }
            }
        }

        void dijkstra() {
            distances = new int[numofVertices];
            previousVertices = new int[numofVertices];

            for (int index = 0; index < numofVertices; index++) {
                distances[index] = Integer.MAX_VALUE;
            }

            distances[originVertex] = 0;
            previousVertices[originVertex] = -1;

            Set<Integer> visitedSet = new HashSet<Integer>();

            PriorityQueue<VertexNode> priorityQueue = new PriorityQueue<VertexNode>(numofVertices, new VertexNode());

            priorityQueue.add(new VertexNode(originVertex, 0));

            while (!priorityQueue.isEmpty()) {
                int currentVertex = priorityQueue.remove().vertex;
                visitedSet.add(currentVertex);
                for (int index = 0; index < adjacencyList.get(currentVertex).size(); index++) {
                    int adjacentVertex = adjacencyList.get(currentVertex).get(index).vertex;
                    if (!visitedSet.contains(adjacentVertex)) {
                        int edgeWeight = adjacencyList.get(currentVertex).get(index).cost; // Edge cost from currentVertex to adjacentVertex
                        if (distances[currentVertex] + edgeWeight < distances[adjacentVertex]) {
                            distances[adjacentVertex] = distances[currentVertex] + edgeWeight;
                            previousVertices[adjacentVertex] = currentVertex;
                            priorityQueue.add(new VertexNode(adjacentVertex, distances[adjacentVertex]));
                        }
                    }
                }
            }
        }

        public void printPath() {
            System.out.print("Calculating shortest path from the source vertex to all other vertices:\n");
            for (int index = 0; index < numofVertices; index++) {
                if (index != originVertex) {
                    int vertex = index;
                    String path = "";
                    while (previousVertices[vertex] != -1) {
                        path = " -> " + vertexLabels.charAt(vertex) + path;
                        vertex = previousVertices[vertex];
                    }
                    System.out.print("\n" + vertexLabels.charAt(originVertex) + path + " Path cost: ");
                    System.out.print(distances[index]);
                }
            }
        }

        public class VertexNode implements Comparator<VertexNode> {
            public int vertex;
            public int cost;

            public VertexNode() {
            }

            public VertexNode(int vertex, int cost) {
                this.vertex = vertex;
                this.cost = cost;
            }

            public int compare(VertexNode n1, VertexNode n2) {
                if (n1.cost < n2.cost)
                    return -1;
                if (n1.cost > n2.cost)
                    return 1;
                return 0;
            }
        }
    }

    // Applying Kruskal’s algorithm for MST
    static class Kruskal_Algorithm {

        FileReader reader = null;
        Scanner scanner = null;
        int mstEdgeIndex;
        int numofVertices, numofEdges;
        Edge mstEdges[];
        String vertexLabels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        Edge allEdges[];

        public void fileReader(String filePath) throws Exception {
            reader = new FileReader(new File(filePath));
            scanner = new Scanner(reader);
            numofVertices = scanner.nextInt();
            numofEdges = scanner.nextInt();
            System.out.println("Total Vertices Count: " + numofVertices);
            System.out.println("Total Edges Count: " + numofEdges);

            scanner.next().charAt(0); // Omit this character
            allEdges = new Edge[numofEdges];
            for (int i = 0; i < numofEdges; ++i) {
                allEdges[i] = new Edge();
            }

            for (int index = 0; index < numofEdges; index++) {
                int sourceVertex = vertexLabels.indexOf(scanner.next().charAt(0));
                int destinationVertex = vertexLabels.indexOf(scanner.next().charAt(0));
                int edgeWeight = scanner.nextInt();
                allEdges[index].src = sourceVertex;
                allEdges[index].dest = destinationVertex;
                allEdges[index].weight = edgeWeight;
            }
        }

        // Locate the set for the given vertex
        int findSet(Subset subsets[], int vertex) {
            if (subsets[vertex].parent != vertex)
                subsets[vertex].parent = findSet(subsets, subsets[vertex].parent);
            return subsets[vertex].parent;
        }

        // Join two sets by rank to maintain tree balance
        void union(Subset subsets[], int x, int y) {
            int xroot = findSet(subsets, x);
            int yroot = findSet(subsets, y);
            if (subsets[xroot].rank < subsets[yroot].rank)
                subsets[xroot].parent = yroot;
            else if (subsets[xroot].rank > subsets[yroot].rank)
                subsets[yroot].parent = xroot;
            else {
                subsets[yroot].parent = xroot;
                subsets[xroot].rank++;
            }
        }

        void kruskal() {
            mstEdges = new Edge[numofVertices];
            mstEdgeIndex = 0;
            int index;
            for (index = 0; index < numofVertices; ++index) {
                mstEdges[index] = new Edge();
            }

            // Arrange edges by weight in non-decreasing order
            Arrays.sort(allEdges);

            Subset subsets[] = new Subset[numofVertices];

            // Set up an individual set for each vertex
            for (index = 0; index < numofVertices; ++index) {
                subsets[index] = new Subset();
            }

            // Set up hierarchy with parent and rank for subsets
            for (int v = 0; v < numofVertices; ++v) {
                subsets[v].parent = v;
                subsets[v].rank = 0;
            }

            index = 0;

            // Continue looping until MST reaches numVertices - 1 edges
            while (mstEdgeIndex < numofVertices - 1) {
                Edge nextEdge = allEdges[index++];
                int x = findSet(subsets, nextEdge.src);
                int y = findSet(subsets, nextEdge.dest);
                if (x != y) {
                    mstEdges[mstEdgeIndex++] = nextEdge;
                    union(subsets, x, y);
                }
            }
        }

        public void display() {
            int totalCost = 0;
            System.out.println("\nnDisplaying Minimum Spanning Tree: \n");
            for (int i = 0; i < mstEdgeIndex; ++i) {
                System.out.println(vertexLabels.charAt(mstEdges[i].src) + " -> " +
                        vertexLabels.charAt(mstEdges[i].dest) + " Cost: " + mstEdges[i].weight);
                totalCost += mstEdges[i].weight;
            }

            System.out.println("\nTotal cost of the Minimum Spanning Tree (MST): " + totalCost);
        }

        // Defines a graph edge as a class
        class Edge implements Comparable<Edge> {
            int src, dest, weight;

            // Organize edges based on weight
            public int compareTo(Edge compareEdge) {
                return this.weight - compareEdge.weight;
            }
        }

        class Subset {
            int parent, rank;
        }
    }

    // Applying DFS for topological order generation
    static class DFS_Algorithm {
        FileReader reader = null;
        Scanner scanner = null;
        int numofVertices, numofEdges;
        boolean isDirectedGraph = false;
        String vertexLabels = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
        List<List<Integer>> adjacencyList;
        Set<String> cycles = new HashSet<>();
        List<String> topoSortStack = new ArrayList<>();
        boolean isCyclic = false;

        final int WHITE = 0;
        final int GRAY = 1;
        final int BLACK = 2;

        public void fileReader(String filePath) throws Exception {
            reader = new FileReader(new File(filePath));
            if (reader != null) {
                scanner = new Scanner(reader);

                String currentLine = null;

                if (scanner.hasNextLine()) {
                    currentLine = scanner.nextLine();
                    String[] tokens = currentLine.trim().split("\\s+");

                    numofVertices = Integer.parseInt(tokens[0]);
                    numofEdges = Integer.parseInt(tokens[1]);
                    String graphTypeStr = tokens[2];

                    if (graphTypeStr.equals("D")) {
                        isDirectedGraph = true;
                        System.out.println("Type of Graph: Directed");
                    } else {
                        System.out.println("Type of Graph: Undirected");
                    }

                    System.out.println("Total Vertices Count: " + numofVertices);
                    System.out.println("Total Edges Count: " + numofEdges);

                    adjacencyList = new ArrayList<>(numofVertices);
                    for (int i = 0; i < numofVertices; i++) {
                        adjacencyList.add(new ArrayList<>());
                    }

                    for (int i = 0; i < numofEdges; i++) {
                        if (scanner.hasNextLine()) {
                            currentLine = scanner.nextLine();
                            tokens = currentLine.trim().split("\\s+");

                            int u = vertexLabels.indexOf(tokens[0]);
                            int v = vertexLabels.indexOf(tokens[1]);

                            adjacencyList.get(u).add(v);
                            if (!isDirectedGraph) {
                                adjacencyList.get(v).add(u);
                            }
                        } else {
                            throw new Exception("Insufficient edges provided in the input file");
                        }
                    }

                    // Source node is optional; not utilized in current DFS
                    if (scanner.hasNextLine()) {
                        currentLine = scanner.nextLine().trim();
                        if (!currentLine.isEmpty()) {
                           // This DFS implementation does not utilize an optional source node
                        }
                    }
                }
            }
        }

        public void dfs() {
            int[] color = new int[numofVertices];
            int[] parent = new int[numofVertices];
            for (int i = 0; i < numofVertices; i++) {
                color[i] = WHITE;
                parent[i] = -1;
            }
            for (int u = 0; u < numofVertices; u++) {
                if (color[u] == WHITE) {
                    dfsVisit(u, color, parent);
                }
            }
        }

        private void dfsVisit(int u, int[] color, int[] parent) {
            color[u] = GRAY;
            for (int v : adjacencyList.get(u)) {
                if (color[v] == WHITE) {
                    parent[v] = u;
                    dfsVisit(v, color, parent);
                } else if (color[v] == GRAY) {
                    // Cycle identified through back edge
                    if (isDirectedGraph || parent[u] != v) {
                        isCyclic = true;
                        storeCycle(u, v, parent);
                    }
                }
            }
            color[u] = BLACK;
            topoSortStack.add(vertexLabels.substring(u, u + 1)); // Store the vertex's label
        }

        private void storeCycle(int start, int end, int[] parent) {
            List<Integer> cycleVertices = new ArrayList<>();
            cycleVertices.add(end);
            int current = start;
            while (current != end && current != -1) {
                cycleVertices.add(current);
                current = parent[current];
            }
            cycleVertices.add(end); // to complete the cycle
            Collections.reverse(cycleVertices);
            StringBuilder cycleStr = new StringBuilder();
            for (int v : cycleVertices) {
                cycleStr.append(vertexLabels.charAt(v)).append(" -> ");
            }
            cycleStr.delete(cycleStr.length() - 4, cycleStr.length()); // remove the last ' -> '
            if (!cycles.contains(cycleStr.toString())) {
                cycles.add(cycleStr.toString());
            }
        }

        public void printCycles() {
            System.out.println("Cycles detected in the graph:");
            int idx = 1;
            for (String cycle : cycles) {
                String[] vertices = cycle.split(" -> ");
                System.out.println("Cycle " + idx + " (length " + (vertices.length - 1) + "): " + cycle);
                idx++;
            }
        }

        public void printTopologicalSort() {
            System.out.println("Order of Topological Sorting:");
            Collections.reverse(topoSortStack);
            System.out.println(String.join(" -> ", topoSortStack));
        }
    }

	public static void main(String[] args) throws Exception {
		
		long executionStartTime, executionEndTime;
		Dijkstra_Algorithm dijkstra = new Dijkstra_Algorithm();
		Kruskal_Algorithm kruskal = new Kruskal_Algorithm();
        DFS_Algorithm dfs = new DFS_Algorithm();
		
		Scanner scan = new Scanner(System.in);
        String filePath = "1"; //default input file 
        int inputFile;
        
    
        System.out.println("Please choose an option (1-3) to execute the desired algorithm:\n" +
                   "1: Dijkstra's Algorithm\n" +
                   "2: Kruskal's Algorithm\n" +
                   "3: DFS - Topological Sorting and Cycle Detection");

        int option = Integer.parseInt(scan.nextLine());
        
        switch (option) {
		case 1:
        System.out.println("Please choose an option (1-4) to select your desired input file:\n" +
        "1: Input File 1 (Undirected graph)\n" +
        "2: Input File 2 (Directed graph)\n" +
        "3: Input File 3 (Undirected graph)\n" +
        "4: Input File 4 (Undirected graph)");
			inputFile = Integer.parseInt(scan.nextLine());
			filePath = getInputFilePathDijkstra(inputFile);
			
			// Applying Dijkstra's Algorithm to directed or undirected graphs
			System.out.println("DJIKSTRA'S ALGORITHM\n");
			dijkstra.fileReader(filePath);
			executionStartTime = System.nanoTime();
			dijkstra.dijkstra();
			executionEndTime = System.nanoTime();
			dijkstra.printPath();
			System.out.println("\n\nExecution time of Dijkstra's Algorithm (in nanoseconds):"+(executionEndTime-executionStartTime));
			break;
		case 2:
        System.out.println("Please select an option (1-4) to choose your preferred input file:\n" +
        "1: Input File 1 (Undirected graph).\n" +
        "2: Input File 2 (Undirected graph).\n" +
        "3: Input File 3 (Undirected graph).\n" +
        "4: Input File 4 (Undirected graph).");
			inputFile = Integer.parseInt(scan.nextLine());
			filePath = getInputFilePathKruskal(inputFile);
			
			// Applying Kruskal's Algorithm to undirected graphs
			System.out.println("KRUSKAL'S ALGORITHM\n");
			kruskal.fileReader(filePath);
			executionStartTime = System.nanoTime();
			kruskal.kruskal();
			executionEndTime = System.nanoTime();
			kruskal.display();
			System.out.println("\nExecution time of Kruskal's Algorithm (in nanoseconds): "+(executionEndTime-executionStartTime));
			break; 
        case 3:
        System.out.println("Please select an option (1-4) to choose your preferred input file:\n" +
                           "1: Input File 1 (Directed, acyclic graph).\n" +
                           "2: Input File 2 (Directed, acyclic graph).\n" +
                           "3: Input File 3 (Directed, cyclic graph).\n" +
                           "4: Input File 4 (Directed, cyclic graph).");
        
            inputFile = Integer.parseInt(scan.nextLine());
            filePath = getInputFilePathDFS(inputFile);

            // Applying DFS to detect cycles and achieve topological sort
            System.out.println("DFS ALGORITHM\n");
            dfs.fileReader(filePath);
            executionStartTime = System.nanoTime();
            dfs.dfs();
            executionEndTime = System.nanoTime();
            if (dfs.isCyclic) {
                dfs.printCycles();
            } else {
                dfs.printTopologicalSort();
            }
            System.out.println("\nDFS Algorithm execution time (in nanoseconds): " + (executionEndTime - executionStartTime));
            break;
        default:
            System.out.println("Invalid input! Please try again!");
            break;
		}
	}
	
	static String getInputFilePathDijkstra(int inputFile) {
		String filePath;
		switch (inputFile) {
		case 1:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dijkstra_graphs/Undirected_1.txt";
			break;
		case 2:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dijkstra_graphs/directed_2.txt";
			break;
		case 3:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dijkstra_graphs/Undirected_3.txt";
			break;
		case 4:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dijkstra_graphs/Undirected_4.txt";
			break;
		default:
			System.out.println("You selected an invalid option. Default file selected is 1.\n");
			filePath = "/Users/gnanaprakashnarayanairivisetty/AlgoDS_Project-2/dijkstra_graphs/Undirected_graph_1.txt";
			break;
		}
		return filePath;
	}
	
	static String getInputFilePathKruskal(int inputFile) {
		String filePath;
		switch (inputFile) {
		case 1:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/kruskal_graphs/Undirected_1.txt";
			break;
		case 2:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/kruskal_graphs/Undirected_2.txt";
			break;
		case 3:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/kruskal_graphs/Undirected_3.txt";
			break;
		case 4:
			filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/kruskal_graphs/Undirected_4.txt";
			break;
		default:
			System.out.println("You selected an invalid option. Default file selected is 1.\n");
			filePath = "/Users/gnanaprakashnarayanairivisetty/AlgoDS_Project-2/kruskal_graphs/Undirected_graph_1.txt";
			break;
		}
		return filePath;
	}
    static String getInputFilePathDFS(int inputFile) {
        String filePath;
        switch (inputFile) {
            case 1:
                filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dfs_graphs/directed_1_acyclic.txt";
                break;
            case 2:
                filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dfs_graphs/directed_2_acyclic.txt";
                break;
            case 3:
                filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dfs_graphs/directed_3_cyclic.txt";
                break;
            case 4:
                filePath = "C:/Users/ADIL DRAGSTER/OneDrive/Desktop/Assgmts/Algo DS/Adil/Project_2/dfs_graphs/directed_4_cyclic.txt";
                break;
            default:
                System.out.println("You selected an invalid option. Default file selected is 1.\n");
                filePath = "/Users/gnanaprakashnarayanairivisetty/AlgoDS_Project-2/dfs_graphs/directed_graph_1(acyclic).txt";
                break;
        }
        return filePath;
    }
}
