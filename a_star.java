import java.util.*;
import java.nio.file.Path;
import java.util.Random;


public class a_star {

    
    public static class Edge {
            private Map<Integer, Integer> nodes = new HashMap<>();
            private int weight;
            private  int edgeId;
            private List<Integer> channels;

            public Edge(int node1, int node2, int weight, int edgeId, int P) {
                this.nodes.put(node1, node2);
                this.nodes.put(node2, node1);
                this.weight = weight;
                this.edgeId = edgeId;
                this.channels = new ArrayList<>(Collections.nCopies(P, -1));
            }

            public int getWeight() {
                return weight;
            }

            public int getEdgeId() {
                return edgeId;
            }

            public List<Integer> getChannels() {
                return channels;
            }

            public void setChannel(int channel, int serviceId) {
                this.channels.set(channel, serviceId);
            }

            public int getChannel(int channel) {
                return this.channels.get(channel);
            }

            public ArrayList<Integer> getNodes() {
                return new ArrayList<>(nodes.keySet());
            }

            public int getOtherNode(int node) {
                return nodes.get(node);
            }

            
        }
    public static class Graph {
        private int N;
        private int P;
        private Map<Integer, List<Integer>> adjList;
        private List<Edge> edges;
        private int edgeIdCounter;

        public Graph(int N, int P) {
            this.N = N;
            this.P = P;
            this.adjList = new HashMap<>();
            this.edges = new ArrayList<>();
            this.edgeIdCounter = 0;
            for (int i = 0; i < N; i++) {
                adjList.put(i, new ArrayList<>());
            }
        }

        public int addEdge(int node1, int node2, int weight) {
            int edgeId = edgeIdCounter;
            List<Integer> channels = new ArrayList<>();
            for (int i = 0; i < P; i++) {
                channels.add(-1);
            }
            Edge edge = new Edge(node1, node2, weight, edgeIdCounter, P);
            adjList.get(node1).add(edgeId);
            adjList.get(node2).add(edgeId);
            edges.add(edge);
            edgeIdCounter++;
            return edgeId;
        }

        public List<Integer> getNeighbors(int node) {
            Set<Integer> uniqueNeighbors = new HashSet<>();
            for (int edgeId : adjList.get(node)) {
                uniqueNeighbors.add(edges.get(edgeId).getOtherNode(node));
            }
            return new ArrayList<>(uniqueNeighbors);
        }

        public Edge getEdgeById(int edgeId) {
            return edges.get(edgeId);
        }

        public List<Integer> getNodesById(int edgeId) {
            Edge edge = getEdgeById(edgeId);
            return edge.getNodes();
        }

        public void addService(int serviceId, List<Integer> path, int channel) {
            for (int edgeId : path) {
                Edge edge = getEdgeById(edgeId);
                List<Integer> channels =  edge.getChannels();
                if (channels.get(channel) == -1) {
                    channels.set(channel, serviceId);
                } else {
                    System.out.println("Error: channel already used");
                    return;
                }
            }
        }

        public int getNodeConnected(int edge1Id, int edge2Id) {
            List<Integer> points1 = getNodesById(edge1Id);
            List<Integer> points2 = getNodesById(edge2Id);
            for (int p : points1) {
                for (int p2 : points2) {
                    if (p == p2) {
                        return p;
                    }
                }
            }
            return -1;
        }

        public void printGraph() {
            for (Edge edge : edges) {
                System.out.println(edge);
            }
        }


        public Map<Integer, Integer> dijkstraClassic(int source) {
            Map<Integer, Integer> distances = new HashMap<>();
            Map<Integer, Integer> previousNodes = new HashMap<>();
            PriorityQueue<Map.Entry<Integer, Integer>> priorityQueue = new PriorityQueue<>(Map.Entry.comparingByValue());

            for (int i = 0; i < N; i++) {
                distances.put(i, Integer.MAX_VALUE);
                previousNodes.put(i, null);
            }

            distances.put(source, 0);
            priorityQueue.add(new AbstractMap.SimpleEntry<>(source, 0));

            while (!priorityQueue.isEmpty()) {
                Map.Entry<Integer, Integer> entry = priorityQueue.poll();
                int current_node = entry.getKey();
                int current_distance = entry.getValue();

                if (current_distance > distances.get(current_node)) {
                    continue;
                }

                for (int edgeId : adjList.get(current_node)) {
                    Edge neighborData = edges.get(edgeId);
                    int neighbor =  neighborData.getOtherNode(current_node);
                    int weight =  neighborData.getWeight();
                    int distance = current_distance + weight;

                    if (distance < distances.get(neighbor)) {
                        distances.put(neighbor, distance);
                        previousNodes.put(neighbor, edgeId);
                        priorityQueue.add(new AbstractMap.SimpleEntry<>(neighbor, distance));
                    }
                }
            }

            return distances;
        }

        public List<Tuple3<Integer, List<Integer>, List<Boolean>>> dijkstraNPaths(int start, int end, int numPaths) {
            List<Tuple3<Integer, List<Integer>, List<Boolean>>> shortestPaths = new ArrayList<>();
            PriorityQueue<PathEntry> heap = new PriorityQueue<>(Comparator.comparingInt(PathEntry::getKey));
            heap.add(new PathEntry(0, start, new ArrayList<>(), new ArrayList<>(Collections.nCopies(P, true))));

            while (!heap.isEmpty()) {
                PathEntry entry = heap.poll();
                int currentCost = entry.getKey();
                int previousNode = entry.getValue();
                List<Integer> currentPath = new ArrayList<>(entry.getPath());
                List<Boolean> currentChannels = new ArrayList<>(entry.getChannels());

                if (previousNode == end) {
                    shortestPaths.add(new Tuple3<>(currentCost, currentPath, currentChannels));

                    if (shortestPaths.size() >= numPaths) {
                        break;
                    }
                }

                for (int edgeId : adjList.get(previousNode)) {
                    Edge neighborData = edges.get(edgeId);
                    int currentNeighbor = neighborData.getOtherNode(previousNode);
                    int weight = neighborData.getWeight();

                    List<Integer> channels = neighborData.getChannels();
                    List<Boolean> newChannels = new ArrayList<>();
                    for (int i = 0; i < P; i++) {
                        newChannels.add(currentChannels.get(i) && (channels.get(i) == -1));
                    }

                    if (newChannels.contains(true)) {
                        if (!currentPath.contains(edgeId)) {
                            List<Integer> newPath = new ArrayList<>(currentPath);
                            newPath.add(edgeId);
                            heap.add(new PathEntry(currentCost + weight, currentNeighbor, newPath, newChannels));
                        }
                    }
                }
            }

            return shortestPaths;
        }

    }

    public static class PathEntry implements Comparable<PathEntry> {
        private int key;
        private int value;
        private List<Integer> path;
        private List<Boolean> channels;

        public PathEntry(int key, int value, List<Integer> path, List<Boolean> channels) {
            this.key = key;
            this.value = value;
            this.path = path;
            this.channels = channels;
        }

        public int getKey() {
            return key;
        }

        public int getValue() {
            return value;
        }

        public List<Integer> getPath() {
            return path;
        }

        public List<Boolean> getChannels() {
            return channels;
        }

        @Override
        public int compareTo(PathEntry other) {
            return Integer.compare(this.key, other.key);
        }
    }

    public static class Tuple3<A, B, C> {
        private A first;
        private B second;
        private C third;

        public Tuple3(A first, B second, C third) {
            this.first = first;
            this.second = second;
            this.third = third;
        }

        public A getFirst() {
            return first;
        }

        public B getSecond() {
            return second;
        }

        public C getThird() {
            return third;
        }

        @Override
        public String toString() {
            return "(" + first + ", " + second + ", " + third + ")";
        }
    }

    public static class Tuple2<A, B> {
        private A first;
        private B second;

        public Tuple2(A first, B second) {
            this.first = first;
            this.second = second;
        }

        public A getFirst() {
            return first;
        }

        public B getSecond() {
            return second;
        }

        @Override
        public String toString() {
            return "(" + first + ", " + second + ")";
        }
    }

    public static List<Integer> sortService(List<Integer> servicesCost) {
        List<Map.Entry<Integer, Integer>> sortedServices = new ArrayList<>();
        for (int i = 0; i < servicesCost.size(); i++) {
            sortedServices.add(new AbstractMap.SimpleEntry<>(i, servicesCost.get(i)));
        }
        sortedServices.sort(Map.Entry.comparingByValue());
        List<Integer> sortedIndices = new ArrayList<>();
        for (Map.Entry<Integer, Integer> entry : sortedServices) {
            sortedIndices.add(entry.getKey());
        }
        return sortedIndices;
    }

    public static float costFunction(int distance, List<Integer> path, int D, Graph G, int P) {
        int alpha = 10000;
        int beta = 1;
        int gamma = 100;

        int maxNbChannels = Integer.MIN_VALUE;
        for (int edgeId : path) {
            List<Integer> channels = (List<Integer>) G.getEdgeById(edgeId).getChannels();
            int nbChannelsOfEdge = 0;
            for (int i=0; i<P;i++) {
                if (channels.get(i) > -1)
                    nbChannelsOfEdge += 1;
            }
            if (nbChannelsOfEdge > maxNbChannels) {
                maxNbChannels = nbChannelsOfEdge;
            }
        }

        return alpha * maxNbChannels + beta * path.size() + gamma * distance / D;
    }


    public static List<Integer> positionAmplifiers(List<Integer> path, Graph G, int D) {
        List<Integer> amplifiers = new ArrayList<>();
        int distance = 0;
        for (int i = 0; i < path.size(); i++) {
            int edge = path.get(i);
            int edgeDistance = (int) G.getEdgeById(edge).getWeight();

            if (distance + edgeDistance > D) {
                // Distance of an edge exceeds D, add an amplifier, reset distance to 0
                // We don't need to consider the case where i == 0 since the edge distance is less than D
                int node = G.getNodeConnected(edge, path.get(i - 1));
                amplifiers.add(node);
                distance = edgeDistance;
            } else {
                distance += edgeDistance;
            }
        }
        return amplifiers;
    }

    public static Tuple3<Integer, List<Integer>, Integer> verifyChannelForPath(List<Integer> path, Graph G, int P) {
        int[] channelsPath = new int[P];
        int lenPath = path.size();
        for (int i = 0; i < lenPath; i++) {
            Edge edge = G.getEdgeById(path.get(i));
            List<Integer> channels =  edge.getChannels();
            for (int j = 0; j < P; j++) {
                if (channels.get(j)==-1) {
                    channelsPath[j]++;
                }
            }
        }

        int channel = -1;
        int max_channel_posible = -1;
        for (int i=0; i<P;i++){
            if (channelsPath[i]>max_channel_posible) {
                max_channel_posible = channelsPath[i];
                channel = i;
            }
        }
        List<Integer> edgesToAdd = new ArrayList<>();

        for (int i = 0; i < lenPath; i++) {
            Edge edge = G.getEdgeById(path.get(i));
            List<Integer> channels = edge.getChannels();
            if (channels.get(channel)>-1) {
                edgesToAdd.add(path.get(i));
            }
        }

        return new Tuple3<>(channel, edgesToAdd, lenPath - max_channel_posible);
    }

    public static Tuple3<List<Integer>, List<Integer>, Integer> chooseNewEdge(Graph G, Tuple2<Integer, Integer> service, Graph GInitial, int N, int P) {
        List<Tuple3<Integer, List<Integer>, List<Boolean>>> res = GInitial.dijkstraNPaths(service.getFirst(), service.getSecond(), N);
        int minNumNewEdge = Integer.MAX_VALUE;
        int minDistance = Integer.MAX_VALUE;
        List<Integer> newEdges = new ArrayList<>();
        List<Integer> bestPath = new ArrayList<>();
        int bestChannelsPossible = 0;

        for (Tuple3<Integer, List<Integer>, List<Boolean>> r : res) {
            List<Integer> newPath = new ArrayList<>(r.getSecond());
            int distance = r.getFirst();
            Tuple3<Integer, List<Integer>, Integer> verifyResult = verifyChannelForPath(newPath, G, P);
            int numNewEdge = verifyResult.getThird();
            List<Integer> newEdgeIndices = verifyResult.getSecond();
            int numChannelsPossible = verifyResult.getFirst();

            if (numNewEdge < minNumNewEdge || (numNewEdge == minNumNewEdge && distance < minDistance)) {
                minNumNewEdge = numNewEdge;
                minDistance = distance;
                newEdges = newEdgeIndices;
                bestPath = newPath;
                bestChannelsPossible = numChannelsPossible;
            }
        }

        Map<Integer, Integer> edgeReplace = new HashMap<>();
        for (int i : newEdges) {
            Edge oldEdge = G.getEdgeById(i);
            List<Integer> nodes = G.getNodesById(i);
            int newId = G.addEdge(nodes.get(0), nodes.get(1), (int) oldEdge.getWeight());
            edgeReplace.put(i, newId);
        }

        for (int i = 0; i < bestPath.size(); i++) {
            if (edgeReplace.containsKey(bestPath.get(i))) {
                bestPath.set(i, edgeReplace.get(bestPath.get(i)));
            }
        }

        return new Tuple3<>(newEdges, bestPath, bestChannelsPossible);
    }

    public static Graph createGraph(int N, int M, int P, List<Tuple3<Integer,Integer,Integer>> edges) {
        Graph graph = new Graph(N, P);
        for (int i = 0; i < M; i++) {
            Tuple3<Integer,Integer,Integer> edge = edges.get(i);
            int s = edge.getFirst();
            int t = edge.getSecond();
            int d = edge.getThird();
            graph.addEdge(s, t, d);
        }
        return graph;
    }


    public static void main(String[] args) {
        // Scanner scanner = new Scanner(System.in);

        // // 读取 N，M，T，P，D
        // int N = scanner.nextInt();
        // int M = scanner.nextInt();
        // int T = scanner.nextInt();
        // int P = scanner.nextInt();
        // int D = scanner.nextInt();

        // // 读取边信息
        // List<Tuple3<Integer,Integer,Integer>> edges = new ArrayList<>();
        // for (int i = 0; i < M; i++) {
        //     int ci = scanner.nextInt();
        //     int si = scanner.nextInt();
        //     int ti = scanner.nextInt();
        //     int di = scanner.nextInt();
        //     edges.add(new Tuple3<>(si, ti, di));
        // }

        // // 读取服务信息
        // List<Tuple2<Integer, Integer>> services = new ArrayList<>();
        // for (int i = 0; i < T; i++) {
        //     int Sj = scanner.nextInt();
        //     int Tj = scanner.nextInt();
        //     services.add(new Tuple2<>(Sj, Tj));
        // }

        // scanner.close();

        int N = 1000;
        int M = (N*(N/2))/2;
        int T = 100;
        int P = 4;
        int D = 100;

        long startTime = System.currentTimeMillis();

        // Code you want to measure the execution time for
        // ...

        
        List<Tuple3<Integer,Integer,Integer>> edges = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            int ci = i;
            int si = (int) (Math.random() * N);
            int ti = (int) (Math.random() * N);
            int di = (int) (Math.random() * 100);
            edges.add(new Tuple3<>(si, ti, di));
        }

        

        Graph graph = createGraph(N, M, P, edges);
        Graph graphInitial = createGraph(N, M, P, edges);

        List<Tuple2<Integer, Integer>> services = new ArrayList<>();
        for (int i = 0; i < T; i++) {
            int Sj = (int) (Math.random() * N);
            int Tj = (int) (Math.random() * N);
            services.add(new Tuple2<>(Sj, Tj));
        }
        for(int i = 0; i<services.size();i++){
            Tuple2<Integer,Integer> service = services.get(i);
            int start = service.getFirst();
            int end = service.getSecond();
            List<Tuple3<Integer, List<Integer>, List<Boolean>>> dikstra_res = graphInitial.dijkstraNPaths(start, end, 1);
            System.out.println(i);
        }

        // Record the end time
        long endTime = System.currentTimeMillis();

        // Calculate the execution time in milliseconds
        long executionTime = endTime - startTime;

        System.out.println("Execution time: " + executionTime );


 
        // int nbPath1 = 1; // Number of shortest paths to find for the first dijkstra
        // int nbPath2 = 1; // Number of shortest paths to find for the second dijkstra



        // Map<Integer, List<Integer>> newEdge = new HashMap<>();
        // List<List<Object>> pathComplete = new ArrayList<>();

        // List<Integer> servicesCostBrut = new ArrayList<>();

        // for(int i = 0; i<services.size();i++){

        //     Tuple2<Integer,Integer> service = services.get(i);
        //     int start = service.getFirst();
        //     int end = service.getSecond();

        //     List<Tuple3<Integer, List<Integer>, List<Boolean>>> dikstra_res = graphInitial.dijkstraNPaths(start, end, 1);

        //     Tuple3<Integer, List<Integer>, List<Boolean>> r = dikstra_res.get(0);

        //     List<Integer> newPath = new ArrayList<>(r.getSecond());

        //     int distance = r.getFirst();

        //     Tuple3<Integer, List<Integer>, Integer> verifyResult = verifyChannelForPath(newPath, graph, P);

        //     int numNewEdge = verifyResult.getThird();

        //     List<Integer> newEdgeIndices = verifyResult.getSecond();

        //     int numChannelsPossible = verifyResult.getFirst();

        //     Map<Integer, Integer> edgeReplace = new HashMap<>();

        //     for (int j : newEdgeIndices) {

        //         Edge oldEdge = graph.getEdgeById(j);
        //         List<Integer> nodes = graph.getNodesById(j);
        //         int newId = graph.addEdge(nodes.get(0), nodes.get(1), oldEdge.getWeight());
        //         edgeReplace.put(j, newId);
        //         newEdge.put(newId, Arrays.asList(nodes.get(0), nodes.get(1)));

        //     }

        //     for (int j = 0; j < newPath.size(); j++) {
        //         if (edgeReplace.containsKey(newPath.get(j))) {
        //             newPath.set(j, edgeReplace.get(newPath.get(j)));
        //         }
        //     }


        //     List<Integer> posAmp = positionAmplifiers(newPath, graph, D);

        //     List<Object> pathComp = new ArrayList<>();

        //     pathComp.add(numChannelsPossible);
        //     pathComp.add(newPath.size());
        //     pathComp.add(posAmp.size());

        //     for (int j=0; j<newPath.size();j++)
        //         pathComp.add(newPath.get(j));

        //     for (int k=0; k<posAmp.size();k++)
        //         pathComp.add(posAmp.get(k));

        //     pathComplete.add(pathComp);
        //     graph.addService(i, newPath, numChannelsPossible);
        
        
        // }

        // System.out.println(newEdge.size());
        // for (List<Integer> values : newEdge.values()) {
        //     System.out.print(values.get(0)+" "+values.get(1)+"\n");
        // }
        // for (int i=0; i< pathComplete.size();i++) {
        //     for (int j=0;j<pathComplete.get(i).size();j++)
        //         System.out.print(pathComplete.get(i).get(j)+" ");
        //     System.out.println("");
        // }
    }

}


