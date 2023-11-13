import java.util.*;

public class Main {

    /*
     * @param list: a list of integers
     * @return: return index of the minimum element
     */

    public static int find_min_value_index(List<Integer> list) {
        int minIndex = 0;
        int minValue = list.get(0); // Assume the first element is the minimum

        for (int i = 1; i < list.size(); i++) {
            if (list.get(i) < minValue) {
                minValue = list.get(i);
                minIndex = i;
            }
        }
        return minIndex;
    }

    public static Map<Integer, Map<Node, Integer>> matrice_path(Graph graphe) {
        Map<Integer, Map<Node, Integer>> mat_path = new HashMap<>();
        for (int i = 0; i < graphe.N; i++) {
            mat_path.put(i, graphe.dijkstra_for_ranking(i));
        }
        return mat_path;
    }

    public static void main(String[] args) {

        Scanner scanner = new Scanner(System.in);

        // read N，M，T，P，D
        int N = scanner.nextInt();
        int M = scanner.nextInt();
        int T = scanner.nextInt();
        int P = scanner.nextInt();
        int D = scanner.nextInt();

        Graph graph = new Graph(N, P, D);

        for (int i = 0; i < M; i++) {
            int ci = scanner.nextInt();
            int si = scanner.nextInt();
            int ti = scanner.nextInt();
            int di = scanner.nextInt();
            graph.addEdgeFromNodeId(si, ti, di);//0
        }

        // read service
//            LinkedHashMap<Integer, List<Integer>> services = new LinkedHashMap<>();
//            for (int i = 0; i < T; i++) {
//                int Sj = scanner.nextInt();
//                int Tj = scanner.nextInt();
//                services.put(i, List.of(Sj, Tj));
//            }
//            List<List<Integer>> service_ranking = rankingService(services, graph);
//
//            for (int i=0; i<T;i++) {
//                int serv_id = service_ranking.get(i).get(0);
//                int Sj = services.get(serv_id).get(0);
//                int Tj = services.get(serv_id).get(1);
//                graph.dijkstra(Sj, Tj, i);
//            }


        for (int i = 0; i < T; i++) {
            int Sj = scanner.nextInt();
            int Tj = scanner.nextInt();
            graph.dijkstra(Sj, Tj, i);
        }

        scanner.close();
        System.out.println(graph.getResult());

    }

    public static class Node {
        private final int id;
        private int nb_neighbors;
        // TODO
        //private  int heuristic;
        private int P;
        private List<Integer> channels_counter;

        /*
         * @param id: unique id of the node
         * @param P: number of channels
         * @return: a node
         * @description: constructor of a node
         */
        public Node(int id, int P) {
            this.id = id;
            this.nb_neighbors = 0;
            this.channels_counter = new ArrayList<>(Collections.nCopies(P, 0));
            this.P = P;
        }

        public List<Integer> get_channels_counter() {
            return this.channels_counter;
        }

        public void set_channels_counter(List<Integer> channels) {
            this.channels_counter = channels;

        }

        public void init_channel_counter() {
            this.channels_counter = new ArrayList<>(Collections.nCopies(P, 0));
        }

        public int getId() {
            return id;
        }

        public int getNb_neighbors() {
            return nb_neighbors;
        }

        /*
         * @description: increase the number of neighbors
         */
        public void addNeighbor() {
            this.nb_neighbors++;
        }

        /*
         * @description: decrease the number of neighbors
         */
        public void removeNeighbor() {
            this.nb_neighbors--;
        }

        /*
         * @param o: another object
         * @return: whether two nodes are equal
         * @description: check whether two nodes are equal
         */
        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;
            Node myNode = (Node) o;
            return this.id == myNode.id;
        }

        /*
         * @return: hash code of the node
         * @description: get the hash code of the node
         */
        @Override
        public int hashCode() {
            return Objects.hash(id);
        }

        /*
         * @return: string representation of the node
         * @description: get the string representation of the node
         */
        @Override
        public String toString() {
            return String.valueOf(id);
        }
    }

    public static class Edge {
        private final Map<Node, Node> nodes = new HashMap<>();
        private final int weight;
        private final int id;
        private final List<Integer> channels;
        private int nb_channels_used;

        /*
         * @param node1: one end of the edge
         * @param node2: the other end of the edge
         * @param weight: weight of the edge
         * @param edgeId: unique id of the edge
         * @param P: number of channels
         * @return: an edge
         * @description: constructor of an edge
         */
        public Edge(Node node1, Node node2, int weight, int edgeId, int P) {
            this.nodes.put(node1, node2);
            this.nodes.put(node2, node1);
            this.weight = weight;
            this.id = edgeId;
            this.nb_channels_used = 0;
            this.channels = new ArrayList<>(Collections.nCopies(P, -1));
        }

        public int getWeight() {
            return weight;
        }

        public int getId() {
            return id;
        }

        public List<Integer> getChannels() {
            return channels;
        }

        /*
         * @param channel: channel id
         * @description: set the channel to be used, if the channel is already used, throw an exception
         */
        public void setChannels(int channel) {
            if (channels.get(channel) == -1) {
                this.channels.set(channel, 1);
                nb_channels_used++;
                if (nb_channels_used == 0) {
                    for (Node n : getNodes()) {
                        n.removeNeighbor();
                    }
                }
            } else {
                throw new UnsupportedOperationException("Channel already used");
            }
        }

        /*
         * @param channel: channel id
         * @return: whether the channel is used, 1 for used, -1 for not used
         */
        public int getChannel(int channel) {
            return this.channels.get(channel);
        }

        public int getNbChannelUsed() {
            return this.nb_channels_used;
        }

        /*
         * @return: the two nodes of the edge
         */
        public ArrayList<Node> getNodes() {
            return new ArrayList<>(nodes.keySet());
        }

        /*
         * @param node: one end of the edge
         * @return: the other end of the edge
         */
        public Node getOtherNode(Node node) {
            return nodes.get(node);
        }


        /*
         * @return: string representation of the edge
         * @description: get the string representation of the edge
         */
        @Override
        public String toString() {
            List<Node> nodes = getNodes();
            return nodes.get(0) + " " + nodes.get(1);
        }
    }

    public static class Service {
        private final int id;
        private List<Edge> path;
        private List<Integer> amplifiers;

        private int channelId;

        /*
         * @param id: unique id of the service
         * @param path: path of the service
         * @param amplifiers: amplifiers of the service
         * @param channelId: channel id of the service
         */
        public Service(int id, List<Edge> path, List<Integer> amplifiers, int channelId) {
            this.id = id;
            this.path = path;
            this.amplifiers = amplifiers;
            this.channelId = channelId;
        }

        public int getId() {
            return id;
        }

        public List<Edge> getPath() {
            return path;
        }

        /*
         * @param channelId: channel id of the service
         */
        public void setPath(List<Edge> path) {
            this.path = path;
        }

        public List<Integer> getAmplifiers() {
            return amplifiers;
        }

        /*
         * @param amplifiers: positions of amplifiers of the service
         */
        public void setAmplifiers(List<Integer> amplifiers) {
            this.amplifiers = amplifiers;
        }

        public int getPathLength() {
            return path.size();
        }

        public int getChannelId() {
            return channelId;
        }
    }

    public static class Result {
        private final List<Edge> addedEdges;
        private final List<Service> services;

        /*
         * @description: constructor of the result, it contains the added edges and the services
         */
        public Result() {
            this.addedEdges = new ArrayList<>();
            this.services = new ArrayList<>();
        }

        /*
         * @param edge: an edge
         * @description: add an edge to the result(addedEdges)
         */
        public void addEdge(Edge edge) {
            addedEdges.add(edge);
        }

        /*
         * @param service: a service
         * @description: add a service to the result(services)
         */
        public void addService(Service service) {
            services.add(service);
        }

        /*
         * @return: string representation of the result
         * @description: get the string representation of the result
         */
        @Override
        public String toString() {
            StringBuilder output = new StringBuilder();
            int nb_addedEdges = addedEdges.size();
            output.append(nb_addedEdges).append("\n");
            for (Edge edge : addedEdges) {
                output.append(edge.toString()).append('\n');
            }

            // add service by son order id
//                Comparator<Service> idComparator = Comparator.comparingInt(Service::getId);
//                Collections.sort(services, idComparator);
            int nb_amp = 0;
            for (Service service : services) {
                output.append(service.getChannelId()).append(" ");
                output.append(service.getPathLength()).append(" ");
                List<Integer> amplifiers = service.getAmplifiers();
                nb_amp += amplifiers.size();
                output.append(amplifiers.size()).append(" ");
                for (Edge edge : service.getPath()) {
                    output.append(edge.getId()).append(" ");
                }

                for (int i = 0; i < amplifiers.size(); i++) {
                    if (i != amplifiers.size() - 1) {
                        output.append(amplifiers.get(i)).append(' ');
                    } else {
                        output.append(amplifiers.get(i));
                    }
                }

                output.append("\n");
            }
            return output.toString();
        }
    }

    public static class Graph {
        private final int D;
        private final Map<Node, List<Edge>> adjList;
        private final List<Node> nodes;
        private final List<Edge> edges;
        private int N;
        private int P;
        private int edgeIdCounter;
        private Result result;

        /*
         * @param N: number of nodes
         * @param P: number of channels
         * @param D: maximum distance between two amplifiers
         * @description: constructor of the graph
         */
        public Graph(int N, int P, int D) {
            this.N = N;
            this.P = P;
            this.D = D;

            this.nodes = new ArrayList<>();
            this.edges = new ArrayList<>();
            this.adjList = new HashMap<>();

            for (int i = 0; i < N; i++) {
                Node n = new Node(i, P);
                this.nodes.add(n);
                adjList.put(n, new ArrayList<>());
            }
            this.edgeIdCounter = 0;
            this.result = new Result();
        }

        /*
         * @param node1: one end of the edge
         * @param node2: the other end of the edge
         * @param weight: weight of the edge
         * @return: the edge
         * @description: add an edge to the graph (used for adding new edges when solving the problem)
         */
        public Edge addEdge(Node node1, Node node2, int weight) {
            node1.addNeighbor();
            node2.addNeighbor();
            Edge edge = new Edge(node1, node2, weight, edgeIdCounter, P);
            adjList.get(node1).add(edge);
            adjList.get(node2).add(edge);
            edges.add(edge);
            edgeIdCounter++;
            return edge;
        }

        /*
         * @param node1: one end of the edge
         * @param node2: the other end of the edge
         * @param weight: weight of the edge
         * @description: add an edge to the graph (used for reading input)
         */
        public void addEdgeFromNodeId(int node1, int node2, int weight) {
            Node n1 = nodes.get(node1);
            Node n2 = nodes.get(node2);
            addEdge(n1, n2, weight);
        }

        /*
         * @param node: a node
         * @return: the neighbors of the node
         * @description: get the neighbors of the node
         */
        public List<Node> getNeighbors(Node node) {
            Set<Node> uniqueNeighbors = new HashSet<>();
            for (Edge edge : adjList.get(node)) {
                uniqueNeighbors.add(edge.getOtherNode(node));
            }
            return new ArrayList<>(uniqueNeighbors);
        }

        /*
         * @param edgeId: unique id of the edge
         * @return: the edge
         * @description: get the edge by its id
         */
        public Edge getEdgeById(int edgeId) {
            return edges.get(edgeId);
        }

        public Result getResult() {
            return result;
        }

        /*
         * @param edge: an edge
         * @return: the two nodes of the edge
         * @description: get the two nodes of the edge
         */
        public List<Node> getNodesById(Edge edge) {
            return edge.getNodes();
        }

        /*
         * @param serviceId: unique id of the service
         * @param path: path of the service
         * @param channel: channel id of the service
         * @description: add a service to the result
         */
        public void addService(int serviceId, List<Edge> path, int channel) {

            for (Edge edge : path) {
                edge.setChannels(channel);
            }
        }

        //        public int getNodeConnected(int edge1Id, int edge2Id) {
//            List<Integer> points1 = getNodesById(edge1Id);
//            List<Integer> points2 = getNodesById(edge2Id);
//            for (int p : points1) {
//                for (int p2 : points2) {
//                    if (p == p2) {
//                        return p;
//                    }
//                }
//            }
//            return -1;
//        }
        public void printGraph() {
            for (Edge edge : edges) {
                System.out.println(edge);
            }
        }

        /*
         * @param previousEdge: a map containing the previous edge of each node
         * @param start: start node
         * @param end: end node
         * @param service_id: unique id of the service
         * @description: reconstruct the path of the service from node start to node end
         */
        private void reconstructPath(Map<Node, Edge> previousEdge, Node start, Node end, int service_id) {
            List<Edge> path = new ArrayList<>();
            Edge edge = previousEdge.get(end);
            List<Integer> added_amplifiers = new ArrayList<>();
            List<Integer> count_channels = new ArrayList<>(Collections.nCopies(P, 0));

            Node previousNode = end;
            int distance = 0;
            while (previousNode != start) {

                int w = edge.getWeight();
                distance += w;
                if (distance > D) {
                    added_amplifiers.add(previousNode.getId());
                    distance = w;
                }
                previousNode = edge.getOtherNode(previousNode);

                List<Integer> c = edge.getChannels();
                for (int i = 0; i < P; i++) {
                    count_channels.set(i, count_channels.get(i) + c.get(i));
                }
                path.add(edge);
                edge = previousEdge.get(previousNode);
            }
            int channel_used = find_min_value_index(count_channels);
            Collections.reverse(path);

            for (int i = 0; i < path.size(); i++) {
                Edge path_edge = path.get(i);
                if (path_edge.getChannel(channel_used) == -1) {
                    path_edge.setChannels(channel_used);
                } else {
                    List<Node> n = path_edge.getNodes();
                    path_edge = addEdge(n.get(0), n.get(1), path_edge.getWeight());
                    path_edge.setChannels(channel_used);
                    path.set(i, path_edge);
                    result.addEdge(path_edge);
                }
            }

            Collections.reverse(added_amplifiers);
            Service service = new Service(service_id, path, added_amplifiers, channel_used);
            result.addService(service);
        }

        public void dijkstra_old(int s, int e, int service_id) {
            Node start = nodes.get(s);
            Node end = nodes.get(e);

            Map<Node, Integer> distance = new HashMap<>();
            Map<Node, Edge> previousEdge = new HashMap<>();
            PriorityQueue<Node> minHeap = new PriorityQueue<>(Comparator.comparingInt(distance::get));

            for (Node node : nodes) {
                distance.put(node, node.equals(start) ? 0 : Integer.MAX_VALUE);
                minHeap.offer(node);
            }

            while (!minHeap.isEmpty()) {
                Node current = minHeap.poll();

                if (current.equals(end)) {
                    reconstructPath(previousEdge, start, current, service_id);
                    return;
                }

                for (Edge edge : adjList.get(current)) {
                    Node neighbor = edge.getOtherNode(current);
                    int newDistance = distance.get(current) + edge.getWeight();
                    if (newDistance < distance.get(neighbor)) {
                        distance.put(neighbor, newDistance);
                        previousEdge.put(neighbor, edge);
                        minHeap.remove(neighbor);
                        minHeap.offer(neighbor);
                    }
                }
            }
        }

        public Map<Node, Integer> dijkstra_for_ranking(int s) {
            Node start = nodes.get(s);

            Map<Node, Integer> distance = new HashMap<>();
            PriorityQueue<Node> minHeap = new PriorityQueue<>(Comparator.comparingInt(distance::get));

            for (Node node : nodes) {
                distance.put(node, node.equals(start) ? 0 : Integer.MAX_VALUE);
                minHeap.offer(node);
            }

            while (!minHeap.isEmpty()) {
                Node current = minHeap.poll();


                for (Edge edge : adjList.get(current)) {
                    Node neighbor = edge.getOtherNode(current);
                    int newDistance = distance.get(current) + 1;
                    if (newDistance < distance.get(neighbor)) {
                        distance.put(neighbor, newDistance);
                        minHeap.remove(neighbor);
                        minHeap.offer(neighbor);
                    }
                }
            }
            return distance;
        }

        public void dijkstra(int s, int e, int service_id) {
            Node start = nodes.get(s);
            Node end = nodes.get(e);

            Map<Node, Integer> distance = new HashMap<>();
            Map<Node, Edge> previousEdge = new HashMap<>();//important
            PriorityQueue<Node> minHeap = new PriorityQueue<>(Comparator.comparingInt(distance::get));

            for (Node node : nodes) {
                distance.put(node, node.equals(start) ? 0 : Integer.MAX_VALUE);
                minHeap.offer(node);
            }

            while (!minHeap.isEmpty()) {
                Node current = minHeap.poll();

                if (current.equals(end)) {
                    // condition stop
                    reconstructPath(previousEdge, start, end, service_id);
                    return;
                }

                // compare the available channels with the previous edge
                List<Integer> previous_channels = new ArrayList<>(Collections.nCopies(P, 0));
                Edge prev_edge = previousEdge.get(current);
                if (prev_edge != null) {
                    previous_channels = prev_edge.getChannels();
                }

                for (Edge edge : adjList.get(current)) {
                    Node neighbor = edge.getOtherNode(current);

                    // count the number of not available channels in the previous edge and current edge
                    int nb_not_available_channels = 0;
                    List<Integer> current_channels = edge.getChannels();
                    List<Integer> counter_channels = new ArrayList<>(Collections.nCopies(P, 0));
                    for (int i = 0; i < P; i++) {
                        int r = current_channels.get(i) + previous_channels.get(i);
                        if (r > -2)
                            counter_channels.set(i, 1);
                    }

                    List<Integer> counter_channels_neig = new ArrayList<>(Collections.nCopies(P, 0));
                    int nb_neig = adjList.get(neighbor).size();
                    for (Edge edge_neig : adjList.get(neighbor)) {
                        if (edge_neig == null) {
                            nb_not_available_channels = 10;
                        } else {
                            List<Integer> l = edge_neig.getChannels();
                            for (int k = 0; k < P; k++) {
                                if (l.get(k) == 1)
                                    counter_channels_neig.set(k, counter_channels_neig.get(k) + 1);
                            }
                        }
                    }

                    for (int k = 0; k < P; k++) {
                        if (counter_channels.get(k) == 0 && counter_channels_neig.get(k) >= nb_neig - 1) {
                            nb_not_available_channels++;
                        }
                        if (counter_channels.get(k) > 0)
                            nb_not_available_channels += 1;

                    }

                    int somme = 0;
                    for (int k = 0; k < P; k++) {
                        somme += counter_channels_neig.get(k);
                    }
                    int newDistance = distance.get(current) + (edge.getNbChannelUsed()) + nb_not_available_channels + neighbor.getNb_neighbors();

//                        if(somme>0) {
//                            newDistance += 10000 / somme / P;
//                        }

                    if (newDistance < distance.get(neighbor)) {
                        // update node (channels)
                        distance.put(neighbor, newDistance);
                        previousEdge.put(neighbor, edge);
                        minHeap.remove(neighbor);
                        minHeap.offer(neighbor);
                    }
                }
            }
        }
    }
}


