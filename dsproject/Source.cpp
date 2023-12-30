#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include<algorithm>
#include <iomanip>
#include <queue>
#include <SFML/Graphics.hpp>

using namespace std;

class Random {
public:
    Random() {
        srand(time(nullptr));
    }

    double nextDouble() {
        return static_cast<double>(rand()) / RAND_MAX;
    }
};

struct Node {
    int nodeNum;
    double x, y;
};

class WirelessNetwork {
private:
    int numNodes;
    double size;
    int totaldeg, maxdeg;
    float avgdeg;
    vector<Node> nodes;
    vector<vector<bool>> adjacencyMatrix;

public:
    WirelessNetwork(int numNodes, double size) : numNodes(numNodes), size(size) {
        Random random;

        // Generate random points within the square
        for (int i = 0; i < numNodes; i++) {
            Node node;
            node.nodeNum = i;
            node.x = random.nextDouble() * size;
            node.y = random.nextDouble() * size;
            nodes.push_back(node);
        }

        // Create the adjacency matrix
        adjacencyMatrix.resize(numNodes, vector<bool>(numNodes, false));

        // Create edges within distance 1 unit
        for (int i = 0; i < numNodes; i++) {
            for (int j = i + 1; j < numNodes; j++) {
                if (calculateDistance(nodes[i], nodes[j]) <= 1.0) {
                    adjacencyMatrix[i][j] = true;
                    adjacencyMatrix[j][i] = true;
                }
            }
        }
    }

    // Default constructor
    WirelessNetwork() : WirelessNetwork(500, 10.0) {} // Using default parameters
    double get_maxdeg() {
        return maxdeg;
    }
    double get_avgdeg() {
        return avgdeg;
    }
    void degreeCalc() {
        int temp;
        totaldeg = 0;
        avgdeg = 0.0;
        maxdeg = 0;
        for (int i = 0; i < numNodes; i++) {
            temp = 0;
            for (int j = 0; j < numNodes; j++) {
                if (adjacencyMatrix[i][j] == 1) {
                    totaldeg += 1;
                    temp += 1;
                }
            }
            if (temp > maxdeg) {
                maxdeg = temp;
            }
        }
        avgdeg = (float)totaldeg / numNodes;
    }

    double calculateDistance(const Node& node1, const Node& node2) {
        double dx = node1.x - node2.x;
        double dy = node1.y - node2.y;
        return sqrt(dx * dx + dy * dy);
    }

    double calculateDistanceint(int s, int t) {
        double dx = nodes[s].x - nodes[t].x;
        double dy = nodes[s].y - nodes[t].y;
        return sqrt(dx * dx + dy * dy);
    }

    vector<int> getNeighborsx(int nodeNum) {
        vector<int> neighbors;
        for (int i = 0; i < numNodes; i++) {
            if (adjacencyMatrix[nodeNum][i]) {
                neighbors.push_back(i);
            }
        }
        return neighbors;
    }

    void printAdjacencyList() {
        degreeCalc();
        for (int i = 0; i < numNodes; i++) {
            cout << "Node " << i << " connected to: ";
            for (int j = 0; j < numNodes; j++) {
                if (adjacencyMatrix[i][j]) {
                    cout << j << " ";
                }
            }
            cout << endl;
        }
        cout << "Total degree: " << totaldeg << endl;
        cout << "Max degree: " << maxdeg << endl;
        cout << "Avg Degree: " << avgdeg << endl;
    }

    void deleteEdge(int u, int v) {
        adjacencyMatrix[u][v] = 0;
        adjacencyMatrix[v][u] = 0;
    }

    vector<Node> compassRouting(int s, int t, int depthcontrol = 100) {
        vector<Node> neighb = getNeighbors(s);
        float angle, minangle = 90;
        double d1, d2, dmin;
        if (s == t) {
            return { nodes[s] };
        }
        else {
            int nextVertex = -1;

            vector<double> v1;
            v1.push_back(nodes[t].x - nodes[s].x);
            v1.push_back(nodes[t].y - nodes[s].y);

            for (int i = 0; i < neighb.size(); i++) {
                vector<double> v2 = {};
                v2.push_back(neighb[i].x - nodes[s].x);
                v2.push_back(neighb[i].y - nodes[s].y);

                angle = angleBetweenVectors(v1, v2);

                if (angle < minangle) {
                    nextVertex = neighb[i].nodeNum;
                    minangle = angle;
                }
            }
            if (nextVertex == -1 || depthcontrol <= 0) {
                return {};
            }
            vector<Node> path = compassRouting(nextVertex, t, depthcontrol - 1);
            path.insert(path.begin(), nodes[s]);
            return path;
        }
    }

    vector<Node> getNeighbors(int s) {
        vector<Node> neighbors;
        for (int i = 0; i < numNodes; i++) {
            if (adjacencyMatrix[s][i]) {
                neighbors.push_back(nodes[i]);
            }
        }
        return neighbors;
    }

    double dotProduct(const vector<double>& vector1, const vector<double>& vector2) {
        double dotProduct = 0.0;
        for (int i = 0; i < vector1.size(); i++) {
            dotProduct += vector1[i] * vector2[i];
        }
        return dotProduct;
    }

    double magnitude(const vector<double>& vector) {
        double magnitude = 0.0;
        for (double value : vector) {
            magnitude += value * value;
        }
        return sqrt(magnitude);
    }

    double angleBetweenVectors(const vector<double>& vector1, const vector<double>& vector2) {
        double dot = dotProduct(vector1, vector2);
        double mag1 = magnitude(vector1);
        double mag2 = magnitude(vector2);

        return acos(dot / (mag1 * mag2));
    }


    void sortNeighbors(vector<int>& neighbors, const Node& currentNode) {
        int n = neighbors.size();

        for (int i = 0; i < n - 1; i++) {
            int minIndex = i;

            for (int j = i + 1; j < n; j++) {
                double distA = calculateDistance(currentNode, nodes[neighbors[j]]);
                double distB = calculateDistance(currentNode, nodes[neighbors[minIndex]]);

                if (distA < distB) {
                    minIndex = j;
                }
            }

            if (minIndex != i) {
                swap(neighbors[i], neighbors[minIndex]);
            }
        }
    }

    float routeDist(const vector<Node>& path) {
        float dist = 0;
        for (int i = 0; i < path.size() - 1; i++) {
            Node current = path[i];
            Node next = path[i + 1];
            dist += calculateDistance(current, next);
        }
        return dist;
    }

    //void TopologyControl() {
    //    vector<vector<bool>> updatedAdjacencyMatrix(numNodes, vector<bool>(numNodes, false));

    //    for (int i = 0; i < numNodes; i++) {
    //        vector<int> neighbors_u;    // Set of neighbors of u [N(u)]
    //        vector<int> neighbors_v;    // Set of neighbors of v ∈ N(u)
    //        vector<int> temp;           // Temporary vector: it'll be a copy of neighbors of u.

    //        double dist1, dist2, dist3; // dist1: |u,w|; w ∈ N(u) & N(v)
    //        // dist2: |u,v|
    //        // dist3: |v,w|

    //        neighbors_u = getNeighborsx(i);
    //        sortNeighbors(neighbors_u, nodes[i]);
    //        temp = getNeighborsx(i);
    //        sortNeighbors(temp, nodes[i]);
    //        if (!neighbors_u.empty()) {
    //            int maxEdges = min(static_cast<int>(neighbors_u.size()), 7);
    //            for (int j = 0; j < maxEdges; j++) {
    //                neighbors_v = getNeighborsx(neighbors_u[j]);

    //                if (!neighbors_v.empty()) {
    //                    for (int k = 0; k < temp.size(); k++) {
    //                        for (int l = 0; l < neighbors_v.size(); l++) {
    //                            if (temp[k] == neighbors_v[l]) {
    //                                dist1 = calculateDistance(nodes[i], nodes[neighbors_v[l]]);    // dist1: |u,w|
    //                                dist2 = calculateDistance(nodes[i], nodes[neighbors_u[j]]);    // dist2: |u,v|
    //                                dist3 = calculateDistance(nodes[neighbors_u[j]], nodes[neighbors_v[l]]);    // dist3: |v,w|

    //                                if ((dist1 < dist2) && (dist3 < dist2)) {
    //                                    updatedAdjacencyMatrix[i][neighbors_u[j]] = true;
    //                                }
    //                            }
    //                        }
    //                    }
    //                }
    //            }
    //        }
    //    }

    //    adjacencyMatrix = updatedAdjacencyMatrix;
    //}

    void XTCTopologyControl() {
        vector<vector<bool>> updatedAdjacencyMatrix(numNodes, vector<bool>(numNodes, false));

        for (int i = 0; i < numNodes; i++) {
            vector<pair<double, int>> neighborDistances;
            for (int j = 0; j < numNodes; j++) {
                if (i != j) {
                    double distance = calculateDistance(nodes[i], nodes[j]);
                    neighborDistances.push_back({ distance, j });
                }
            }

            sort(neighborDistances.begin(), neighborDistances.end());

            for (int j = 0; j < 6; j++) {
                int neighbor = neighborDistances[j].second;
                updatedAdjacencyMatrix[i][neighbor] = true;
                updatedAdjacencyMatrix[neighbor][i] = true;
            }
        }

        // Remove extra edges if any node has more than 6 neighbors
        for (int i = 0; i < numNodes; i++) {
            vector<int> neighbors = getNeighborsx(i);
            if (neighbors.size() > 6) {
                sortNeighbors(neighbors, nodes[i]);
                for (int j = 6; j < neighbors.size(); j++) {
                    int neighbor = neighbors[j];
                    updatedAdjacencyMatrix[i][neighbor] = false;
                    updatedAdjacencyMatrix[neighbor][i] = false;
                }
            }
        }

        adjacencyMatrix = updatedAdjacencyMatrix;
    }
    bool isConnected() {
        // Mark all nodes as not visited
        vector<bool> visited(numNodes, false);

        // Create a queue for BFS
        queue<int> bfsQueue;

        // Start with the first node
        int startNode = 0;

        // Mark the start node as visited and enqueue it
        visited[startNode] = true;
        bfsQueue.push(startNode);

        while (!bfsQueue.empty()) {
            // Dequeue a node from the queue
            int currentNode = bfsQueue.front();
            bfsQueue.pop();

            // Get all adjacent nodes of the dequeued node
            for (int neighbor = 0; neighbor < numNodes; neighbor++) {
                // If the neighbor is connected and has not been visited, mark it as visited and enqueue it
                if (adjacencyMatrix[currentNode][neighbor] && !visited[neighbor]) {
                    visited[neighbor] = true;
                    bfsQueue.push(neighbor);
                }
            }
        }

        // Check if all nodes have been visited
        for (bool visit : visited) {
            if (!visit)
                return false;
        }
        return true;
    }

    void drawGraph() {
        sf::RenderWindow window(sf::VideoMode(800, 800), "Wireless Network Graph");

        while (window.isOpen()) {
            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            window.clear();

            // Draw edges
            for (int i = 0; i < numNodes; i++) {
                for (int j = i + 1; j < numNodes; j++) {
                    if (adjacencyMatrix[i][j]) {
                        sf::Vertex line[] = {
                            sf::Vertex(sf::Vector2f(nodes[i].x * 80, nodes[i].y * 80)),
                            sf::Vertex(sf::Vector2f(nodes[j].x * 80, nodes[j].y * 80))
                        };
                        window.draw(line, 2, sf::Lines);
                    }
                }
            }

            // Draw nodes
            for (const Node& node : nodes) {
                sf::CircleShape circle(5);
                circle.setPosition(node.x * 80 - 5, node.y * 80 - 5);
                circle.setFillColor(sf::Color::Red);
                window.draw(circle);
            }

            window.display();
        }
    }

};

int main() {
    //cout << "\n                                           EXPERIMENT # 1\n" << endl;
    //cout << "                                         Before Topology Control\n" << endl;
    //WirelessNetwork network[10] = { WirelessNetwork(500,10.0),WirelessNetwork(550,10.0),WirelessNetwork(600,10.0),WirelessNetwork(650,10.0),WirelessNetwork(700,10.0),WirelessNetwork(750,10.0),WirelessNetwork(800,10.0),WirelessNetwork(850,10.0),WirelessNetwork(900,10.0),WirelessNetwork(950,10.0) };
    //cout << setw(10) << "Network" << setw(15) << "Max Degree" << setw(15) << "Avg Degree" << endl;
    //int j = 500;
    //for (int i = 0; i < 10; i++) {
    //    network[i].degreeCalc();
    //    cout << setw(10) << j << setw(15) << network[i].get_maxdeg() << setw(15) << network[i].get_avgdeg() << endl;
    //    j += 50;
    //}
    //cout << "\n                                          After Topology Control\n" << endl;
    //cout << setw(10) << "Network" << setw(15) << "Max Degree" << setw(15) << "Avg Degree" << endl;
    //int k = 500;
    //for (int i = 0; i < 10; i++) {
    //    network[i].XTCTopologyControl();
    //    network[i].degreeCalc();
    //    cout << setw(10) << k << setw(15) << network[i].get_maxdeg() << setw(15) << network[i].get_avgdeg() << endl;
    //    k += 50;
    //}
    ////Experiment # 2
    //cout << "\n                                         EXPERIMENT # 2\n" << endl;
    //cout << "                                            Network G" << endl;
    //srand(time(nullptr));
    //WirelessNetwork network_G(1000, 10);
    //cout << left << setw(10) << "Source" << left << setw(12) << "Destination" << left << setw(18) << "Path Distance" << left << setw(12) << "Path Found" << left << setw(20) << "Path" << endl;
    //cout << string(102, '-') << endl;
    //for (int i = 0; i < 10; i++) {
    //    int src = rand() % 500;
    //    int dest = rand() % 500;
    //    vector<Node> path = network_G.compassRouting(src, dest);
    //    cout << left << setw(10) << src
    //        << left << setw(12) << dest;
    //    if (path.size() < 25) {
    //        cout << left << setw(18) << network_G.routeDist(path)
    //            << left << setw(12) << "Found";
    //        for (int i = 0; i < path.size(); i++) {
    //            cout << path[i].nodeNum;
    //            if ((i + 1) != path.size()) {
    //                cout << "->";
    //            }
    //        }
    //        cout << endl;
    //    }
    //    else {
    //        cout << left << setw(18) << network_G.routeDist(path)
    //            << left << setw(12) << "Not Found" << endl;
    //    }
    //    // cout << left << setw(20);
    //}
    ////Experiment 3
    //cout << "\n                                         EXPERIMENT # 3\n" << endl;
    //cout << "                                             Network H\n" << endl;
    //network_G.XTCTopologyControl();//topology control to get Network H
    ////WirelessNetwork network_H = network_G;
    //cout << left << setw(10) << "Source" << left << setw(12) << "Destination" << left << setw(18) << "Path Distance" << left << setw(12) << "Path Found" << left << setw(20) << "Path" << endl;
    //cout << string(102, '-') << endl;
    //for (int i = 0; i < 10; i++) {
    //    int src = rand() % 500;
    //    int dest = rand() % 500;
    //    vector<Node> path = network_G.compassRouting(src, dest);
    //    cout << left << setw(10) << src
    //        << left << setw(12) << dest;
    //    if (path.size() < 25) {
    //        cout << left << setw(18) << network_G.routeDist(path)
    //            << left << setw(12) << "Found";
    //        for (int i = 0; i < path.size(); i++) {
    //            cout << path[i].nodeNum;
    //            if ((i + 1) != path.size()) {
    //                cout << "->";
    //            }
    //        }
    //        cout << endl;
    //    }
    //    else {
    //        cout << left << setw(18) << network_G.routeDist(path)
    //            << left << setw(12) << "Not Found" << endl;
    //    }
    //}

    WirelessNetwork network(500, 10.0);
    cout << "Initial Adjacency List:" << endl;
    network.printAdjacencyList();
    vector<Node> path = network.compassRouting(1, 25);
    cout << "Path:";
    for (int i = 0; i < path.size(); i++) {
        cout << path[i].nodeNum;
        if ((i + 1) != path.size()) {
            cout << "->";
        }
    }
    cout << endl << "Path Distance: " << network.routeDist(path) << endl;
    if (network.isConnected()) {
        cout << "Connected" << endl;
    }
    else {
        cout << "Not connected" << endl;
    }

    network.drawGraph();

    network.XTCTopologyControl();

    cout << "Adjacency List after XTC Topology Control:" << endl;
    network.printAdjacencyList();

    path.clear();
    path = network.compassRouting(1, 25);
    cout << "Path:";
    for (int i = 0; i < path.size(); i++) {
        //cout << path[i].nodeNum;
        if ((i + 1) != path.size()) {
            //cout << "->";
        }
    }
    cout << endl << "Path Distance: " << network.routeDist(path) << endl;
    if (network.isConnected()) {
        cout << "Connected" << endl;
    }
    else {
        cout << "Not connected" << endl;
    }
    network.drawGraph();
    system("pause");
    return 0;
}