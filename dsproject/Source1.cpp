//#include <iostream>
//#include <vector>
//#include <random>
//#include <cmath>
//#include <cstdlib> //random function // rand // srand
//#include <ctime>
//#include<algorithm>
//#include <iomanip> //set width
//#include <Windows.h>
//
//using namespace std;
//
//class Random {
//public:
//    Random() {
//        srand(time(nullptr));  //srand function , time seed, argument nullptr
//    }
//
//    double nextDouble() {
//        return static_cast<double>(rand()) / RAND_MAX; //floating point number btw zero and 1 generate
//    }
//};
//
//struct Node {
//    int nodeNum;   //vertex number
//    double x, y;   //2D
//};
//
//class WirelessNetwork {
//private:
//    int numNodes;    // no of vertices
//    double size;
//    int totaldeg, maxdeg;
//    float avgdeg;
//    vector<Node> nodes;
//    vector<vector<bool>> adjacencyMatrix;
//
//public:
//    WirelessNetwork(int numNodes, double size) : numNodes(numNodes), size(size) {
//        Random random;
//
//        // Generate random points within the square
//        for (int i = 0; i < numNodes; i++) {
//            Node node;
//            node.nodeNum = i;
//            node.x = random.nextDouble() * size;
//            node.y = random.nextDouble() * size;
//            nodes.push_back(node);
//        }
//
//        // Create the adjacency matrix
//        adjacencyMatrix.resize(numNodes, vector<bool>(numNodes, false));
//
//        // Create edges within distance 1 unit
//        for (int i = 0; i < numNodes; i++) {
//            for (int j = i + 1; j < numNodes; j++) {
//                if (calculateDistance(nodes[i], nodes[j]) <= 1.0) {
//                    adjacencyMatrix[i][j] = true;
//                    adjacencyMatrix[j][i] = true;
//                }
//            }
//        }
//    }
//
//    // Default constructor
//    WirelessNetwork() : WirelessNetwork(500, 10.0) {} // Using default parameters
//    double get_maxdeg() {
//        return maxdeg;
//    }
//    double get_avgdeg() {
//        return avgdeg;
//    }
//    void degreeCalc() {
//        int temp;
//        totaldeg = 0;
//        avgdeg = 0.0;
//        maxdeg = 0;
//        for (int i = 0; i < numNodes; i++) {
//            temp = 0;            // so that for one node we cal the max degree of the current node
//            for (int j = 0; j < numNodes; j++) {
//                if (adjacencyMatrix[i][j] == 1) {   // if edge exists then add 1 to degree
//                    totaldeg += 1;
//                    temp += 1;
//                }
//            }
//            if (temp > maxdeg) {
//                maxdeg = temp;
//            }
//        }
//        avgdeg = (float)totaldeg / numNodes;
//    }
//
//    double calculateDistance(const Node& node1, const Node& node2) {
//        double dx = node1.x - node2.x;
//        double dy = node1.y - node2.y;
//        return sqrt(dx * dx + dy * dy);  //distance formula underoot of  x^2 + y^2 
//    }
//
//    double calculateDistanceint(int s, int t) {
//        double dx = nodes[s].x - nodes[t].x;
//        double dy = nodes[s].y - nodes[t].y;
//        return sqrt(dx * dx + dy * dy);
//    }
//
//    vector<int> getNeighborsx(int nodeNum) {
//        vector<int> neighbors;
//        for (int i = 0; i < numNodes; i++) {
//            if (adjacencyMatrix[nodeNum][i]) {   // if edge exists between current node and node i
//                neighbors.push_back(i);
//            }
//        }
//        return neighbors;   //returns the vector that contains the neighbors of the current node
//    }
//
//    void printAdjacencyList() {
//        degreeCalc();
//        for (int i = 0; i < numNodes; i++) {
//            cout << "Node " << i << " connected to: ";
//            for (int j = 0; j < numNodes; j++) {
//                if (adjacencyMatrix[i][j]) {
//                    cout << j << " ";
//                }
//            }
//            cout << endl;
//        }
//        cout << "Total degree: " << totaldeg << endl;
//        cout << "Max degree: " << maxdeg << endl;
//        cout << "Avg Degree: " << avgdeg << endl;
//    }
//
//    void deleteEdge(int u, int v) {
//        adjacencyMatrix[u][v] = 0;
//        adjacencyMatrix[v][u] = 0;
//    }
//
//    vector<Node> compassRouting(int s, int t, int depthcontrol = 100) {   //s is source node, t is destination node
//        vector<Node> neighb = getNeighbors(s);   //current nodes neighbours
//        float angle, minangle = 90;      //stores angle between vectors
//        if (s == t) {                        //if source=dest
//            return { nodes[s] };
//        }
//        else {
//            int nextVertex = -1;
//
//            vector<double> v1;
//            v1.push_back(nodes[t].x - nodes[s].x);     //distance between source node and destination node
//            v1.push_back(nodes[t].y - nodes[s].y);
//
//            for (int i = 0; i < neighb.size(); i++) {    // loop iterates for total no of neighbours 
//                vector<double> v2 = {};
//                v2.push_back(neighb[i].x - nodes[s].x);        //distance between source node and neighbor
//                v2.push_back(neighb[i].y - nodes[s].y);
//
//                angle = angleBetweenVectors(v1, v2);
//
//                if (angle < minangle) {
//                    nextVertex = neighb[i].nodeNum;
//                    minangle = angle;
//                }
//            }
//            if (nextVertex == -1 || depthcontrol <= 0) {
//                return {};      //indicates failure in finding the path
//            }
//            vector<Node> path = compassRouting(nextVertex, t, depthcontrol - 1);  // new source node is next vertex and depth control is subtracted by 1
//            path.insert(path.begin(), nodes[s]);  //inserts node[s] at beiginning of the path
//            return path;
//        }
//    }
//
//    vector<Node> getNeighbors(int s) {
//        vector<Node> neighbors;
//        for (int i = 0; i < numNodes; i++) {
//            if (adjacencyMatrix[s][i]) {
//                neighbors.push_back(nodes[i]);
//            }
//        }
//        return neighbors;
//    }
//
//    double dotProduct(const vector<double>& vector1, const vector<double>& vector2) {
//        double dotProduct = 0.0;
//        for (int i = 0; i < vector1.size(); i++) {
//            dotProduct += vector1[i] * vector2[i];   //a.b 
//        }
//        return dotProduct;
//    }
//
//    double magnitude(const vector<double>& vector) {
//        double magnitude = 0.0;
//        for (double value : vector) {
//            magnitude += value * value;   // |a|
//        }
//        return sqrt(magnitude);
//    }
//
//    double angleBetweenVectors(const vector<double>& vector1, const vector<double>& vector2) {
//        double dot = dotProduct(vector1, vector2);
//        double mag1 = magnitude(vector1);
//        double mag2 = magnitude(vector2);
//
//        return acos(dot / (mag1 * mag2));    // cosine inverse of a.b/ |a||b|
//    }
//
//
//    void sortNeighbors(vector<int>& neighbors, const Node& currentNode) {    //selection sort
//        int n = neighbors.size();
//
//        for (int i = 0; i < n - 1; i++) {
//            int minIndex = i;
//
//            for (int j = i + 1; j < n; j++) {
//                double distA = calculateDistance(currentNode, nodes[neighbors[j]]);
//                double distB = calculateDistance(currentNode, nodes[neighbors[minIndex]]);
//
//                if (distA < distB) { //this means that the neighbour at index is closer to current node than the one at i
//                    minIndex = j;
//                }
//            }
//
//            if (minIndex != i) {
//                swap(neighbors[i], neighbors[minIndex]);
//            }
//        }
//    }
//
//    float routeDist(const vector<Node>& path) { // distance between the first and last node of the path
//        float dist = 0;
//        for (int i = 0; i < path.size() - 1; i++) {     // It stops at path.size() - 1 because the distance between
//            //the last node and its subsequent node does not need to be calculated (as there is no subsequent node).
//            Node current = path[i];
//            Node next = path[i + 1];
//            dist += calculateDistance(current, next);
//        }
//        return dist;
//    }
//
//
//    void XTCTopologyControl() {
//        vector<vector<bool>> updatedAdjacencyMatrix(numNodes, vector<bool>(numNodes, false));
//
//        for (int i = 0; i < numNodes; i++) {           //distance between current node and its neighbours
//            vector<pair<double, int>> neighborDistances;
//            for (int j = 0; j < numNodes; j++) {
//                if (i != j) {
//                    double distance = calculateDistance(nodes[i], nodes[j]);
//                    neighborDistances.push_back({ distance, j });
//                }
//            }
//
//            sort(neighborDistances.begin(), neighborDistances.end());   //the distances are sorted 
//
//            for (int j = 0; j < 6; j++) {                    //first six are added to the updated adjacency matrix
//                int neighbor = neighborDistances[j].second;
//                updatedAdjacencyMatrix[i][neighbor] = true;
//                updatedAdjacencyMatrix[neighbor][i] = true;
//            }
//        }
//
//        // Remove extra edges if any node has more than 6 neighbors
//        for (int i = 0; i < numNodes; i++) {
//            vector<int> neighbors = getNeighborsx(i);
//            if (neighbors.size() > 6) {
//                sortNeighbors(neighbors, nodes[i]);
//                for (int j = 6; j < neighbors.size(); j++) {
//                    int neighbor = neighbors[j];
//                    updatedAdjacencyMatrix[i][neighbor] = false;
//                    updatedAdjacencyMatrix[neighbor][i] = false;
//                }
//            }
//        }
//        //dfcvfvvd
//        adjacencyMatrix = updatedAdjacencyMatrix;
//
//        for (int i = 0; i < numNodes; i++) {   //reassuring that neighbours donot exceed six
//            int count = 0;
//            for (int j = 0; j < numNodes; j++) {
//                if (updatedAdjacencyMatrix[i][j] == true) {
//                    count++;
//                }
//            }
//            if (count > 6) {
//                vector<int> neighbors = getNeighborsx(i);
//                for (int j = 6; j < count; j++) {
//                    int neighbor = neighbors[j];
//                    updatedAdjacencyMatrix[i][neighbor] = false;
//                    updatedAdjacencyMatrix[neighbor][i] = false;
//                    /* cout << " " << i << " disconnected to " << neighbor << endl;*/
//                }
//
//            }
//            /* cout << i << " -> " << count << endl;*/
//        }
//        /*   cout << "aamfbjdbvjdsbvjhsvbjhv " << endl;*/
//        for (int i = 0; i < numNodes; i++) {   //output
//            int count1 = 0;
//            for (int j = 0; j < numNodes; j++) {
//                if (updatedAdjacencyMatrix[i][j] == true) {
//                    count1++;
//                }
//            }
//            /*  cout << i << " -> " << count1 << endl;*/
//        }
//
//        adjacencyMatrix = updatedAdjacencyMatrix;
//    }
//};
//
//int main() {
//    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);  // interacts with console
//    CONSOLE_SCREEN_BUFFER_INFO consoleInfo;   // stores info about current console buffer
//    WORD originalAttributes;    //info about original console attributes
//
//    // Get the current console attributes
//    GetConsoleScreenBufferInfo(hConsole, &consoleInfo);    //stores the current console info in the consoleInfo buffer into hconsole
//    originalAttributes = consoleInfo.wAttributes;
//
//    // Set the console text color to green
//   /* SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
//
//    cout << "\n                                           EXPERIMENT # 1\n" << endl;
//    cout << "                                         Before Topology Control\n" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    cout << string(102, '-') << endl;
//    WirelessNetwork network[10] = { WirelessNetwork(500,10.0),WirelessNetwork(550,10.0),WirelessNetwork(600,10.0),WirelessNetwork(650,10.0),WirelessNetwork(700,10.0),WirelessNetwork(750,10.0),WirelessNetwork(800,10.0),WirelessNetwork(850,10.0),WirelessNetwork(900,10.0),WirelessNetwork(950,10.0) };
//    /* for (int i = 0; i < 3; i++) {
//         network[i].printAdjacencyList();
//  }*/
//  /*   SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
//    cout << setw(10) << "Network" << setw(15) << "Max Degree" << setw(15) << "Avg Degree" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    int j = 500;
//    for (int i = 0; i < 10; i++) {
//        network[i].degreeCalc();
//        cout << setw(10) << j << setw(15) << network[i].get_maxdeg() << setw(15) << network[i].get_avgdeg() << endl;
//        j += 50;
//    }
//    /* SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
//
//    cout << "\n                                          After Topology Control\n" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    cout << string(102, '-') << endl;
//    /* SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
//
//    cout << setw(10) << "Network" << setw(15) << "Max Degree" << setw(15) << "Avg Degree" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    int k = 500;
//
//    for (int i = 0; i < 10; i++) {
//        network[i].XTCTopologyControl();
//        network[i].degreeCalc();
//        cout << setw(10) << k << setw(15) << network[i].get_maxdeg() << setw(15) << network[i].get_avgdeg() << endl;
//        k += 50;
//    }
//
//
//    //Experiment # 2
//  /*  SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
//
//    cout << "\n                                         EXPERIMENT # 2\n" << endl;
//    cout << "                                            Network G" << endl;
//    srand(time(nullptr));
//    WirelessNetwork network_G(1000, 10);
//    /* SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
//    cout << left << setw(10) << "Source" << left << setw(12) << "Destination" << left << setw(18) << "Path Distance" << left << setw(12) << "Path Found" << left << setw(20) << "Path" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    cout << string(102, '-') << endl;
//
//    for (int i = 0; i < 10; i++) {
//        int src = rand() % 1000;
//        int dest = rand() % 1000;
//
//        vector<Node> path = network_G.compassRouting(src, dest);
//
//        cout << left << setw(10) << src
//            << left << setw(12) << dest;
//
//        if (path.size() < 25) {
//            cout << left << setw(18) << network_G.routeDist(path)
//                << left << setw(12) << "Found";
//            for (int i = 0; i < path.size(); i++) {
//                cout << path[i].nodeNum;
//                if ((i + 1) != path.size()) {
//                    cout << "->";
//                }
//            }
//            cout << endl;
//        }
//        else {
//            cout << left << setw(18) << network_G.routeDist(path)
//                << left << setw(12) << "Not Found" << endl;
//        }
//        // cout << left << setw(20);
//
//    }
//    //Experiment 3
//  /*  SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_BLUE);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN | FOREGROUND_BLUE | FOREGROUND_INTENSITY);
//
//    cout << "\n                                         EXPERIMENT # 3\n" << endl;
//    cout << "                                             Network H\n" << endl;
//    network_G.XTCTopologyControl();//topology control to get Network H
//    //WirelessNetwork network_H = network_G;
//   /* SetConsoleTextAttribute(hConsole, FOREGROUND_GREEN);*/
//    SetConsoleTextAttribute(hConsole, FOREGROUND_RED | FOREGROUND_GREEN | FOREGROUND_INTENSITY);
//    cout << left << setw(10) << "Source" << left << setw(12) << "Destination" << left << setw(18) << "Path Distance" << left << setw(12) << "Path Found" << left << setw(20) << "Path" << endl;
//    SetConsoleTextAttribute(hConsole, originalAttributes);
//    cout << string(102, '-') << endl;
//    for (int i = 0; i < 10; i++) {
//        int src = rand() % 1000;
//        int dest = rand() % 1000;
//
//        vector<Node> path = network_G.compassRouting(src, dest);
//
//        cout << left << setw(10) << src
//            << left << setw(12) << dest;
//
//        if (path.size() < 25) {
//            cout << left << setw(18) << network_G.routeDist(path)
//                << left << setw(12) << "Found";
//            for (int i = 0; i < path.size(); i++) {
//                cout << path[i].nodeNum;
//                if ((i + 1) != path.size()) {
//                    cout << "->";
//                }
//            }
//            cout << endl;
//        }
//        else {
//            cout << left << setw(18) << network_G.routeDist(path)
//                << left << setw(12) << "Not Found" << endl;
//        }
//    }
//    system("pause");
//    return 0;
//}