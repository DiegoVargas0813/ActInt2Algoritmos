#include <iostream>
#include <vector>
#include <limits.h>
#include <algorithm>
#include <queue>
#include <cmath>
using namespace std;

const int VAL_INFINITO = INT_MAX; //Valor "infinito"

// Algoritmo de Prim. Complejidad: O(V^2).
vector<pair<int, int>> primMST(const vector<vector<int>>& GRAFO, int N) {
    vector<int> pesoMin(N, VAL_INFINITO); //Vector de pesos minimos        
    vector<bool> enMST(N, false);    //Vector de nodos visitados
    vector<int> nodoAnterior(N, -1);     //Vector de padres
    
    //Nodo raiz
    pesoMin[0] = 0;               
    nodoAnterior[0] = -1;              

    for (int count = 0; count < N - 1; ++count) { //Recorremos todos los nodos. O(V)
        int u = -1;

        // Encontramos el nodo no visitado con el peso menor
        for (int i = 0; i < N; ++i) {
            if (!enMST[i] && (u == -1 || pesoMin[i] < pesoMin[u])) { //Si no está en el MST y el peso es menor o es el primer nodo
                u = i; //Actualizamos el nodo
            }
        }

        enMST[u] = true; //Visitado

        // Actualizamos el peso de los nodos adjacentes
        for (int v = 0; v < N; ++v) {
            if (GRAFO[u][v] && !enMST[v] && GRAFO[u][v] < pesoMin[v]) { //Si no está en el MST y el peso es menor
                pesoMin[v] = GRAFO[u][v];
                nodoAnterior[v] = u;
            }
        }
    }

    //Lista de arcos para el MST
    vector<pair<int, int>> rFinalMST;
    for (int i = 1; i < N; ++i) {
        rFinalMST.push_back({nodoAnterior[i], i});
    }
    return rFinalMST;
}

// Algoritmo Held-Karp para el TSP usando programación dinámica. Complejiad: O(2^N * N^2).
int tsp(int posNodo, int BitMaskVisitado, int N, const vector<vector<int>>& GRAFO, vector<vector<int>>& dpTabla) {
    if (BitMaskVisitado == (1 << N) - 1) {
        return GRAFO[posNodo][0]; 
    }
    if (dpTabla[posNodo][BitMaskVisitado] != -1) {
        return dpTabla[posNodo][BitMaskVisitado];
    }
    
    int costoMin = VAL_INFINITO;
    for (int city = 0; city < N; ++city) {
        if ((BitMaskVisitado & (1 << city)) == 0) { 
            int nuevoCosto = GRAFO[posNodo][city] + tsp(city, BitMaskVisitado | (1 << city), N, GRAFO, dpTabla);
            costoMin = min(costoMin, nuevoCosto);
        }
    }
    return dpTabla[posNodo][BitMaskVisitado] = costoMin;
}

//Reconstruir la ruta óptima
void reconstructRoute(int posNodo, int BitMaskVisitado, int N, const vector<vector<int>>& GRAFO, vector<vector<int>>& dpTabla, vector<int>& path) {
    if (BitMaskVisitado == (1 << N) - 1) {
        return;
    }

    int costoMin = VAL_INFINITO;
    int nextCity = -1;

    for (int city = 0; city < N; ++city) {
        if ((BitMaskVisitado & (1 << city)) == 0) {  
            int nuevoCosto = GRAFO[posNodo][city] + dpTabla[city][BitMaskVisitado | (1 << city)];
            if (nuevoCosto < costoMin) {
                costoMin = nuevoCosto;
                nextCity = city;
            }
        }
    }
    
    path.push_back(nextCity);
    reconstructRoute(nextCity, BitMaskVisitado | (1 << nextCity), N, GRAFO, dpTabla, path);
}

//Búsqueda en anchura para encontrar un camino de aumento. Complejidad: O(V^2).
bool bfs(const vector<vector<int>>& capacidad, vector<vector<int>>& flow, vector<int>& nodoAnterior, int N, int source, int sink) {
    vector<bool> BitMaskVisitado(N, false);
    queue<int> q;
    q.push(source);
    BitMaskVisitado[source] = true;

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v = 0; v < N; ++v) {
            if (!BitMaskVisitado[v] && capacidad[u][v] - flow[u][v] > 0) {
                q.push(v);
                BitMaskVisitado[v] = true;
                nodoAnterior[v] = u;
                if (v == sink) return true;
            }
        }
    }
    return false;
}

// Algoritmo de Ford-Fulkerson para calcular el flujo máximo. Complejidad: O(V * E^2).z
int fordFulkerson(const vector<vector<int>>& capacidad, int N, int source, int sink) {
    vector<vector<int>> flow(N, vector<int>(N, 0)); 
    vector<int> nodoAnterior(N);

    int maxFlujo = 0;

    // Mientras haya un camino de aumento
    while (bfs(capacidad, flow, nodoAnterior, N, source, sink)) {
        // Encuentra el flujo máximo posible en el camino de aumento
        int pathFlow = VAL_INFINITO; 
        for (int v = sink; v != source; v = nodoAnterior[v]) { //Recorremos el camino de aumento
            int u = nodoAnterior[v]; //Nodo anterior
            pathFlow = min(pathFlow, capacidad[u][v] - flow[u][v]); //Flujo máximo
        }

        // Actualiza los flujos a lo largo del camino de aumento
        for (int v = sink; v != source; v = nodoAnterior[v]) {
            int u = nodoAnterior[v];
            flow[u][v] += pathFlow;
            flow[v][u] -= pathFlow; 
        }

        maxFlujo += pathFlow;
    }
    return maxFlujo;
}

// distancia entre dos puntos
double calcularDistancia(double x1, double y1, double x2, double y2) {
    return sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
}

//central más cercana
int encontrarCentralMasCercana(double x, double y, const vector<pair<double, double>>& centrales) {
    int centralMasCercana = -1;
    double distanciaMinima = VAL_INFINITO;

    for (int i = 0; i < centrales.size(); ++i) {
        double dist = calcularDistancia(x, y, centrales[i].first, centrales[i].second);
        if (dist < distanciaMinima) {
            distanciaMinima = dist;
            centralMasCercana = i;
        }
    }
    return centralMasCercana;
}

int main() {
    int N;
    cout << "Ingresa el numero de colonias: ";
    cin >> N;

    vector<vector<int>> GRAFO(N, vector<int>(N));
    cout << "Ingresa la matriz de distancias entre las colonias:\n";
    for (int i = 0; i < N; ++i) {
        cout << "Colonia " << char('A' + i) << ":\n";

        for (int j = 0; j < N; ++j) {
            cout << "Distancia de " << char('A' + i) << " a " << char('A' + j) << ": ";
            cin >> GRAFO[i][j];
        }
    }

    //atriz de capacidades
    vector<vector<int>> capacidad(N, vector<int>(N));
    cout << "Ingresa la matriz de capacidades de transmisión de datos entre colonias:\n";
    for (int i = 0; i < N; ++i) {
        cout << "Colonia " << char('A' + i) << ":\n";
        for (int j = 0; j < N; ++j) {
            cout << "Capacidad de " << char('A' + i) << " a " << char('A' + j) << ": ";
            cin >> capacidad[i][j];
        }
    }


    //coordenadas de las centrales
    vector<pair<double, double>> centrales(N);
    cout << "Ingresa las coordenadas de las centrales (x, y):\n";
    for (int i = 0; i < N; ++i) {
        cout << "Central " << i + 1 << ": ";
        cin >> centrales[i].first >> centrales[i].second;
    }

    //MST usando el algoritmo de Prim
    vector<pair<int, int>> rFinalMST = primMST(GRAFO, N);

    //lista de arcos que forman el MST
    cout << "\nCableado optimo entre colonias:\n";
    for (const auto& edge : rFinalMST) {
        char u = 'A' + edge.first;
        char v = 'A' + edge.second;
        cout << "(" << u << ", " << v << ")\n";
    }

    //ruta óptima usando TSP
    vector<vector<int>> dpTabla(N, vector<int>(1 << N, -1));
    int minPathCost = tsp(0, 1, N, GRAFO, dpTabla);

    // Reconstruimos la ruta del TSP
    vector<int> path = {0}; 
    reconstructRoute(0, 1, N, GRAFO, dpTabla, path);

    // Mostramos la ruta óptima
    cout << "\nRuta optima para distribucion:\n";
    for (int i = 0; i < path.size(); ++i) {
        cout << char('A' + path[i]) << " -> ";
    }
    cout << "A" << endl; 
    cout << "Costo minimo de la ruta: " << minPathCost << endl;

    // Calculamos el flujo máximo entre la colonia 0 (inicio) y la colonia N-1 (final)
    int maxFlujo = fordFulkerson(capacidad, N, 0, N - 1);

    cout << "\nEl flujo maximo de informacion entre la colonia A (origen) y la colonia " 
         << char('A' + N - 1) << " (destino) es: " << maxFlujo << endl;

    // Calcular la cobertura de cada central
    cout << "\nPoligonos de cobertura:\n";
    for (int i = 0; i < N; ++i) {
        cout << "Central " << i + 1 << ": ";
        for (int j = 0; j < N; ++j) {
            int centralMasCercana = encontrarCentralMasCercana(centrales[j].first, centrales[j].second, centrales);
            cout << "(" << centrales[j].first << ", " << centrales[j].second << ")";
            if (j != N - 1) cout << ", ";
        }
        cout << "\n";
    }

    return 0;
}