#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include <chrono>
#include <omp.h>
#include <math.h>
#define N 3618 //number of nodes, predefined at compilation time
#define MC 50 //times to run the Monte Carlo simulation

std::map<int, std::vector<int> > getNetwork(std::string filename);
void printAdjacencyList(std::map<int, std::vector<int> > G);
void infectInitialNodes(bool infected[N], float initialRho);
float infectionProcess(std::map<int, std::vector<int> > G, bool infected[N], float recoveryNumber, int totalN, float beta);
float getInfectedFraction(bool infected[N]);
float mean(float arr[MC]);
float stdev(float arr[MC]);
float MMCA(unsigned adMatrix[N][N], float initialRho, float beta, int time, float mu);
float q(unsigned adMatrix[N][N], float p[N], int i, float beta);
void loadAdjacencyMatrix(unsigned adMatrix[N][N], std::string filename);

float mean(float arr[MC]){
    float meanValue = 0; int i;
    for (i = 0; i < MC; i++){
        meanValue += arr[i];
    }
    meanValue /= MC;
    return meanValue;
}

float stdev(float arr[MC]){
    float meanValue = mean(arr);
    float stdValue = 0; int i;
    for (i = 0; i < MC; i++){
        stdValue += (arr[i]-meanValue)*(arr[i]-meanValue);
    }
    stdValue /= MC;
    return sqrt(stdValue);
}

void printAdjacencyList(std::map<int, std::vector<int> > graph){
    int i, j;
    for(i = 0; i < graph.size(); i++){
        std::cout << i << ": ";
        for (j = 0; j < graph[i].size(); j++){
            std::cout << graph[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

void loadAdjacencyMatrix(unsigned adMatrix[N][N], std::string filename)
{
    std::ifstream f;
    f.open(filename);
    std::string s;
    int i, j;
    for(i = 0; i < N; i++)
    {
        std::getline(f, s);
        for(j = 0; j < N; j++)
        {
            adMatrix[i][j] = s[j] - 48; //Read from ASCII
        }
    }

    f.close();
}

std::map<int, std::vector<int> > getNetwork(std::string filename){
    std::ifstream f; //Declaration of input file
    std::string s; //Auxiliar string to read the file content
    f.open(filename); //Open the network file

    std::getline(f, s); //Get the first line, i.e., number of nodes

    std::getline(f, s); //Get the second line, i.e., number of edges
    int E = std::stoi(s); // We get the number of edges in variable E

    std::map<int, std::vector<int> > network;
    //We define the network as a map object, a pair of key-value
    //where the key is a integer, the node number, and the value is a vector
    //containing all the edges of the node. In other words, we are
    //creating an adjacency list with a map where the keys are the nodes
    //and the values are the neighbours

    int i, j;
    for (i = 0; i < N; i++)
    {
         network.insert(std::pair<int,std::vector<int> >(i, std::vector<int>()));
         // We insert as many elements as there are nodes, nodes being the key
         // and an empty vector as value, to push back the edges of the node
    }

    std::string aux;;
    int firstNode, secondNode; //Nodes composing a single edge of the network
    for (i = 0; i < E; i++)
    {
        std::getline(f, s);
        aux = "";
        for (j = 0; j < s.length(); j++)
        {
            if (s[j] != ' ')
            {
                aux += s[j];;
            }
            if (s[j] == ' ')
            {
                firstNode = std::stoi(aux); //Transform to int the first node
                aux = "";
            }
        }
        secondNode = std::stoi(aux); //Then the second node composing the edge

        //We append each node to the vector of the other node
        network[firstNode].push_back(secondNode);
        network[secondNode].push_back(firstNode);
    }
    f.close(); //Close the network file

    return network;
}

void infectInitialNodes(bool infected[N], float initialRho)
{
    //function to infect an initial number of nodes based on probability initialRho
    int i; float r;
    for (i = 0; i < N; i++)
    {
        r = (float) rand()/RAND_MAX; //rng
        if (r < initialRho) //if the random number is lesser, it gets infected
        {
            infected[i] = true;
        }
        else
        {
            infected[i] = false;
        }
    }
}

float getInfectedFraction(bool infected[N]) //get fraction of infection
{
    int i; float rho = 0;
    for (i = 0; i < N; i++){
        if (infected[i] == true)
        {
            rho += 1; //sum all the infected nodes
        }
    }
    rho = rho/N; //get fraction
    return rho;
}

float infectionProcess(std::map<int, std::vector<int> > graph,
    bool infected[N], float recoveryNumber, int totalN, float beta)
{
    int t, i, j; //loop indexes
    float r, rhoAverage = 0; //random rumber
    float copy[N];

    //First part of the iteration: change state of infected with RNG
    for(t = 0; t < totalN; t++) //loop for discrete time steps
    {
        for(i = 0; i < N; i++)
        {
            if (infected[i] == true)
            {
                r = (float) rand()/RAND_MAX;
                if (r < recoveryNumber)
                {
                    copy[i] = false;
                }
                else
                {
                    copy[i] = true;
                }
            }
            else
            {
                copy[i] = false;
            }
        }
        //Second part is to iterate through all the nodes, checking if
        //they have any infected neighbours, and then generating a random number
        //to infect the node
        for (i = 0; i < N; i++)
        {
            if (infected[i] == false) //Check if the node can be infected
            {
                for (j = 0; j < graph[i].size(); j++) //iterate through all its neighbours
                {
                    if (infected[graph[i][j]] == true) //check if neighbour is infected
                    {
                        r = (float) rand()/RAND_MAX; //rng
                        if (r < beta) //check if the neighbours will infect the node
                        {
                            copy[i] = true; //infection takes place
                            break; //we stop iterating through the other neighbours
                        }
                    }
                }
            }
        }

        //update infected array
        for (i = 0; i < N; i++)
        {
            infected[i] = copy[i];
        }

        if (t >= 500) //we compute the infected fraction only for t >= 500
        {
            rhoAverage += getInfectedFraction(infected);
        }
    }
    rhoAverage /= 500; //we get the average over the stable phase (we consider half)
    return rhoAverage;
}

float q(unsigned adMatrix[N][N], float p[N], int i, float beta)
{
    int j; float product = 1.0;
    for(j = 0; j < N; j++)
    {
        product *= (1-p[j]*beta*adMatrix[i][j]); //adMatrix[i][j] = r_ij
    }
    return product;
}

float MMCA(unsigned adMatrix[N][N], float initialRho, float beta, int time, float mu)
{
    int i, t;
    float p[N], copy[N];
    for (i = 0; i < N; i++)
    {
        p[i] = initialRho;
        copy[i] = initialRho;
    }

    #pragma omp parallel for private(t, i)
    for (t = 0; t < time; t++)
    {
        for (i = 0; i < N; i++)
        {
            //We apply the time evolution formula for the probability of infection
            p[i] = (1-q(adMatrix, copy, i, beta))*(1-copy[i])+(1-mu)*copy[i];
        }
        for (i = 0; i < N; i++)
        {
            copy[i] = p[i];
        }
    }

    float average = 0; //We compute the average to return
    for (i = 0; i < N; i++)
    {
        average += p[i];
    }
    return average/N;
}

int main ()
{
    auto start = std::chrono::high_resolution_clock::now();

    //Parameters of the model
    float initialRho = 0.05; //initial fraction of infected nodes
    int totalN = 1000; //Total number of simulation steps
    float recoveryNumber = .9; //fraction of nodes that recover per time step


    std::map<int, std::vector<int> > graph;
    //The network or graph G we will work with
    //We have chosen std::map as the data structure to contain the information
    //of the network in the form of an adjacency list

    std::string filename = "airports_UW.txt"; //The path to the file to open
    graph = getNetwork(filename); //Load the content of the network in graph
    //printAdjacencyList(graph); //We can print the adjacency list to check the graph

    bool infected[N]; //array containing information about the state of a node
    static unsigned adMatrix[N][N]; //adjacency matrix for MMCA computing
    loadAdjacencyMatrix(adMatrix, "airports_UWAM.txt");

    srand(time(NULL)); //Seed for random number generation (RNG)

    //Monte carlo simulations!
    //We run MC repetitions of the process of network infection for each betaValue
    //Then we export the results to a file for outside plotting

    float betaValue = 0.01; //We start off with a value of beta = 0.01
    float rhoArray[MC], MMCAresult;

    std::ofstream outputFile; //We declare the output file
    outputFile.open("airports_UW_mu9.txt"); //We open the file to write in
    int i;

    while (betaValue <= 1.0) //Beta = 1 will be the last value
    {

        //parallelization of the monte carlo simulations
        #pragma omp parallel for private(i, infected)
        for(i = 0; i < MC; i++) //We run MC iterations of the same conditions
        {
            infectInitialNodes(infected, initialRho); //infection of initial nodes
            rhoArray[i] = infectionProcess(graph, infected, recoveryNumber, totalN, betaValue);
            //simulation takes place
        }
        MMCAresult = MMCA(adMatrix, initialRho, betaValue, totalN, recoveryNumber);
        outputFile << betaValue << " " << mean(rhoArray) << " " << stdev(rhoArray) << " " << MMCAresult << "\n";
        std::cout << betaValue << " " << mean(rhoArray) << " " << stdev(rhoArray) << " " << MMCAresult << "\n";

        //We print beta and rho into separate columns
        betaValue += 0.01; //increase the value of beta and repeat
    }

    outputFile.close();

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Execution time: " << duration.count() << " miliseconds" << std::endl;

    return 0;
}
