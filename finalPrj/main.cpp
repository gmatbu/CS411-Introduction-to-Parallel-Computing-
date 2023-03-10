#include <stdlib.h>
#include <unordered_map>
#include <vector>
#include <string.h>
#include <sys/time.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <iostream>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;
using std::pair;
using std::string;
using std::stringstream;
using std::unordered_map;
using std::vector;

typedef struct GraphNodes
{
    int degree;
    int *neighbors;
} GraphNode;

void printArrToFile(int, double *, ofstream&);
vector<string> split(string,char);
void getTopFive(int, double*, ofstream&);
int getMax(ifstream&);
void addEdge(unordered_map<int, GraphNode> &, int v, vector<int> &);
void buildAdjacencyList(ifstream&, unordered_map<int, GraphNode> &);
void printAdjacency(unordered_map<int, GraphNode> &);
double getRandom(int);
int getRandomInt(int, int, int);
void getPageranks(unordered_map<int, GraphNode>&, double*, int, int, double);
void getPageranksParallel(unordered_map<int, GraphNode>&, double*, int, int, double);

int main(int argc, char const *argv[])
{
    int p, max, walks;
    double damp, start_time, total_time;
    long double total;
    double * pageranks;
    string fname, out;
    stringstream stream;
    ifstream infile;
    ofstream outfile;
    unordered_map<int, GraphNode> adj;
    vector<string> splits;

    fname = argv[1];
    p = std::stoi(argv[2]);
    walks = std::stoi(argv[3]);
    damp = std::stod(argv[4]);

    splits = split(fname, '/');
    out = splits.back();
    splits = split(out, '.');
    out = splits.front();
    stream << out << "_";
    stream << p << "_"; 
    stream << walks << "_";
    stream << ((int) (damp*100.0));
    stream >> out;
    
    outfile.open((out + ".csv"));
    infile.open(fname);

    omp_set_num_threads(p);

    max = getMax(infile);
    adj.rehash(max);
    buildAdjacencyList(infile, adj);

    pageranks = ((double*) malloc(max * sizeof(double)));

    start_time = omp_get_wtime();
    getPageranksParallel(adj, pageranks, max, walks, damp);
    // getPageranks(adj, pageranks, max, walks, damp);
    total_time = omp_get_wtime() - start_time;

    total = ((long double) adj.size()) * ((long double) walks);
    omp_set_num_threads(8);

    # pragma omp parallel for 
    for (int i = 0; i < max; i++)
    {
        pageranks[i] = ((double) (((long double) pageranks[i]) / total));
    }

    getTopFive(max, pageranks, outfile);
    // printArrToFile(max, pageranks, outfile);

    cout << split(out,'_').front() << ",";
    cout << walks << ",";
    cout << damp << ",";
    cout << p << ",";
    cout << total_time;
    cout << endl;

    if (infile.is_open())
        infile.close();
    if (outfile.is_open())
        outfile.close();

    return 0;
}

// sourced from https://thispointer.com/how-to-split-a-string-in-c/
vector<string> split(string str, char delimeter)
{
    stringstream stream(str);
    string item;
	vector<string> splits;
    while (getline(stream, item, delimeter))
		splits.push_back(item);
	return splits;
}

void printArrToFile(int size, double* arr, ofstream  &outfile)
{
    for (int i = 0; i < size; i++)
        outfile << i << "," << arr[i] << endl;
}

void getTopFive(int size, double * pageranks, ofstream &outfile)
{
    int top[5] = {-1, -1, -1, -1, -1};
    for (int i = 0; i  < size; i++)
    {
        int j,next,temp;
        for (j = 0; j < 5; j++)
        {
            if (top[j] == -1 || pageranks[i] > pageranks[top[j]])
            {
                next = pageranks[i];
                top[j] = i;
                break;
            }
        }
        while (++j < 5 && top[j] != -1)
        {
            temp = pageranks[j];
            pageranks[j] = next;
            next = temp;
        }
    }
    outfile << "vertex,rank" << endl; 
    for (int i = 0; i < 5; i++)
    {
        outfile << top[i] << ",";
        outfile << pageranks[top[i]];
        outfile << endl;
    }
}

int getMax(ifstream& infile)
{
    infile.clear();
    infile.seekg(0, std::ios::beg); 
    int to, from, max, local_max;
    max = 0;
    vector<int> neighbors;
    for (string line; getline(infile, line);)
    {
        if (line.at(0) == '#')
            continue;
        stringstream stream(line);
        if (stream >> from >> to)
        {
            local_max = std::max(to, from);
            if (local_max > max)
                max = local_max;
        }
    }
    return max;
}

void printAdjacency(unordered_map<int, GraphNode> &adj)
{
    for (auto const &pair : adj)
    {
        std::cout << pair.first << ": ";
        for (int i = 0; i < pair.second.degree; i++)
            cout << pair.second.neighbors[i] << ", ";
        cout << endl;
    }
}

void addEdge(unordered_map<int, GraphNode> &adj, int v, vector<int> &neighbors)
{
    GraphNode node = {
        ((int)neighbors.size()),
        ((int *)malloc(sizeof(int *) * neighbors.size()))};
    copy(neighbors.begin(), neighbors.end(), node.neighbors);
    adj[v] = node;
    neighbors.clear();
}

void buildAdjacencyList(ifstream &infile, unordered_map<int, GraphNode> &adj)
{
    infile.clear();
    infile.seekg(0, std::ios::beg); 
    int to, from;
    int prev = -1;
    vector<int> neighbors;
    for (string line; getline(infile, line);)
    {
        if (line.at(0) == '#')
            continue;
        stringstream stream(line);
        if (stream >> from >> to)
        {
            if (from != prev && prev != -1)
                addEdge(adj, prev, neighbors);
            if (adj.find(to) == adj.end())
                adj[to] = {0, ((int *)malloc(sizeof(int)))};
            prev = from;
            neighbors.insert(neighbors.begin(), to);
        }
    }
    if (adj.find(prev) == adj.end() || adj.at(prev).degree == 0)
        addEdge(adj, prev, neighbors);
}

double getRandom(int seed)
{
    int result;
    double random;
    struct drand48_data randBuffer;
    struct timeval time;
    gettimeofday(&time, NULL);
    seed = seed * time.tv_usec;
    srand48_r(seed, &randBuffer);
    drand48_r(&randBuffer, &random);
    return random;
}

int getRandomInt(int seed, int lower, int upper)
{
    double random;
    int result;
    random = getRandom(seed);
    random *= (100.0 * ((float)upper));
    result = (((int)random) % (upper - lower)) + lower;
    result = (((int)random) % (upper - lower)) + lower;
    return result;
}

void getPageranks(unordered_map<int, GraphNode> &adj, double* pageranks, int max, int walk, double damp)
{
    int k, next, current;
    double rand;

    for (int i = 0; i < max; i++)
    {
        if (adj.find(i) == adj.end())
            continue;
        current = i;
        k = walk;
        while (--k >= 0)
        {
            pageranks[current] += 1;
            rand = getRandom(13);
            if (rand < damp)
            {
                next = getRandomInt(k * 137, 1, max) - 1;
                while (adj.find(next) == adj.end())
                    next = getRandomInt(k * 137, 1, max) - 1;
                current = next;
            }
            else
            {
                int next = getRandomInt(97, 0, adj.at(current).degree + 1);
                if (next < adj.at(current).degree)
                    current = adj.at(current).neighbors[next]; 
            }
        }
    }
}

void getPageranksParallel(unordered_map<int, GraphNode> &adj, double* pageranks, int max, int walk, double damp)
{
    int k, current, next;
    double rand;
    
    #pragma omp parallel for private(k, current, next, rand) 
    for (int i = 0; i <= max; i++)
    {
        if (adj.find(i) == adj.end())
            continue;
        current = i;
        k = walk;
        while (--k >= 0)
        {
            #pragma omp atomic
                pageranks[current] += 1;

            rand = getRandom(13);
            if (rand < damp)
            {
                next = getRandomInt(k * 137, 1, max) - 1;
                while (adj.find(next) == adj.end())
                    next = getRandomInt(k * 137, 1, max) - 1;
                current = next;
            }
            else
            {
                next = getRandomInt(97, 0, adj.at(current).degree + 1);
                if (next < adj.at(current).degree)
                    current = adj.at(current).neighbors[next];
            }
        }
    }
}