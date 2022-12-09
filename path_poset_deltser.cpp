// Author: Jason Smith
// Department of Physics and Mathematics, Nottingham Trent University
// Date: Feb 2022

//Run with:
//./path_poset num_vertices num_threads edges
//where edges is the address of a text file containing a list of edges
//For example:
//./path_poset 100 8 test.edges

#include <vector>
#include <thread>
#include <deque>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <algorithm>
#include <functional>
#include <sstream>
#include <fstream>
#include <string>
#include <math.h>
#include <numeric>
#include <iterator>
#include <boost/functional/hash.hpp>

typedef uint32_t vertex_t;
typedef uint32_t index_t;
typedef std::pair<vertex_t,vertex_t> edge_t;
typedef std::pair<std::vector<std::pair<int,int>>,std::vector<std::pair<int,int>>> relation_t;
typedef std::unordered_map<edge_t,index_t, boost::hash<edge_t>> dict_t;
typedef std::vector<std::vector<std::vector<std::unordered_set<index_t>>>> paths_t;

//****************************************************************************//
//Functions for parsing input file

inline std::string trim(const std::string& s) {
    auto wsfront = std::find_if_not(s.begin(), s.end(), [](int c) { return std::isspace(c); });
    auto wsback = std::find_if_not(s.rbegin(), s.rend(), [](int c) { return std::isspace(c); }).base();
    return (wsback <= wsfront ? std::string() : std::string(wsfront, wsback));
}
unsigned int string_to_uint(std::string s) { return atoi(s.c_str()); }
std::vector<vertex_t> split(const std::string& s, char delim, const std::function<vertex_t(std::string)>& transform) {
    std::vector<vertex_t> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) elems.push_back(transform(item));
    return elems;
}

std::vector<edge_t> get_edges(vertex_t num_vertices, std::string edge_address) {
    std::vector<edge_t> edges;
    std::string line;
    std::ifstream input_stream(edge_address);
    if (input_stream.fail()) {  std::cerr << "couldn't open file " << edge_address << std::endl; exit(-1); }
    while (not input_stream.eof()) {
        std::getline(input_stream, line);
        line = trim(line);
        if (line.length() == 0) continue;
        std::vector<vertex_t> vertices = split(line, ' ', string_to_uint);
        if (vertices[0] >= num_vertices || vertices[1] >= num_vertices){
            std::cerr << "ERROR: edge " << line << " not allowed as number of vertices = " << num_vertices << std::endl;
            exit(-1);
        }
        edges.push_back(edge_t(vertices[0],vertices[1]));
    }
    return edges;
}

//****************************************************************************//
//Graph class

struct graph_t{
    vertex_t num_vertices;
    std::vector<std::unordered_set<vertex_t>> out_neighbours;

    //Constructor
    graph_t(vertex_t _number_of_vertices, std::vector<edge_t> edges)
        : num_vertices(_number_of_vertices) {
        out_neighbours.assign(_number_of_vertices, std::unordered_set<vertex_t>());
        for(auto edge : edges){
            out_neighbours[edge.first].insert(edge.second);
        }
    }
};

bool is_connected_by_an_edge(vertex_t from, vertex_t to, const graph_t& graph) {
    return graph.out_neighbours[from].find(to) != graph.out_neighbours[from].end();
}

//****************************************************************************//
//Subpath class

struct subpath_t {
    //A structure that stores a subgraph made of a disconnected set of directed paths
    //paths[i] is the i'th path stored as a vector of the vertices in that path
    //the first element of paths should be in increasing order to avoid double counting
    subpath_t* pre_path;
    bool first_subpath;
    std::vector<vertex_t> path;
    uint32_t size;

    //constructor
    subpath_t(vertex_t v, vertex_t w) {
        first_subpath = true;
        path.push_back(v);
        path.push_back(w);
        pre_path = nullptr;
        size = 1;
    }
    subpath_t(subpath_t* p, vertex_t v, vertex_t w){
        pre_path = p;
        first_subpath = false;
        path.push_back(v);
        path.push_back(w);
        size = p->size+1;
    }

    void add_next(vertex_t v){
        path.push_back(v);
        size += 1;
    }
    void print(bool first){
          for (int i = 0; i < path.size(); i++){
              std::cout << path[i] << " ";
          }
         std::cout << ": ";
         if(!first_subpath) pre_path->print(false);
         if (first) std::cout << std::endl;
    }
    vertex_t latest_initial(){
        return path[0];
    }
    vertex_t last(){
        return path.back();
    }
    void remove_last(){
        path.pop_back();
        size -= 1;
    }
    bool visited(vertex_t v){
        if (std::find(path.begin(), path.end(), v) != path.end()) return true;
        if (first_subpath) return false;
        return pre_path->visited(v);
    }
    void add_path_to(std::vector<std::vector<vertex_t>>& addto){
        addto.push_back(std::vector<vertex_t>(path));
        if(!first_subpath) pre_path->add_path_to(addto);
    }
};

//****************************************************************************//
std::unordered_set<index_t> frozenpath(subpath_t& to_copy, dict_t& edge_dict){
    std::vector<std::vector<vertex_t>> path_full;
    std::unordered_set<index_t> path;
    to_copy.add_path_to(path_full);
    for(auto p : path_full){
        for(int i = 0; i < p.size()-1; i++){
            path.insert(edge_dict[std::make_pair(p[i],p[i+1])]);
        }
    }
    return path;
}

//****************************************************************************//
//Functions for constructing paths and threading
void cont_path(subpath_t& current_subpath, const graph_t& graph, int64_t& s, int thread,
    paths_t& the_paths, dict_t& edge_dict);

void new_path(subpath_t& current_subpath, const graph_t& graph, int64_t& s, int thread,
    paths_t& the_paths, dict_t& edge_dict){
    for (int i = current_subpath.latest_initial()+1; i < graph.num_vertices; i++){
        if (!current_subpath.visited(i)){
            for (auto v :  graph.out_neighbours[i]){
                if (!current_subpath.visited(v)){
                    subpath_t new_subpath = subpath_t(&current_subpath,i,v);
                    the_paths[thread][new_subpath.size].push_back(frozenpath(new_subpath,edge_dict));
                    s += pow(-1,new_subpath.size);
                    cont_path(new_subpath, graph, s, thread, the_paths, edge_dict);
               }
            }
        }
    }
}

void cont_path(subpath_t& current_subpath, const graph_t& graph, int64_t& s, int thread,
    paths_t& the_paths, dict_t& edge_dict){
    //get neighbours of current end of path
    new_path(current_subpath, graph, s, thread, the_paths, edge_dict);
    const vertex_t last = current_subpath.last();
    for (auto const & v : graph.out_neighbours[last]){
        if (!current_subpath.visited(v)){
            current_subpath.add_next(v);
            s += pow(-1,current_subpath.size);
            the_paths[thread][current_subpath.size].push_back(frozenpath(current_subpath,edge_dict));
            cont_path(current_subpath, graph, s, thread, the_paths, edge_dict);
            current_subpath.remove_last();
       }
    }
}

void worker_thread(std::vector<edge_t>& edges, const graph_t& graph, int64_t& s, int thread, int num_threads,
    paths_t& the_paths, dict_t& edge_dict) {
    //For every edge starting at a vertex of start_vertices, create a path and call cont_path on it
    for (uint64_t i = thread; i < edges.size(); i += num_threads){
        subpath_t current_subpath = subpath_t(edges[i].first, edges[i].second);
        the_paths[thread][current_subpath.size].push_back(frozenpath(current_subpath, edge_dict));
        s -= 1;
        cont_path(current_subpath, graph, s, thread, the_paths, edge_dict);
    }
}

void combine_paths(paths_t& the_paths){
    for (int i = 1; i < the_paths.size(); i++){
        for (int j = 0; j < the_paths[0].size(); j++){
            the_paths[0][j].insert(std::end(the_paths[0][j]),std::begin(the_paths[i][j]),std::end(the_paths[i][j]));
            the_paths[i][j].clear();
        }
    }
}


std::vector<index_t> path_to_faces(std::unordered_set<index_t>& fpath, std::vector<std::unordered_set<index_t>>& dim_below_list){
    std::vector<index_t> face_list;
    for(auto i : fpath){
        std::unordered_set<index_t> one_less_path(fpath);
        one_less_path.erase(i);
        face_list.push_back(std::find(dim_below_list.begin(), dim_below_list.end(), one_less_path)-dim_below_list.begin());
    }
    return face_list;
}

//****************************************************************************//
//Main function

int main(int argc, char** argv) {
    if (argc < 5) {std::cerr << "Error: Missing input arguments" << std::endl; exit(-1);}
    int num_vertices = atoi(argv[1]);
    int num_threads = atoi(argv[2]);
    std::string edge_address = argv[3];
    std::string out_address = argv[4];

    std::vector<edge_t> edges = get_edges(num_vertices, edge_address);
    // Create a dictionary that maps an edge to location in list
    dict_t edge_dict;
    for (int i = 0; i < edges.size(); i++) edge_dict[edges[i]] = i;

    const graph_t graph(num_vertices, edges);

    std::cout << "Graph Loaded." << std::endl;

    paths_t the_paths;
    the_paths.resize(num_threads);
    for (int i = 0; i < num_threads; i++){ the_paths[i].resize(num_vertices); }

    //run the threads
    std::vector<int64_t> mf(num_threads,0);
    std::vector<std::thread> t(num_threads);
    for (size_t index = 0; index < num_threads - 1; ++index){
        t[index] = std::thread(worker_thread, std::ref(edges), std::ref(graph), std::ref(mf[index]), index, num_threads, std::ref(the_paths), std::ref(edge_dict));
    }

    worker_thread(edges, graph, mf[num_threads-1], num_threads-1, num_threads, the_paths, edge_dict);
    for (size_t i = 0; i < num_threads - 1; ++i) t[i].join();
    std::cout<<"Built all multipaths" <<std::endl;
    combine_paths(the_paths);
    std::cout<<"Combined thread output" <<std::endl;
    int max_dim = 1;
    for (int i = 0; i < the_paths[0].size(); i ++){
        if (the_paths[0][i].size() > 0) max_dim = i;
    }

    //TODO: parallelise this step
    std::vector<std::vector<std::vector<index_t>>> list_by_face;
    list_by_face.resize(max_dim);
    list_by_face[0].resize(edges.size());
    for (int i = 0; i < edges.size(); i++)  list_by_face[0][i] = std::vector<index_t>(1,i);
    for (int i = 1; i < list_by_face.size(); i++){
        list_by_face[i].resize(the_paths[0][i+1].size());
        for (int p = 0; p < the_paths[0][i+1].size(); p++){
            list_by_face[i][p] = path_to_faces(the_paths[0][i+1][p],the_paths[0][i]);;
        }
    }

    std::cout<<"Converted to deltser format" <<std::endl;

    std::ofstream outfile(out_address);
    outfile << "dim 0" << std::endl;
    for (int i = 0; i < edges.size()-1; i++) outfile << "0 ";
    outfile << "0" << std::endl;
    for (int i = 1; i < list_by_face.size(); i++){
        outfile << "dim " << i << std::endl;
        for (auto p : list_by_face[i]){
            for (auto j : p) outfile << j << " ";
            outfile << "0" << std::endl;
        }
    }
    // outfile << "dim 0" << std::endl;
    // for (int i = 0; i < edges.size()-1; i++) outfile << "0 ";
    // outfile << "0" << std::endl;
    // for (int i = 2; i < max_dim+1; i++){
    //     outfile << "dim " << i-1 << std::endl;
    //     for (auto p : the_paths[0][i]){
    //         for (auto j : p) outfile << j << " ";
    //         outfile << "0" << std::endl;
    //     }
    // }
    outfile.close();

    std::cout<<"Output file created." <<std::endl;
}
