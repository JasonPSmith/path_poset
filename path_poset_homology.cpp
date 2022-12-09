// Author: Jason Smith
// Department of Physics and Mathematics, Nottingham Trent University
// Date: Feb 2022

//Run with:
//./path_poset num_vertices num_threads edges
//where edges is the address of a text file containing a list of edges
//For example:
//./path_poset 100 8 test.edges

#define SORT_COLUMNS_BY_PIVOT

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
#include <cassert>
#include <queue>


typedef float value_t;
typedef int16_t coefficient_t;
typedef int64_t index_t;
typedef index_t entry_t;
typedef uint32_t vertex_t;
typedef std::pair<value_t, index_t> filtration_index_t;
typedef std::pair<vertex_t,vertex_t> edge_t;
typedef std::pair<std::vector<std::pair<int,int>>,std::vector<std::pair<int,int>>> relation_t;
typedef std::unordered_map<edge_t,index_t, boost::hash<edge_t>> dict_t;
typedef std::vector<std::vector<std::vector<std::unordered_set<index_t>>>> paths_t;

//-------------------------------------------------------------------------//
//BEGIN deltser

typedef std::deque<index_t> pivot_column_index_t;
const index_t INVALID_INDEX = std::numeric_limits<index_t>::max();

float string_to_float(std::string s) { return atof(s.c_str()); }

const index_t get_index(const entry_t& i) { return i; }
index_t get_coefficient(const entry_t& i) { return 1; }
entry_t make_entry(index_t _index, coefficient_t _value) { return entry_t(_index); }

const entry_t& get_entry(const entry_t& e) { return e; }

value_t get_filtration(const filtration_index_t& i) { return i.first; }
index_t get_index(const filtration_index_t& i) { return i.second; }

class filtration_entry_t : public std::pair<value_t, entry_t> {
public:
	filtration_entry_t() {}
	filtration_entry_t(const entry_t& e) : std::pair<value_t, entry_t>(0, e) {}
    filtration_entry_t(index_t _filtration, index_t _index)
        : std::pair<value_t, entry_t>(_filtration, entry_t(_index)) {}
	filtration_entry_t(value_t _filtration, index_t _index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(_filtration, make_entry(_index, _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index, coefficient_t _coefficient)
	    : std::pair<value_t, entry_t>(get_filtration(_filtration_index),
	                                 make_entry(get_index(_filtration_index), _coefficient)) {}
	filtration_entry_t(const filtration_index_t& _filtration_index)
	    : filtration_entry_t(_filtration_index, 1) {}
};

template <typename Heap> filtration_entry_t pop_pivot(Heap& column) {
    if (column.empty())
        return filtration_entry_t(-1);
    else {
        auto pivot = column.top();
        column.pop();
        while (!column.empty() &&
               get_index(column.top()) == get_index(pivot)) {
            column.pop();

            if (column.empty())
                return filtration_entry_t(-1);
            else {
                pivot = column.top();
                column.pop();
            }
        }
        return pivot;
    }
}

const entry_t& get_entry(const filtration_entry_t& p) { return p.second; }
const index_t get_index(const filtration_entry_t& p) { return get_index(get_entry(p)); }
const coefficient_t get_coefficient(const filtration_entry_t& p) { return get_coefficient(get_entry(p)); }
const value_t& get_filtration(const filtration_entry_t& p) { return p.first; }

struct greater_filtration_or_smaller_index {
	bool operator()(const filtration_index_t a, const filtration_index_t b) {
		return (get_filtration(a) > get_filtration(b)) ||
		       ((get_filtration(a) == get_filtration(b)) && (get_index(a) < get_index(b)));
	}
};

template <typename Entry> struct smaller_index {
	bool operator()(const Entry& a, const Entry& b) { return get_index(a) < get_index(b); }
};

class filtered_union_find {
	std::vector<index_t> parent;
	std::vector<std::vector<index_t>> rank;
	const std::vector<value_t> filtration;

public:
	filtered_union_find(const std::vector<value_t>& _filtration)
	    : rank(_filtration.size()), filtration(_filtration), parent(_filtration.size()) {
		for (index_t i = 0; i < _filtration.size(); ++i){
			parent[i] = i;
			rank[i] = std::vector<index_t>{i};
		}
	}
	index_t find(index_t x) {
		return parent[x];
	}
	value_t link(index_t x, index_t y) {
		x = find(x);
		y = find(y);
		if (x == y) return -1;
		if (filtration[x] < filtration[y] || (filtration[x] == filtration[y] && rank[x].size() > rank[y].size())){
			for(auto i : rank[y]) parent[i] = x;
			rank[x].insert(rank[x].end(),rank[y].begin(),rank[y].end());
			return filtration[y];
		} else {
			for(auto i : rank[x]) parent[i] = y;
			rank[y].insert(rank[y].end(),rank[x].begin(),rank[x].end());
			return filtration[x];
		}
	}
};

template <typename ValueType> class compressed_sparse_matrix {
	std::deque<size_t> bounds;
	std::deque<ValueType> entries;

public:
	size_t size() const { return bounds.size(); }

	void clear() {
		bounds.clear();
		bounds.shrink_to_fit();
		entries.clear();
		entries.shrink_to_fit();
	}

	typename std::deque<ValueType>::const_iterator cbegin(size_t index) const {
		assert(index < size());
		return index == 0 ? entries.cbegin() : entries.cbegin() + bounds[index - 1];
	}

	typename std::deque<ValueType>::const_iterator cend(size_t index) const {
		assert(index < size());
		return entries.cbegin() + bounds[index];
	}

	template <typename Iterator> void append_column(Iterator begin, Iterator end) {
		for (Iterator it = begin; it != end; ++it) { entries.push_back(*it); }
		bounds.push_back(entries.size());
	}

	void append_column() { bounds.push_back(entries.size()); }

	void push_back(ValueType e) {
		assert(0 < size());
		entries.push_back(e);
		++bounds.back();
	}

	void pop_back() {
		assert(0 < size());
		entries.pop_back();
		--bounds.back();
	}

	template <typename Collection> void append_column(const Collection collection) {
		append_column(collection.cbegin(), collection.cend());
	}
};

//END header
//-------------------------------------------------------------------------//
//BEGIN delta_complex




class delta_complex_t;

template <typename t>
std::vector<t> split(const std::string& s, char delim, const std::function<t(std::string)>& transform) {
	std::vector<t> elems;
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) elems.push_back(transform(item));
	return elems;
}

class delta_complex_cell_t {
public:
    std::vector<delta_complex_cell_t*> boundary;
    std::vector<index_t> coboundary;
	std::unordered_set<delta_complex_cell_t*> children;
    index_t location;
    value_t filtration;
    std::vector<index_t> vertices;
    int dimension;
	index_t oldest_coface;

    //initialise class
    delta_complex_cell_t(int _dim, index_t _v, value_t _filt) : dimension(_dim), filtration(_filt), location(_v), oldest_coface(-1) {vertices.push_back(_v);}
    delta_complex_cell_t(int _dim, std::vector<delta_complex_cell_t*> _cells,
        value_t _filt, value_t _loc) : dimension(_dim), boundary(_cells), filtration(_filt), location(_loc), oldest_coface(-1) {}

    void set_children(){
        for ( auto b : boundary){
            children.insert(b);
            children.insert(b->children.begin(),b->children.end());
            b->coboundary.push_back(location);
        }
    }

    void set_vertices(){
        for ( auto c : children ){
            if ( c->dimension == 0 ){
                vertices.push_back(c->vertices.front());
            }
        }
    }
    size_t coboundary_size(){
        return coboundary.size();
    }
    void set_filtration(value_t filter_val){
        filtration = filter_val;
    }
	void compute_oldest_coface(delta_complex_t* complex);
}; //END delta_complex_cell_t

class delta_complex_t {
public:
    //stores all cells as a vector of vectors, each vector being all cells of that dimension
    std::vector<std::vector<delta_complex_cell_t>> cells;

    //initialised with a string s with the address of a file containing all simplices
    //the format of this list is:
    //dim 0, followed by a line with an int for each vertex, the value of which is the filtration
    //after any subsequent dim i each line contained a list of facets of the simplex, followed
    //by a filtration value, the facets are indexed by their position in the list for dim i-1
    delta_complex_t(){};
	delta_complex_t(std::string s){
        std::ifstream infile;
        infile.open(s);

        index_t current_dimension = 0;
        int val = 0;
        std::string line;

        //read through all lines of input file, if a dim i is encountered then set
        //current dimension to i
        while (std::getline(infile,line)){
            if (line.length() == 0) continue;
            if (line[0] == 'd' && line[1] == 'i' && line[2] == 'm') {
                current_dimension = (int)line[4] - '0';
                cells.push_back( std::vector<delta_complex_cell_t>() );
            }
            //for dim 0 create a vertex cell for each entry on that line with filtration
            //value given by the entry
            else if (current_dimension == 0) {
                std::vector<value_t> vertex_filtration = split<value_t>(line, ' ', string_to_float);
                for( auto v : vertex_filtration ){
                    cells.back().push_back(delta_complex_cell_t(0,cells.back().size(),v));
                }
            //for all larger dimensions create a cell by giving the boundary faces and filtration value
        	} else {
        		std::vector<int> faces = split<int>(line, ' ', string_to_float);
                val = faces.back();
                faces.pop_back();
                std::vector<delta_complex_cell_t*> new_face;
                for (auto p : faces){
                    new_face.push_back(&cells[current_dimension-1][p]);
                }
                value_t loc = cells.back().size();
                cells.back().push_back(delta_complex_cell_t(current_dimension,new_face,val,loc));
                cells.back().back().set_children();
                if(current_dimension == 0){
                    cells.back().back().set_vertices();
                }
        	}
        }
	}
	delta_complex_t(std::vector<std::vector<std::vector<index_t>>>& faces){
        int val = 0;
        cells.resize(faces.size());
        //for dim 0 create a vertex cell for each entry on that line with filtration
	    //value given by the entry
        for( auto v : faces[0] ){
            cells[0].push_back(delta_complex_cell_t(0,cells[0].size(),0));
        }
        //read through all lines of input file, if a dim i is encountered then set
        //current dimension to i
        for (index_t current_dimension = 1; current_dimension < faces.size(); current_dimension++){
            for (auto f : faces[current_dimension]){
                val = 0;
                //f.pop_back();
                std::vector<delta_complex_cell_t*> new_face;
                for (auto p : f){
                    new_face.push_back(&cells[current_dimension-1][p]);
                }
                value_t loc = cells[current_dimension].size();
                cells[current_dimension].push_back(delta_complex_cell_t(current_dimension,new_face,val,loc));
                cells[current_dimension].back().set_children();
                if(current_dimension == 1){
                    cells[current_dimension].back().set_vertices();
                }
        	}
        }
	}
    index_t number_of_cells(index_t dimension) const {
        if ( dimension >= cells.size() ) { return 0; }
        return cells[dimension].size();
    }
    bool is_top_dimension(index_t dimension){
        return cells.size() == dimension-1;
    }
    index_t top_dimension(){
        return cells.size()-1;
    }
    value_t filtration(index_t dimension,index_t index){
        return cells[dimension][index].filtration;
    }
	std::vector<value_t> vertex_filtration(){
		std::vector<value_t> out;
		for( auto c : cells[0] ){
			out.push_back(c.filtration);
		}
        return out;
    }
    delta_complex_cell_t* get(index_t dimension,index_t index){
        return &cells[dimension][index];
    }
	void compute_oldest_cofaces(){
		for (auto p : cells){
			for(auto q : p){
				q.compute_oldest_coface(this);
			}
		}
	}
    void print(){
        for(int d = 0; d < cells.size(); d++){
            for(auto f : cells[d]){
                if(d==0){std::cout << f.location << " : ";}
                else{
	                for(auto i : f.boundary){
	                    std::cout << i->location << " ";
	                }
                }
                std::cout << " : ";
				for(auto i : f.coboundary){
					std::cout << i << " ";
				}
                std::cout << " : ";
				for(auto i : f.children){
					std::cout << "("<<i->dimension<<","<<i->location<<")" << " ";
				}
                std::cout << " : " << f.location << " : " << f.dimension << " : " << f.filtration;
                std::cout << std::endl;
            }
        }
    }
//END delta_complex_t
};

void delta_complex_cell_t::compute_oldest_coface(delta_complex_t* complex){
	 value_t oldest = -1;
	for( auto c : coboundary ){
		value_t f = complex->get(dimension+1,c)->filtration;
		if( f > oldest ){
			oldest = f;
			oldest_coface = c;
		}
	}
}

class simplex_coboundary_enumerator {
private:
	index_t  idx_above, dim;
	delta_complex_cell_t* simplex;
	delta_complex_t* complex;

public:
	simplex_coboundary_enumerator(const filtration_entry_t _simplex, index_t _dim,
	                              delta_complex_t* _complex)
	    : idx_above(0), simplex(_complex->get(_dim,_simplex.second)), dim(_dim), complex(_complex) {}

	bool has_next() {
		return idx_above < simplex->coboundary_size();
	}

	filtration_entry_t next() {
        idx_above++;
        return filtration_entry_t(complex->get(dim+1,simplex->coboundary[idx_above-1])->filtration, simplex->coboundary[idx_above-1]);
	}
};

#ifdef SORT_COLUMNS_BY_PIVOT
struct greater_filtration_or_better_pivot_or_smaller_index {
	greater_filtration_or_better_pivot_or_smaller_index(delta_complex_t* _complex, index_t _dimension) : complex(_complex), dimension(_dimension) {}
	bool operator()(filtration_index_t a, filtration_index_t b) const {
		// First order by the filtration value
		if (get_filtration(a) > get_filtration(b)) return true;
		if (get_filtration(a) < get_filtration(b)) return false;

		auto ta = get_coboundary_size_and_gap_and_pivot(get_index(a));
		auto tb = get_coboundary_size_and_gap_and_pivot(get_index(b));

		// Then the number of non-trivial coboundary entries
		if (std::get<0>(ta) < std::get<0>(tb)) return true;
		if (std::get<0>(ta) > std::get<0>(tb)) return false;

		// Then order by the better pivoting
		if (std::get<2>(ta) < std::get<2>(tb)) return true;
		if (std::get<2>(ta) > std::get<2>(tb)) return false;

		if (std::get<1>(ta) > std::get<1>(tb)) return true;
		if (std::get<1>(ta) < std::get<1>(tb)) return false;

		// Finally, order by their indices
		return get_index(a) < get_index(b);
	}

private:
	delta_complex_t* complex;
	index_t dimension;

	// A column is considered to be a better pivot if the jump from pivot to the next
	// non-trivial element is as big as possible. This prevents accidentally inserting
	// non-trivial elements just below the pivot, which sometimes creates very long
	// reduction chains.
	// The second sort criterium is for it to be small because the small pivots will be
	// used the most.
	std::tuple<size_t, size_t, index_t> get_coboundary_size_and_gap_and_pivot(index_t a) const {
		// Look at the first two gaps of the pivot and the next element
		index_t pivot = 0;
		size_t gap_after_pivot = 0;
		simplex_coboundary_enumerator iterator(a,dimension,complex);
		size_t coboundary_size = 0;
		while (iterator.has_next()) {
			coboundary_size++;
			index_t next_index = get_index(iterator.next().second);
			if (next_index > pivot) {
				gap_after_pivot = next_index - pivot;
				pivot = next_index;
			}
		}

		return std::make_tuple(coboundary_size, gap_after_pivot, pivot);
	}
};
#endif

// This class is just an ordinary priority queue, but once the
// queue gets too long (because a lot of faces are inserted multiple
// times) it starts collecting the coefficients and only inserting each
// new face once
template <class Container, class Comparator>
class priority_queue_t : public std::priority_queue<filtration_entry_t, Container, Comparator> {
	std::unordered_map<index_t, bool> coefficients;
	static const filtration_entry_t dummy;
	bool use_dense_version = false;
	size_t dense_threshold;

public:
	priority_queue_t(size_t _dense_threshold)
	    : dense_threshold(_dense_threshold) {}

	void push(const filtration_entry_t& value) {
		if (use_dense_version) {
			// If we already have this value: update the count and don't push it again
			auto p = coefficients.find(get_index(value));
			if (p != coefficients.end()) {
				p->second = !p->second;
				return;
			}
		}

		std::priority_queue<filtration_entry_t, Container, Comparator>::push(value);

		if (use_dense_version) coefficients.insert(std::make_pair(get_index(value), get_coefficient(value)));

		if (!use_dense_version &&
		    std::priority_queue<filtration_entry_t, Container, Comparator>::size() >= dense_threshold)
			use_dense_version = true;
	}

	void pop() {
		// Don't use this, only allow get_pivot
		throw std::exception();
	}

	filtration_entry_t pop_pivot() {
		remove_trivial_coefficient_entries();
		if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
			return dummy;
		else {
			auto pivot = get_top();
			safe_pop();
			while (!std::priority_queue<filtration_entry_t, Container, Comparator>::empty() &&
			       get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()) ==
			           get_index(pivot)) {
				safe_pop();
				remove_trivial_coefficient_entries();

				if (std::priority_queue<filtration_entry_t, Container, Comparator>::empty())
					return dummy;
				else {
					pivot = get_top();
					safe_pop();
				}
			}
			return pivot;
		}
	}

	filtration_entry_t get_pivot() {
		filtration_entry_t result = pop_pivot();
		if (get_index(result) != -1) { push(result); }
		return result;
	}

private:
	inline filtration_entry_t get_top() {
		auto pivot = std::priority_queue<filtration_entry_t, Container, Comparator>::top();
		return pivot;
	}

	inline void safe_pop() {
		if (use_dense_version) {
			auto e =
			    coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			if (e != coefficients.end()) coefficients.erase(e);
		}
		std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
	}

	inline void remove_trivial_coefficient_entries() {
		if (use_dense_version) {
			auto p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			while (p != coefficients.end() && p->second == false) {
				coefficients.erase(p);
				std::priority_queue<filtration_entry_t, Container, Comparator>::pop();
				p = coefficients.find(get_index(std::priority_queue<filtration_entry_t, Container, Comparator>::top()));
			}
		}
	}
};
template <class Container, class Comparator>
const filtration_entry_t priority_queue_t<Container, Comparator>::dummy(filtration_entry_t(0.0, -1));



//END delta_complex
//-------------------------------------------------------------------------//
//BEGIN deltser


class deltser {
	delta_complex_t* complex;
	std::ofstream outfile;
	index_t n, dim_max;
	bool python;
	mutable std::vector<filtration_entry_t> coface_entries;
	size_t max_entries;
	std::vector<size_t> skipped_entries;

public:
	std::vector<std::vector<std::pair<value_t,value_t>>> finite_pairs;
	std::vector<std::vector<value_t>> infinite_pairs;
	deltser(	delta_complex_t* _complex, std::string _outname, size_t _max_entries, bool _python)
	    : complex(_complex), n(complex->number_of_cells(0)),
	      dim_max(complex->top_dimension()),
		  max_entries(_max_entries),
		  python(_python) {
				if(!python) outfile.open(_outname);
                skipped_entries.assign(dim_max+1,0);
				infinite_pairs.resize(dim_max+1);
				finite_pairs.resize(dim_max+1);
			  }

	value_t compute_filtration(const index_t index, index_t dim) const {
        return complex->get(dim,index)->filtration;
	}

	void assemble_columns_to_reduce(std::vector<filtration_index_t>& simplices,
	                                std::vector<filtration_index_t>& columns_to_reduce,
	                                pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices);

	void compute_dim_0_pairs(std::vector<filtration_index_t>& edges,
	                         std::vector<filtration_index_t>& columns_to_reduce) {
        //Get all edges and sort them
		filtered_union_find dset(complex->vertex_filtration());
		edges = get_edges();
		std::sort(edges.rbegin(), edges.rend(),
		          greater_filtration_or_smaller_index());

		for (auto e : edges) {
			value_t birth = dset.link(complex->get(1,get_index(e))->vertices[0], complex->get(1,get_index(e))->vertices[1]);
			if (birth == -1) columns_to_reduce.push_back(e);
		}
		std::reverse(columns_to_reduce.begin(), columns_to_reduce.end());

		for (index_t i = 0; i < n; ++i)
			if (dset.find(i) == i){
				infinite_pairs[0].push_back(0);
			}
	}

	template <typename Column, typename Iterator>
	filtration_entry_t add_coboundary_and_get_pivot(Iterator column_begin, Iterator column_end,
	                                              Column& working_coboundary, const index_t& dim,
                                                  std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>&);

	void sort_columns(std::vector<filtration_index_t>& columns_to_reduce, index_t dimension) {
#ifdef SORT_COLUMNS_BY_PIVOT
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
			greater_filtration_or_better_pivot_or_smaller_index(complex,dimension));
#else
		std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	        greater_filtration_or_smaller_index());
#endif
}

	void compute_pairs(std::vector<filtration_index_t>& columns_to_reduce,
	                   pivot_column_index_t& pivot_column_index, index_t dim) {


		std::cout << "Computing Dimension " << dim << std::endl;
        compressed_sparse_matrix<filtration_entry_t> reduction_coefficients;

		std::vector<filtration_entry_t> coface_entries;

		for (index_t index_column_to_reduce = 0; index_column_to_reduce < columns_to_reduce.size();
		     ++index_column_to_reduce) {
			auto column_to_reduce = columns_to_reduce[index_column_to_reduce];
            std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>> reduction_column;

			priority_queue_t<std::vector<filtration_entry_t>,
			                    greater_filtration_or_smaller_index>
			    working_coboundary(columns_to_reduce.size());

			value_t filtration = get_filtration(column_to_reduce);

			index_t index_column_to_add = index_column_to_reduce;

			filtration_entry_t pivot;

            reduction_coefficients.append_column();
            reduction_coefficients.push_back(filtration_entry_t(column_to_reduce, 1));

			while (true) {
                auto reduction_column_begin = reduction_coefficients.cbegin(index_column_to_add);
                auto reduction_column_end = reduction_coefficients.cend(index_column_to_add);

				pivot = add_coboundary_and_get_pivot(
				    reduction_column_begin, reduction_column_end,
				    working_coboundary, dim, reduction_column);
				if (get_index(pivot) > -1) {
                    auto pivot_column_idx = pivot_column_index[get_index(pivot)];

                    if (pivot_column_idx != INVALID_INDEX) {
                        index_column_to_add = pivot_column_idx;
                        continue;
                    } else {
						value_t death = get_filtration(pivot);

                        pivot_column_index[get_index(pivot)] =  index_column_to_reduce;
                            reduction_coefficients.pop_back();
            				while (true) {
            					filtration_entry_t e = pop_pivot(reduction_column);
            					if (get_index(e) == -1) break;
            					reduction_coefficients.push_back(e);
            				}
						break;
					}
				} else if(get_index(pivot) == -1) {
					infinite_pairs[dim].push_back(filtration);
					break;
				}else {
					skipped_entries[dim]++;
					break;
				}
			}
		}
	}

	std::vector<filtration_index_t> get_edges();

	std::vector<value_t> num_infinite_pairs(){
		std::vector<value_t> out;
		for(auto i : infinite_pairs) out.push_back(i.size());
		return out;
	}

	void print_summary(){
		if(!python){
			std::vector<value_t> inf_pairs = num_infinite_pairs();
			outfile << std::endl;
			outfile << "# Betti Numbers:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim H_" << i << " = " << inf_pairs[i];
				if( skipped_entries[i] > 0){
					outfile << " : (" << skipped_entries[i] << " entries skipped)";
				}
				outfile << std::endl;
			}
			outfile << std::endl;
			outfile << "# Cell Counts:" << std::endl;
			for(index_t i = 0; i <= dim_max; i++){
				outfile << "#        dim C_" << i << " = " << complex->number_of_cells(i) << std::endl;
			}
		}
	}

	void compute_barcodes() {
		std::vector<filtration_index_t> simplices, columns_to_reduce;
		compute_dim_0_pairs(simplices, columns_to_reduce);
		for (index_t dim = 1; dim <= dim_max; ++dim) {
			pivot_column_index_t pivot_column_index(complex->number_of_cells(dim + 1), INVALID_INDEX);
			sort_columns(columns_to_reduce,dim);
			compute_pairs(columns_to_reduce, pivot_column_index, dim);
			if (dim < dim_max) {
				assemble_columns_to_reduce(simplices, columns_to_reduce, pivot_column_index,
				                           dim + 1, complex->number_of_cells(dim+1));
			}
		}
		print_summary();
	}
};

template <typename Column, typename Iterator>
filtration_entry_t deltser::add_coboundary_and_get_pivot(
    Iterator column_begin, Iterator column_end,
    Column& working_coboundary, const index_t& dim,
    std::priority_queue<filtration_entry_t, std::deque<filtration_entry_t>, smaller_index<filtration_entry_t>>& reduction_column) {
	index_t iterations = 0;
	for (auto it = column_begin; it != column_end; ++it) {
		filtration_entry_t simplex = *it;

        reduction_column.push(simplex);

		coface_entries.clear();
		simplex_coboundary_enumerator cofaces(simplex, dim, complex);
		while (cofaces.has_next()) {
			filtration_entry_t coface = cofaces.next();

            iterations++;
            working_coboundary.push(coface);
		}
		if (iterations > max_entries) {
			return filtration_entry_t(0,-2);
		}
	}

	return working_coboundary.get_pivot();
}

//returns a vector of all the edges where each is representated as a pair:
//(filtration,index), where the edge is the simplex at complex.get(1,index)
std::vector<filtration_index_t> deltser::get_edges() {
	std::vector<filtration_index_t> edges;
    int n = complex->number_of_cells(1);
	for ( index_t index = 0; index < n; index++) {
		edges.push_back(std::make_pair(complex->get(1,index)->filtration, index));
	}
	return edges;
}

void deltser::assemble_columns_to_reduce(
    std::vector<filtration_index_t>& simplices, std::vector<filtration_index_t>& columns_to_reduce,
    pivot_column_index_t& pivot_column_index, index_t dim, index_t num_simplices ) {

	columns_to_reduce.clear();

	for (index_t index = 0; index < num_simplices; ++index) {
		if (pivot_column_index[index] == INVALID_INDEX) {
			value_t filtration = compute_filtration(index, dim);
			columns_to_reduce.push_back(std::make_pair(filtration, index));
		}
	}

	//std::sort(columns_to_reduce.begin(), columns_to_reduce.end(),
	//          greater_filtration_or_smaller_index());
}


//END deltser


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
    //initialise approximate functionality
	size_t max_entries = std::numeric_limits<size_t>::max();
	if(argc > 5) max_entries = atoi(argv[5]);

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


    delta_complex_t complex(list_by_face);
	complex.compute_oldest_cofaces();

    //create deltser object and compute persistent homology
	deltser(&complex,out_address,max_entries,false).compute_barcodes();

    std::cout<<"Finished" <<std::endl;

}
