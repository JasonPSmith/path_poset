build: path_poset

path_poset: path_poset.cpp
	@echo "Compiling \"path_poset\"." && g++ -std=c++11 -pthread -O3 path_poset.cpp -o path_poset
