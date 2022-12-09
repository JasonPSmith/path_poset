build: path_poset path_poset_sage path_poset_homology

path_poset: path_poset.cpp
	@echo "Compiling \"path_poset\"." && g++ -std=c++11 -pthread -O3 path_poset.cpp -o path_poset

path_poset_sage: path_poset_sage.cpp
	@echo "Compiling \"path_poset_sage\"." && g++ -std=c++11 -pthread -O3 path_poset_sage.cpp -o path_poset_sage

path_poset_deltser: path_poset_deltser.cpp
	@echo "Compiling \"path_poset_deltser\"." && g++ -std=c++11 -pthread -O3 path_poset_deltser.cpp -o path_poset_deltser

path_poset_homology: path_poset_homology.cpp
	@echo "Compiling \"path_poset_homology\"." && g++ -std=c++14 -pthread -O3 path_poset_homology.cpp -o path_poset_homology
