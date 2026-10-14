# C++ implementation of the Certified Cut Index and its tools. `make` builds everything into bin/.
CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -DNDEBUG -Wno-deprecated-declarations -pthread
# the third-party headers rely on includes that libc++ (clang) pulls in transitively and libstdc++ (GCC) does not
KIR = -Ithird_party/kirsanov -include cstring -include memory
VTPI = -Ithird_party/vtp -include cstring -include memory
BIN = bin
OBJ = $(BIN)/obj

all: $(BIN)/solver_check $(BIN)/cci_bench $(BIN)/cci_exp $(BIN)/proximity

# competitors (their code, not part of `all`): the UP-Oracle authors' oracles behind a static driver
UPO ?= ../external/UpdatedStructureTerrainCode/src
SHIM ?= ../external/weighted_distance_oracle/shim
$(BIN)/upo_static: tools/upo_static.cpp
	@mkdir -p $(BIN); $(CXX) -O3 -std=c++14 -w -I$(UPO) -I$(SHIM) $< -o $@

$(OBJ)/solver_mmp.o: src/solver_mmp.cpp src/solver.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -std=c++14 $(KIR) -c $< -o $@
$(OBJ)/solver_vtp.o: src/solver_vtp.cpp src/solver.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -std=c++14 $(VTPI) -c $< -o $@
$(BIN)/solver_check: tools/solver_check.cpp $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o -o $@

$(OBJ)/solver_steiner.o: src/solver_steiner.cpp src/solver.h src/steiner.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -c $< -o $@
$(OBJ)/certified.o: src/certified.cpp src/certified.h src/terrain.h src/graph.h src/solver.h src/steiner.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)/cci_bench: tools/cci_bench.cpp $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o -o $@

$(BIN)/cci_exp: tools/cci_exp.cpp $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o -o $@

$(BIN)/proximity: tools/proximity.cpp $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o $(OBJ)/solver_steiner.o -o $@

# DGG baseline (authors' code, ported in ../external/dgg_port; not part of `all`)
DGG ?= ../external/dgg_port
DGGSHIM ?= ../external/dgg_shim
dgg: $(BIN)/dgg_precompute $(BIN)/dgg_bench
$(BIN)/dgg_precompute: $(wildcard $(DGG)/DGG_precompute/*.cpp) $(wildcard $(DGG)/ICH/*.cpp)
	@mkdir -p $(BIN); $(CXX) -O2 -std=c++14 -w -include cfloat -include climits -include cmath -include cstring \
	  -I$(DGG) -I$(DGG)/ICH -I$(DGG)/DGG_precompute -I$(DGGSHIM) $^ -o $@
$(BIN)/dgg_bench: tools/dgg_bench.cpp
	@mkdir -p $(BIN); $(CXX) -O3 -std=c++14 -w -include limits -include cmath -I$(DGG) $< -o $@

clean:
	rm -rf $(BIN)
