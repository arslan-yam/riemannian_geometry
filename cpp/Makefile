# C++ implementation of the Certified Cut Index and its tools. `make` builds everything into bin/.
CXX ?= c++
CXXFLAGS ?= -O3 -std=c++17 -DNDEBUG -Wno-deprecated-declarations -pthread
KIR = -Ithird_party/kirsanov
VTPI = -Ithird_party/vtp
BIN = bin
OBJ = $(BIN)/obj

all: $(BIN)/solver_check $(BIN)/cci_bench $(BIN)/cci_exp

# competitors (their code, not part of `all`): the UP-Oracle authors' oracles behind a static driver
UPO ?= ../external/UpdatedStructureTerrainCode/src
SHIM ?= ../external/weighted_distance_oracle/shim
$(BIN)/upo_static: tools/upo_static.cpp
	@mkdir -p $(BIN); $(CXX) -O3 -std=c++14 -w -I$(UPO) -I$(SHIM) $< -o $@

$(OBJ)/solver_mmp.o: src/solver_mmp.cpp src/solver.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -std=c++14 $(KIR) -c $< -o $@
$(OBJ)/solver_vtp.o: src/solver_vtp.cpp src/solver.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -std=c++14 $(VTPI) -c $< -o $@
$(BIN)/solver_check: tools/solver_check.cpp $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o -o $@

$(OBJ)/certified.o: src/certified.cpp src/certified.h src/terrain.h src/graph.h src/solver.h
	@mkdir -p $(OBJ); $(CXX) $(CXXFLAGS) -c $< -o $@
$(BIN)/cci_bench: tools/cci_bench.cpp $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o -o $@

$(BIN)/cci_exp: tools/cci_exp.cpp $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o src/mesh_io.h
	$(CXX) $(CXXFLAGS) $< $(OBJ)/certified.o $(OBJ)/solver_mmp.o $(OBJ)/solver_vtp.o -o $@

clean:
	rm -rf $(BIN)
