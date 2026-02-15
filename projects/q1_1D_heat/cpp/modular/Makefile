CXX := g++
CXXFLAGS := -O2 -std=c++17 -Iinclude

OBJ := build/main.o build/helpers.o build/methods.o build/driver.o

# ------------------------------------------------------------------
# Build executable
# ------------------------------------------------------------------
1D_heat: $(OBJ)
	$(CXX) $(CXXFLAGS) -o 1D_heat $(OBJ)

# ------------------------------------------------------------------
# Object files
# ------------------------------------------------------------------
build/main.o: src/main.cpp include/headers.h include/helpers.h include/driver.h | build
	$(CXX) $(CXXFLAGS) -c src/main.cpp -o build/main.o

build/helpers.o: src/helpers.cpp include/helpers.h include/headers.h | build
	$(CXX) $(CXXFLAGS) -c src/helpers.cpp -o build/helpers.o

build/methods.o: src/methods.cpp include/methods.h include/headers.h | build
	$(CXX) $(CXXFLAGS) -c src/methods.cpp -o build/methods.o

build/driver.o: src/driver.cpp include/driver.h include/helpers.h include/methods.h include/headers.h | build
	$(CXX) $(CXXFLAGS) -c src/driver.cpp -o build/driver.o

# ------------------------------------------------------------------
# Ensure build directory exists
# ------------------------------------------------------------------
build:
	mkdir -p build

# ------------------------------------------------------------------
# Run
# ------------------------------------------------------------------
run: 1D_heat
	./1D_heat

# ------------------------------------------------------------------
# Clean
# ------------------------------------------------------------------
clean:
	rm -f data/*.txt plot/*.png build/*.o 1D_heat
