CXX      = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17 -Iinclude
LDFLAGS  = -lm

TARGET = heat_2d
PLOT   = plot_heat_2d.py

SRC = \
  src/main.cpp \
  src/drivers.cpp \
  src/grid_init.cpp \
  src/explicit.cpp \
  src/implicit_adi.cpp \
  src/thomas.cpp \
  src/norms.cpp \
  src/io.cpp \
  src/plot.cpp \
  src/dirs.cpp

OBJ = $(SRC:src/%.cpp=build/%.o)

.PHONY: default build sim convergence all run clean clean_outputs clean_all

default:
	@echo "make sim | make convergence | make all | make run | make clean | make clean_outputs | make clean_all"

build: $(TARGET)

$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

build/%.o: src/%.cpp include/heat_2d.hpp | build_dir
	$(CXX) $(CXXFLAGS) -c $< -o $@

build_dir:
	mkdir -p build

sim: build $(PLOT)
	./$(TARGET) sim

convergence: build $(PLOT)
	./$(TARGET) convergence

all: build $(PLOT)
	./$(TARGET) all

run: build $(PLOT)
	./$(TARGET)

clean:
	rm -f $(TARGET)
	rm -f heat2d
	rm -f a.out
	rm -f *~
	rm -f include/*~
	rm -f src/*~
	rm -rf build

clean_outputs:
	-rm -f data/*.dat
	-rm -f plot/*.png
	-rm -f plot/*.png~

clean_all: clean clean_outputs
	-rm -rf gnuplot_scripts