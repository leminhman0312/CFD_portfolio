CXX      = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17 -Iinclude
LDFLAGS  = -lm

TARGET = heat_2d

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

.PHONY: default all sim convergence clean

default:
	@echo "make sim | make convergence | make all | make clean"

$(TARGET): $(OBJ)
	$(CXX) -o $@ $^ $(LDFLAGS)

build/%.o: src/%.cpp include/heat_2d.hpp | build
	$(CXX) $(CXXFLAGS) -c $< -o $@

build:
	mkdir -p build

sim: $(TARGET)
	./$(TARGET) sim

convergence: $(TARGET)
	./$(TARGET) convergence

all: $(TARGET)
	./$(TARGET) all

clean:
	rm -f $(TARGET)
	rm -rf build
	-rm -f data/*.dat
	-rm -f plot/*.png
