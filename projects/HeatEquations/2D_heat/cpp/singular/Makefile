CXX      = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17
LDFLAGS  = -lm

TARGET = heat2d
SRC    = heat_2d.cpp

.PHONY: default build sim convergence all run clean clean_outputs

default: build

build: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

sim: build
	./$(TARGET) sim

convergence: build
	./$(TARGET) convergence

all: build
	./$(TARGET) all

run: build
	./$(TARGET)

clean:
	rm -f $(TARGET)

clean_outputs:
	-rm -f data/*.dat
	-rm -f plot/*.png
