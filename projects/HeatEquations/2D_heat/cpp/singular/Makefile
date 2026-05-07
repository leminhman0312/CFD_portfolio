CXX      = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++17
LDFLAGS  = -lm

TARGET = heat2d
SRC    = heat_2d.cpp
PLOT   = plot_heat_2d.py

.PHONY: default build sim convergence all run clean clean_outputs clean_all

default: build

build: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

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
	rm -f a.out
	rm -f *~

clean_outputs:
	-rm -f data/*.dat
	-rm -f plot/*.png
	-rm -f plot/*.png~

clean_all: clean clean_outputs
	-rm -rf gnuplot_scripts
