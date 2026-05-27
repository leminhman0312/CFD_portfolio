.PHONY: run plot clean

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

TARGET = laplace
SRC = laplace.cpp

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)

plot:
	python3 plot_streamline.py

clean:
	rm -f $(TARGET)
	rm -f data/psi.dat
	rm -f plot/streamline.png
