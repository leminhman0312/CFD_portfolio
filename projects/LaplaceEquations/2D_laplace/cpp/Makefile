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
	python plot_streamline.py point_gs
	python plot_streamline.py point_sor
	python plot_streamline.py line_gs
	python plot_streamline.py line_sor

clean:
	rm -f $(TARGET)
	rm -f data/*.dat
	rm -f plot/*.png
