.PHONY: run plot plot_streamline plot_iteration \
        plot_history plot_history_compare \
        plot_grid_refinement clean

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall

TARGET = laplace
SRC = laplace.cpp

run: $(TARGET)
	./$(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) $(SRC) -o $(TARGET)



plot: plot_streamline plot_iteration plot_history plot_history_compare plot_grid_refinement

plot_streamline:
	python plot_streamline.py point_gs
	python plot_streamline.py point_sor
	python plot_streamline.py line_gs
	python plot_streamline.py line_sor

plot_iteration:
	python plot_iteration.py

plot_history:
	python plot_history.py point_gs
	python plot_history.py point_sor
	python plot_history.py line_gs
	python plot_history.py line_sor

plot_history_compare:
	python plot_compare.py

plot_grid_refinement:
	python plot_grid_refinement.py

clean:
	rm -f $(TARGET)
	rm -f data/*.dat
	rm -f plot/*.png
