FC       = gfortran
FFLAGS   = -Wall -Wextra -O2 -std=f2008
LDFLAGS  = -lm

TARGET = heat2d
SRC    = heat_2d.f95

.PHONY: default build sim convergence all animation run clean clean_outputs 

default: build

build: $(TARGET)

$(TARGET): $(SRC)
	$(FC) $(FFLAGS) -o $@ $^ $(LDFLAGS)

sim: build
	./$(TARGET) sim

convergence: build
	./$(TARGET) convergence

all: build
	./$(TARGET) all

animation: build
	./$(TARGET) animation

run: build
	./$(TARGET)

clean:
	rm -f $(TARGET) *.o *.mod

clean_outputs:
	-rm -f data/*.dat
	-rm -f plot/*.png
	-rm -f animation/png/*.png
	-rm -f animation/dat/*.dat

