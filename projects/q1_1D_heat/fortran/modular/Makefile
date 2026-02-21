FC      := gfortran
FFLAGS  := -O2 -std=f2008 -Wall -Wextra

INC_DIR := include
SRC_DIR := src
BLD_DIR := build

EXE := 1D_heat

OBJS := \
  $(BLD_DIR)/headers.o \
  $(BLD_DIR)/helpers.o \
  $(BLD_DIR)/methods.o \
  $(BLD_DIR)/driver.o \
  $(BLD_DIR)/main.o

.PHONY: all run clean veryclean dirs

all: dirs $(EXE)

dirs:
	mkdir -p $(BLD_DIR) data plot

$(EXE): $(OBJS)
	$(FC) $(FFLAGS) $(OBJS) -o $@

$(BLD_DIR)/headers.o: $(INC_DIR)/headers.f95 | dirs
	$(FC) $(FFLAGS) -J$(BLD_DIR) -I$(BLD_DIR) -c $< -o $@

$(BLD_DIR)/helpers.o: $(INC_DIR)/helpers.f95 $(BLD_DIR)/headers.o | dirs
	$(FC) $(FFLAGS) -J$(BLD_DIR) -I$(BLD_DIR) -c $< -o $@

$(BLD_DIR)/methods.o: $(INC_DIR)/methods.f95 $(BLD_DIR)/headers.o | dirs
	$(FC) $(FFLAGS) -J$(BLD_DIR) -I$(BLD_DIR) -c $< -o $@

$(BLD_DIR)/driver.o: $(SRC_DIR)/driver.f95 $(BLD_DIR)/headers.o $(BLD_DIR)/helpers.o $(BLD_DIR)/methods.o | dirs
	$(FC) $(FFLAGS) -J$(BLD_DIR) -I$(BLD_DIR) -c $< -o $@

$(BLD_DIR)/main.o: $(SRC_DIR)/main.f95 $(BLD_DIR)/headers.o $(BLD_DIR)/helpers.o $(BLD_DIR)/driver.o | dirs
	$(FC) $(FFLAGS) -J$(BLD_DIR) -I$(BLD_DIR) -c $< -o $@

run: all
	./$(EXE)

clean:
	rm -rf $(BLD_DIR) $(EXE)

veryclean: clean
	rm -rf data plot
