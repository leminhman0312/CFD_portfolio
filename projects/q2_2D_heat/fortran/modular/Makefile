FC      = gfortran
FFLAGS  = -Wall -Wextra -O2 -std=f2008 -J build -I build
LDFLAGS = -lm

TARGET = heat2d
SRC    = src
BUILD  = build
BIN    = bin

OBJS = \
  $(BUILD)/kinds.o \
  $(BUILD)/dir_utils.o \
  $(BUILD)/norms.o \
  $(BUILD)/field_init.o \
  $(BUILD)/io_field.o \
  $(BUILD)/plotting.o \
  $(BUILD)/tridiag.o \
  $(BUILD)/solver_explicit.o \
  $(BUILD)/solver_adi.o \
  $(BUILD)/driver.o \
  $(BUILD)/animation.o \
  $(BUILD)/main.o

.PHONY: default build sim convergence all animation run clean clean_outputs dirs

default: build

dirs:
	mkdir -p $(BUILD) $(BIN) data plot gnuplot_scripts animation

build: dirs $(BIN)/$(TARGET)

$(BIN)/$(TARGET): $(OBJS)
	$(FC) -o $@ $^ $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.f95 | dirs
	$(FC) $(FFLAGS) -c $< -o $@

# Module dependency ordering
$(BUILD)/dir_utils.o: $(BUILD)/kinds.o
$(BUILD)/norms.o: $(BUILD)/kinds.o
$(BUILD)/field_init.o: $(BUILD)/kinds.o
$(BUILD)/io_field.o: $(BUILD)/kinds.o
$(BUILD)/plotting.o: $(BUILD)/kinds.o
$(BUILD)/tridiag.o: $(BUILD)/kinds.o

$(BUILD)/solver_explicit.o: $(BUILD)/kinds.o
$(BUILD)/solver_adi.o: $(BUILD)/kinds.o $(BUILD)/tridiag.o $(BUILD)/dir_utils.o $(BUILD)/io_field.o $(BUILD)/plotting.o

$(BUILD)/driver.o: $(BUILD)/kinds.o $(BUILD)/norms.o $(BUILD)/io_field.o $(BUILD)/plotting.o $(BUILD)/solver_explicit.o $(BUILD)/solver_adi.o
$(BUILD)/animation.o: $(BUILD)/kinds.o $(BUILD)/solver_adi.o
$(BUILD)/main.o: $(BUILD)/kinds.o $(BUILD)/dir_utils.o $(BUILD)/field_init.o $(BUILD)/driver.o $(BUILD)/animation.o

sim: build
	./$(BIN)/$(TARGET) sim

convergence: build
	./$(BIN)/$(TARGET) convergence

all: build
	./$(BIN)/$(TARGET) all

animation: build
	./$(BIN)/$(TARGET) animation

run: build
	./$(BIN)/$(TARGET)

clean:
	rm -rf $(BUILD) $(BIN)/$(TARGET)

clean_outputs:
	-rm -f data/*.dat
	-rm -f plot/*.png
	-rm -f animation/*.mp4
	-rm -f animation/dat/*.dat
	-rm -f animation/png/*.png
