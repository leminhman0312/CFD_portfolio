SHELL := /bin/bash

FC      := gfortran
FFLAGS  := -Wall -Wextra -O3 -std=f2008 -J build -I build
LDFLAGS := -lm

TARGET  := heat2d
SRC     := src
BUILD   := build
BIN     := bin
EXE     := $(BIN)/$(TARGET)

OUTDIRS := data plot gnuplot_scripts animation animation/dat animation/png animation/movie

OBJS := \
  $(BUILD)/kinds.o \
  $(BUILD)/io_dirs.o \
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

.PHONY: default dirs build sim convergence animation all run clean clean_outputs

default: build

dirs:
	@mkdir -p "$(BUILD)" "$(BIN)"
	@mkdir -p $(OUTDIRS)

build: dirs $(EXE)

$(EXE): $(OBJS)
	$(FC) -o "$@" $^ $(LDFLAGS)

$(BUILD)/%.o: $(SRC)/%.f95 | dirs
	$(FC) $(FFLAGS) -c "$<" -o "$@"

sim: build
	"$(EXE)" sim

convergence: build
	"$(EXE)" convergence

animation: build
	"$(EXE)" animation

all: build
	"$(EXE)" all

run: build
	"$(EXE)"

clean_outputs:
	@rm -f -- data/*.dat plot/*.png
	@rm -f -- animation/movie/*.mp4 animation/dat/*.dat animation/png/*.png

clean:
	@set -e; \
	if [ -z "$(BUILD)" ] || [ "$(BUILD)" = "/" ] || [ "$(BUILD)" = "." ]; then \
	  echo "Refusing to clean. BUILD is unsafe: '$(BUILD)'"; exit 2; \
	fi; \
	if [ -z "$(BIN)" ] || [ "$(BIN)" = "/" ] || [ "$(BIN)" = "." ]; then \
	  echo "Refusing to clean. BIN is unsafe: '$(BIN)'"; exit 2; \
	fi; \
	rm -rf -- "$(BUILD)"; \
	rm -f -- "$(EXE)"; \
	$(MAKE) clean_outputs
