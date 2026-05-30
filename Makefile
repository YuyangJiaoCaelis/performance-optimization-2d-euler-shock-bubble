CXX ?= g++
TARGET ?= shock_bubble_solver
SOURCE ?= solver.cpp

CPPFLAGS ?=
CXXFLAGS ?= -O3 -march=native -std=c++17 -Wall -Wextra -Wpedantic
LDFLAGS ?=

all: $(TARGET)

$(TARGET): $(SOURCE)
	@mkdir -p $(dir $(TARGET))
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(SOURCE) -o $(TARGET) $(LDFLAGS)

debug: CXXFLAGS = -O0 -g -std=c++17 -Wall -Wextra -Wpedantic
debug: $(TARGET)

test: $(TARGET)
	@mkdir -p output/smoke
	./$(TARGET) --case planar-shock --nx 40 --ny 16 --final-time 1e-6 --output-prefix output/smoke/planar
	python3 scripts/validate_outputs.py output/smoke/planar_0001.csv --expect-y-invariant

clean:
	rm -f $(TARGET)
	rm -rf build output

.PHONY: all debug test clean
