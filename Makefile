CXX = g++
CXXFLAGS = -g -Wall -Wextra -std=c++17

objects := tde_lcurve_run.o tde_lcurve.o tde_lcurve_assistant.o Trapezoidal_rule/trapezoidal.o

obj_with_headers := tde_lcurve.o tde_lcurve_assistant.o Trapezoidal_rule/trapezoidal.o

all: $(objects) TiDE

TiDE: $(objects)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(obj_with_headers): %.o: %.cpp %.hpp
	$(CXX) $(CXXFLAGS) -c -o $@ $<

tde_lcurve_run.o: tde_lcurve_run.cpp
	$(CXX) $(CXXFLAGS) -c -o $@ $^

.PHONY: clean

clean:
	-rm TiDE $(objects)
