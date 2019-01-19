CXX=g++
CXXFLAGS=-std=c++14 -I src/ -I extern/ -Wall -Wextra
LDFLAGS=


SOURCES:=$(wildcard src/*.cc)
OBJECTS:=$(SOURCES:.cc=.o)
EXECUTABLE=freddie
DEBUGGABLE=$(EXECUTABLE)_dbg
.PHONY: clean all

$(EXECUTABLE): $(EXECUTABLE).cc $(OBJECTS) Makefile
	$(CXX) $(OBJECTS) -O3 $< $(CXXFLAGS) -o $@

$(DEBUGGABLE): $(EXECUTABLE).cc $(OBJECTS) Makefile
	$(CXX) $(OBJECTS) -g -O1 $< $(CXXFLAGS) -o $@

%.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) -O3 -c $< -o $@
	echo $@

clean:
	rm -f $(EXECUTABLE) $(DEBUGGABLE) $(OBJECTS) $(SOURCES:.cc=.h.gch)

all: clean $(EXECUTABLE)
