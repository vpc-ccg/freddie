CXX=g++
CXXFLAGS=-std=c++14 -I src/ -I extern/ -Wall -Wextra
LDFLAGS=


SOURCES:=$(wildcard src/*.cc)
OBJECTS:=$(SOURCES:.cc=.o)
EXECUTABLE=freddie
DBG_OBJECTS:=$(SOURCES:.cc=.dbg)
DEBUGGABLE=$(EXECUTABLE)_dbg
.PHONY: clean all

$(EXECUTABLE): $(EXECUTABLE).cc $(OBJECTS) Makefile
	$(CXX) $(OBJECTS) -O3 $< $(CXXFLAGS) -o $@

$(DEBUGGABLE): $(EXECUTABLE).cc $(DBG_OBJECTS) Makefile
	$(CXX) $(CXXFLAGS) $(DBG_OBJECTS) -g -O1 $<  -o $@

%.o: %.cc Makefile
	$(CXX) $(CXXFLAGS) -O3 -c $< -o $@

%.dbg: %.cc Makefile
	$(CXX) $(CXXFLAGS) -g -O1 -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(DEBUGGABLE) $(OBJECTS) $(SOURCES:.cc=.h.gch) $(DBG_OBJECTS)

all: clean $(EXECUTABLE)
