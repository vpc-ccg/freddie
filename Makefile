CXX      ?=g++
CXXFLAGS  =-O3 -std=c++11
CXXDBG  =  -g -std=c++11
LDFLAGS   =

all: clean freddie

SOURCES   = $(wildcard *.cc) $(wildcard *.h)
# OBJECTS   = $(SOURCES:.cc=.o)

freddie:
	$(CXX) $(SOURCES) -o $@ ${LDFLAGS} $(CXXFLAGS)

freddie_gdb:
	$(CXX) $(SOURCES) -o $@ ${LDFLAGS} $(CXXDBG)

clean:
	@rm -f $(OBJECTS) freddie
