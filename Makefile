CXX      ?=g++
CXXFLAGS  =-O3 -std=c++14
LDFLAGS   =

all: clean freddie

SOURCES   = $(wildcard *.cc) $(wildcard *.h)
# OBJECTS   = $(SOURCES:.cc=.o)

freddie:
	$(CXX) $(SOURCES) -o $@ ${LDFLAGS}

clean:
	@rm -f $(OBJECTS) freddie
