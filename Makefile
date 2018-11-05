CXX      ?= g++
CXXFLAGS  = -O3 -std=c++14
LDFLAGS   =

all: clean freddie 

SOURCES   = $(wildcard *.cc)
OBJECTS   = $(SOURCES:.cc=.o)

freddie: $(OBJECTS)
	$(CXX) $(OBJECTS) -o $@ ${LDFLAGS}

clean:
	@rm -f $(OBJECTS) freddie
