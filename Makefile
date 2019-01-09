CXX?=g++
CXXFLAGS=-std=c++11 -Wall -Wextra
LDFLAGS=


SOURCES= commandline_flags.cc dag_align.cc freddie.cc
HEADERS= commandline_flags.h dag_align.h global.h

# OBJECTS   = $(SOURCES:.cc=.o)
.PHONY: clean

freddie: $(SOURCES) $(HEADERS)
	$(CXX) $(SOURCES) $(HEADERS) -o $@ $(LDFLAGS) $(CXXFLAGS) -O3

freddie_gdb: $(SOURCES) $(HEADERS)
	$(CXX) $(SOURCES) $(HEADERS) -o $@ $(LDFLAGS) $(CXXFLAGS) -g -O0

clean:
	rm -f freddie_gdb freddie

all: clean freddie
