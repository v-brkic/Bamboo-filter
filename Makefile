CXX      := g++
CXXFLAGS := -std=c++17 -O2 -Wall -Wextra -pedantic

all: bamboo_filter          # default cilj

bamboo_filter: main.cpp
	$(CXX) $(CXXFLAGS) $< -o $@

clean:
	$(RM) bamboo_filter
