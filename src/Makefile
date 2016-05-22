CXX = g++
CXXFLAGS = -std=c++11 -O0 -g
CPPFLAGS = -Iinclude
all: hxx_gen-1 \
     hxx_viterbi-1 \
     hxx_baumwelch-1 \
     hxx_forward-1

hxx_baumwelch-1: hxx_baumwelch-1.cpp hxx_forward.o hxx_gen.o hxx_forwardbackward.o hxx_matrices.o hxx_baumwelch.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

hxx_viterbi-1: hxx_viterbi-1.cpp hxx_viterbi.o hxx_forward.o hxx_gen.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

hxx_forward-1: hxx_forward-1.cpp hxx_forward.o hxx_gen.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

hxx_gen-1: hxx_gen-1.cpp hxx_gen.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf *.o
