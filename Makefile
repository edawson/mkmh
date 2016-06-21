CXX:=g++
CXXFLAGS:= -O3 -mtune=native -std=c++11
LD_LIB_FLAGS:=
LD_INC_FLAGS:=

libmkmh: mkmh.o Makefile
	ar -rs $@ $<

mkmh.o: mkmh.cpp mkmh.hpp
	$(CXX) $(CXXFLAGS) -c $< $(LD_LIB_FLAGS) $(LD_INC_FLAGS)

.PHONY: clean clobber

clean:
	$(RM) *.o

clobber: clean
	$(RM) libmkmh
