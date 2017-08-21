#This is the makefile for debug

CXX = g++ 

CPPFLAGS = -I/home/zhang/lib/include  #-O2

cpp-files = LinearSolver.cpp MatVec.cpp DataStruct.cpp Driver.cpp \
            Bioprecip.cpp Build.cpp IO.cpp fish_subroutine.cpp

Tcpp-files =  LinearSolver.cpp MatVec.cpp DataStruct.cpp Driver.cpp \
              TestFast.cpp Build.cpp IO.cpp fish_subroutine.cpp

Diff-files = MatVec.cpp DataStruct.cpp DiffEP_1D.cpp  IO.cpp LinearSolver.cpp

o-files = $(cpp-files: .cpp=.o)
To-files = $(Tcpp-files: .cpp=.o)
Diff-ofiles = $(Diff-files: .cpp=.o)

Bioprecip:$(o-files)
	$(CXX) $(CPPFLAGS) -o $@ $^ /usr/lib/libfishpak_db_cpp.a /usr/lib/libf2c.a -lm
TestFast:$(To-files)
	$(CXX) $(CPPFLAGS) -o $@ $^ /usr/lib/libfishpak_db_cpp.a /usr/lib/libf2c.a -lm
DiffEP_1D:$(Diff-ofiles)
	$(CXX) $(CPPFLAGS) -o $@ $^ -lm
Bioprecip.o:Bioprecip.cpp Driver.h MatVec.h DataStruct.h LinearSolver.h
	$(CXX) $(CPPFLAGS) -c Bioprecip.cpp -o $@
Driver.o:Driver.cpp Driver.h 
	$(CXX) $(CPPFLAGS) -c Driver.cpp -o $@
LinearSolver.o:LinearSolver.h LinearSolver.cpp
	$(CXX) $(CPPFLAGS) -c LinearSolver.cpp -o $@
MatVec.o:MatVec.h MatVec.cpp
	$(CXX) $(CPPFLAGS) -c MatVec.cpp -o $@
DataStruct.o:DataStruct.h DataStruct.cpp MatVec.h
	$(CXX) $(CPPFLAGS) -c DataStruct.cpp -o $@
Build.o:Build.h Build.cpp
	$(CXX) $(CPPFLAGS) -c Build.cpp -o $@
IO.o:IO.h IO.cpp
	$(CXX) $(CPPFLAGS) -c IO.cpp -o $@
fish_subroutine.o:fish_subroutine.h fish_subroutine.cpp
	$(CXX) $(CPPFLAGS) -c fish_subroutine.cpp -o $@
TestFast.o: TestFast.cpp Build.h IO.h DataStruct.h
	$(CXX) $(CPPFLAGS) -c TestFast.cpp -o $@
DiffEP_1D.o:IO.h MatVec.h
	$(CXX) $(CPPFLAGS) -c DiffEP_1D.cpp -o $@
