include ../make.math

CCC = g++
EXECS = tester

#INTERNAL DEFINITIONS
INCLUDES = /home/cipia/NSCPP/include
CCFLAGS =  -lm -Wall -g -Wno-reorder

OBJS = \
prueba.o \
pressureCorrection.o \
mgmres.o\
Mesh.o\
MathCore.o\
MathInterpolation.o\
Cells.o\
Momento.o\
RungeKutta.o

$(EXECS): $(OBJS) /home/cipia/NSCPP/include/*.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -o $@ $(OBJS) 
prueba.o: prueba.cpp /home/cipia/NSCPP/include/*.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
MathCore.o: MathCore.cpp ../include/MathCore.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
pressureCorrection.o: pressureCorrection.cpp ../include/pressureCorrection.h 
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
mgmres.o: mgmres.cpp ../include/mgmres.hpp
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
Mesh.o: Mesh.cpp ../include/Mesh.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
RungeKutta.o: RungeKutta.cpp ../include/RungeKutta.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
MathInterpolation.o: MathInterpolation.cpp ../include/MathInterpolation.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
Cells.o: Cells.cpp ../include/Cells.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
Momento.o: Momento.cpp ../include/Momento.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
clean:
	rm *.o
