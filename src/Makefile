include ../make.math

CCC = g++
EXECS = tester

#INTERNAL DEFINITIONS
INCLUDES = /home/cipia/NSCPP/include $(INCLUDEMATH) 
CCFLAGS =  $(CFMATH) -lm

OBJS = \
prueba.o \
pressureCorrection.o \
Mesh.o

$(EXECS): $(OBJS) /home/cipia/NSCPP/include/*.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -o $@ $(OBJS) $(LIBMATH)
prueba.o: prueba.cpp /home/cipia/NSCPP/include/*.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
pressureCorrection.o: pressureCorrection.cpp ../include/pressureCorrection.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
Mesh.o: Mesh.cpp ../include/Mesh.h
	$(CCC) $(CCFLAGS) -I$(INCLUDES) -c $*.cpp
