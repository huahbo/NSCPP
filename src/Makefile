include ../make.math

CCC = g++
EXECS = run

#INTERNAL DEFINITIONS
INCLUDES = -I /home/cipia/NSCPP/include $(INCLUDEMATH) 
CCFLAGS =  $(CFMATH) -lm

#$(EXECS): prueba.o
#	$(CCC) $(CCFLAGS) $(INCLUDES) -o $@ $< $ $(LIBMATH)
pressureCorrection.o: pressureCorrection.cpp
	$(CCC) $(CCFLAGS) $(INCLUDES) -c $<

