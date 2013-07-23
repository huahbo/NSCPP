include make.math

CCC = g++
EXECS = run

#INTERNAL DEFINITIONS
INCLUDES = -I /home/cipia/NSCPP/include $(INCLUDEMATH) 
CCFLAGS =  $(CFMATH) -lm

$(EXECS): prueba.o
	$(CCC) $(CCFLAGS) $(INCLUDES) -o $@ $< $ $(LIBMATH)
prueba.o: prueba.cpp
	$(CCC) $(CCFLAGS) $(INCLUDES) -c $<

