UNAME_S := $(shell uname -s)
CCC = g++
TARGETDIR=.
CCFLAGS = -O3 -Wall -std=c++11 -Wno-unused-but-set-variable
ifeq ($(UNAME_S),Linux)
    CCFLAGS += -fopenmp
endif
ifeq ($(UNAME_S),Darwin)
    CCFLAGS += -openmp
endif

debug: CCFLAGS += -g 
release: CCFLAGS += -march=native
 
OBJS = \
		$(TARGETDIR)/Main/main.o \
        $(TARGETDIR)/Problem/Problem.o 

$(TARGETDIR)/runTest: $(OBJS)
	$(CCC) $(CCFLAGS) -o $(TARGETDIR)/runTest $(OBJS)
    
$(TARGETDIR)/Problem/Problem.o: $(TARGETDIR)/Problem/Problem.cpp
	$(CCC) $(CCFLAGS) -c $< -o $@

$(TARGETDIR)/Main/main.o: $(TARGETDIR)/Main/main.cpp
	$(CCC) $(CCFLAGS) -c $< -o $@

test: runTest
	./runTest

clean:
	$(RM) \
	$(TARGETDIR)/Main/main.o \
	$(TARGETDIR)/Problem/Problem.o \
	$(TARGETDIR)/runTest
    
debug: clean
	$(MAKE) runTest CCC=g++

rebuild: 
	$(MAKE) debug
