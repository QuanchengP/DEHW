ifdef OS
.PHONY:
INCL = /usr/local/include
LIB = /usr/local/lib
all:
	g++ Test.cpp -o Test.exe -Og -fopenmp -DFILIB_EXTENDED -frounding-math -I$(INCL) -L$(LIB) -lprim
else
.PHONY:
INCL1 = /usr/local/include
INCL2 = /home/quanchengp/Library/petsc/include
INCL3 = /home/quanchengp/Library/petsc/arch-linux-c-debug/include
LIB1 = /usr/local/lib
LIB2 = /home/quanchengp/Library/petsc/arch-linux-c-debug/lib
all:
	g++ Test.cpp -o Test -Og -fopenmp -DFILIB_EXTENDED -frounding-math -I$(INCL1) -I$(INCL2) -I$(INCL3) -L$(LIB1) -L$(LIB2) -Wl,-rpath,$(LIB2) -lprim -lpetsc -lmpi
endif