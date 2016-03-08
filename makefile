CC = gcc
CFLAGS = -O3
LFLAGS = -lm -O3

Objects = lattice2d.o snapshot.o
RandomGen = mt19937-64.o
RandomObjs = iterate.o initialization.o

lattice2d : $(Objects) $(RandomGen) $(RandomObjs)
	$(CC) -o lattice2d $(LFLAGS) $(Objects) $(RandomGen) $(RandomObjs)

$(Objects) : lattice2d.h
$(RandomGen) : mt64.h
$(RandomObjs) : lattice2d.h mt64.h

.PHONY : clean
clean :
	-rm *.o
