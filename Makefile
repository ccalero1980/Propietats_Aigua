.SUFFIXES:
.SUFFIXES: .cpp .o
.c.o:
	gcc -c $*.cpp

membrane: main.o water.o interactions.o heads.o
	gcc -o membrane main.o water.o interactions.o heads.o -lm

main.o water.o interactions.o heads.o: membrane.h

clean:
	rm -f main.o water.o interactions.o heads.o
