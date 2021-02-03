#CC=gcc -Wall -g -pg -fopenmp -std=c99 -std=gnu99
#CC=gcc -fopenmp -std=c99 -std=gnu99 

CC=gcc -g -O2 -std=c99 -std=gnu99  
CC1 = gcc -c
CFLAGS=-Wall -fPIC -shared
LFLAGS=-lm

.PHONY: all

all: exoDiversity 
	@chmod +X preRunChecks.sh
	@sh preRunChecks.sh

exoDiversity: motifAndReadsFunctions.o modelfunctions.o bestModelFunctions.o traindata.o
	$(CC) $(CFLAGS) -o dataStructures.so messages.c motifAndReadsFunctions.c modelfunctions.c bestModelFunctions.c traindata.c  $(LFLAGS)
	echo "python $(CURDIR)/mainFunc.py "\"\$$@\" > exoDiversity; chmod +x exoDiversity
	cat cst1 > cstructures.py; echo "datalib = cdll.LoadLibrary(\"$(CURDIR)/dataStructures.so\")" >> cstructures.py; cat cst2 >> cstructures.py
clean:
	rm -f *.o *.pyc
