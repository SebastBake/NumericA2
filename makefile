#***************************************************************************
#
#   File        : makefile
#   Student Id  : 757931
#   Name        : Sebastian Baker
#
#***************************************************************************


CC=gcc
CFLAGS=-Wall # -Werror #-g -v -o0
OUT=bin/main
IN=input_files/
SRC=src/thomas_alg.o src/newton_raphson.o src/tasks.o src/main.o
RUN=$(OUT) $(IN)in_shock.csv $(IN)in_linalsys.csv $(IN)in_interp.csv 5 $(IN)in_heateqn.csv

compile: $(SRC)
	$(CC) $(SRC) $(CFLAGS) -o $(OUT);

clean:
	rm -f src/*.o;

test: clean compile
	$(RUN)

val: clean compile
	valgrind --leak-check=full --track-origins=yes $(RUN);