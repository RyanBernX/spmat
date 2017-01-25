CC=gcc
OBJ=main.o sparse_matrix.o
LIBS=-lm

%.o: %.c
	$(CC) -c -g $<

.PHONY: exec default clean

default: exec

exec: $(OBJ)
	$(CC) -o main.x $(OBJ) $(LIBS)

clean:
	rm -f $(OBJ) main.x
