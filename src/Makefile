CC = gcc
CFLAGS = -g -Wall -Wvla -Werror -Wno-error=unused-variable

DEPS = types.h simulation.h
OBJ = simulation.o

all: $(OBJ)

%.o : %.c $(DEPS)
	$(CC) $(CFLAGS) -o $@ -c $<
