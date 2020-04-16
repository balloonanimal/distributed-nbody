CC = mpicc
CFLAGS_WARN = -Wall -Wvla -Werror -Wno-error=unused-variable
CFLAGS = -std=c11 -g -O2 $(CFLAGS_WARN)

EXPERIMENTS = test

LIB_FLAGS = -I./lib/ -L. -lnbody -Wl,-rpath,./

all: library $(EXPERIMENTS)

library:
	@echo "Building Library"
	make -C ./lib/
	ln -s ./lib/build/libnbody.so .
	@echo "Done building Library"

test:
	$(CC) $(CFLAGS) $(LIB_FLAGS) test.c -o test

clean:
	unlink nbody.so
	make -C ./lib/ clean
	rm -f $(EXPERIMENTS)
