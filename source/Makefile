CC = mpicc
CFLAGS_WARN = -Wall -Wvla -Werror -Wno-error=unused-variable
CFLAGS = -std=c11 -g -O2 $(CFLAGS_WARN)

EXPERIMENTS = test benchmark benchmark_mpi

LIB_FLAGS = -I./lib/ -L. -lnbody -Wl,-rpath,./ -lm

all: $(EXPERIMENTS)

library:
	@echo "Building Library"
	@make -C ./lib/
	-@ln -sf ./lib/build/libnbody.so .
	-@install_name_tool -id "@rpath/libnbody.so" libnbody.so
	@echo "Done building Library"

test: library
	@echo "Building experiment test"
	$(CC) test.c -o test $(CFLAGS) $(LIB_FLAGS)

benchmark: library
	@echo "Building experiment benchmark"
	$(CC) benchmark.c -o benchmark $(CFLAGS) $(LIB_FLAGS)

benchmark_mpi: library
	@echo "Building experiment benchmark_mpi"
	$(CC) benchmark_mpi.c -o benchmark_mpi $(CFLAGS) $(LIB_FLAGS)

clean:
	-@unlink libnbody.so
	-@make -C ./lib/ clean
	-@rm -f $(EXPERIMENTS)
