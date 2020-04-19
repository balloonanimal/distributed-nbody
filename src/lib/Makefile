CC = mpicc
CFLAGS_WARN = -Wall -Wvla -Werror -Wno-error=unused-variable
CFLAGS = -std=c99 -g -O2 -fPIC -lm $(CFLAGS_WARN)

BUILD_DIR = build
_OBJ = simulation.o mpi_utils.o
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(_OBJ))
HEADERS = simulation.h types.h mpi_utils.h
_LIBFILE = libnbody.so
LIBFILE = $(BUILD_DIR)/$(_LIBFILE)


all: $(BUILD_DIR) $(LIBFILE)

$(BUILD_DIR):
	mkdir $@

$(LIBFILE): $(OBJ)
	$(CC) $(CFLAGS) -shared $(OBJ) -o $@

$(BUILD_DIR)/%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-@rm -f $(BUILD_DIR)/*.o
	-@rm -f $(BUILD_DIR)/*.so
