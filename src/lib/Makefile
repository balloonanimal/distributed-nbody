CC = mpicc
CFLAGS_WARN = -Wall -Wvla -Werror -Wno-error=unused-variable
CFLAGS = -std=c11 -g -O2 -fPIC $(CFLAGS_WARN)

BUILD_DIR = build
_OBJ = simulation.o
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(_OBJ))
HEADERS = simulation.h types.h
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
