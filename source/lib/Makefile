# compiler
CC = mpicc
CFLAGS_WARN = -Wall -Wvla -Werror -Wno-error=unused-variable
CFLAGS = -std=c99 -g -fPIC $(CFLAGS_WARN)

LDLIBS = -lm

# OS specific sharedlib
OS = $(shell uname)
ifeq ($(OS),Darwin)
	LDFLAGS = -dynamiclib
endif
ifeq ($(OS),Linux)
	LDFLAGS = -shared
endif

# build targets
SRC = $(wildcard *.c)
_OBJ = $(SRC:.c=.o)
DEP = $(wildcard *.h)

BUILD_DIR = build
OBJ = $(patsubst %,$(BUILD_DIR)/%,$(_OBJ))

_LIBFILE = libnbody.so
LIBFILE = $(BUILD_DIR)/$(_LIBFILE)


all: $(BUILD_DIR) $(LIBFILE)

$(BUILD_DIR):
	mkdir $@

$(LIBFILE): $(OBJ)
	$(CC) $(CFLAGS) $(LDLIBS) $(LDFLAGS) $(OBJ) -o $@

$(BUILD_DIR)/%.o: %.c $(DEP)
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	-@rm -f $(BUILD_DIR)/*.o
	-@rm -f $(BUILD_DIR)/*.so
