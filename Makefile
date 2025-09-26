CC = cc

SRCDIR = src
BLDDIR = build
BINDIR = bin

SRC = $(wildcard $(SRCDIR)/*.c)
OBJ = $(SRC:$(SRCDIR)/%.c=$(BLDDIR)/%.o)

CFLAGS = -O3 -Wall -Wextra -Wno-unused-but-set-parameter -Wno-missing-braces -Wno-unused-parameter -funroll-loops -march=native
LDFLAGS = -lm

TARGET = $(BINDIR)/a

all: $(TARGET)

$(TARGET): $(OBJ) | $(BINDIR)
	$(CC) $(LDFLAGS) $^ -o $@
$(OBJ): $(SRC) | $(BLDDIR)
	$(CC) $(CFLAGS) -c $(@:$(BLDDIR)/%.o=$(SRCDIR)/%.c) -o $@

$(BLDDIR) $(BINDIR):
	mkdir $@

run: $(TARGET)
	@./$(TARGET) > image.ppm

clean:
	rm -rf $(BLDDIR)