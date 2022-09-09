#directories
INCLUDEDIR = include
SRCDIR = src
OBJDIR = obj
EXECUTABLEDIR = bin

# Compiler
CC = g++
CFLAGS = -c -g -Wall `root-config --cflags` -I$(INCLUDEDIR)

# Linker
LINKER = g++
LDFLAGS = `root-config --libs`

SOURCES := $(wildcard $(SRCDIR)/*.cc)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
TARGETS_SOURCES := $(wildcard $(SRCDIR)/*.cpp)
TARGETS_OBJECTS := $(TARGETS_SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
TARGETS := $(TARGETS_SOURCES:$(SRCDIR)/%.cpp=$(EXECUTABLEDIR)/%)

all: $(TARGETS)

$(TARGETS): $(EXECUTABLEDIR)/%: $(OBJDIR)/%.o $(OBJECTS)
	@echo "\tLINKING "$@
	@$(LINKER) $< $(OBJECTS) $(LDFLAGS) -o $@

$(TARGETS_OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo "\tCOMPILING "$<
	@$(CC) $(CFLAGS) $< -o $@

$(OBJECTS): $(OBJDIR)/%.o: $(SRCDIR)/%.cc
	@echo "\tCOMPILING "$<
	@$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(TARGETS) $(wildcard $(OBJDIR)/*)


