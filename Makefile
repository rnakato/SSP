CC = clang++
CFLAGS  = -std=c++11 -O2 -Wall -W
LDFLAGS = -lz -lgsl -lgslcblas -lboost_thread
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lpthread

SRCDIR = ./src
CMNDIR = ./common
OBJDIR = ./obj
CMNOBJDIR = ./cobj
BINDIR = ./bin

PROGRAMS = ssp
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))

ifdef CLOCK
CFLAGS += -DCLOCK
endif
ifdef DEBUG
CFLAGS += -DDEBUG
endif
ifdef PRINTREAD
CFLAGS += -DPRINTREAD
endif

OBJS = $(OBJDIR)/ssp_main.o $(OBJDIR)/Mapfile.o $(OBJDIR)/ParseMapfile.o $(OBJDIR)/LibraryComplexity.o $(OBJDIR)/ShiftProfile.o $(OBJDIR)/FragmentClusterScore.o
OBJS += $(CMNOBJDIR)/statistics.o $(CMNOBJDIR)/util.o $(CMNOBJDIR)/BoostOptions.o $(CMNOBJDIR)/gzstream.o

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/ssp: $(OBJS)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LDFLAGS)

$(CMNOBJDIR)/gzstream.o: $(CMNDIR)/gzstream.C $(CMNDIR)/gzstream.h
	$(CC) -o $@ -c $< -I. -O

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS)

$(CMNOBJDIR)/%.o: $(CMNDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf $(BINDIR) $(OBJDIR) $(CMNOBJDIR)

HEADS = $(SRCDIR)/ssp_gv.hpp $(SRCDIR)/Mapfile.hpp $(SRCDIR)/ParseMapfile.hpp $(SRCDIR)/LibraryComplexity.hpp $(CMNDIR)/BoostOptions.hpp $(SRCDIR)/MThread.hpp $(SRCDIR)/SeqStats.hpp $(CMNDIR)/BedFormat.hpp $(SRCDIR)/FragmentClusterScore.hpp
HEADS += $(CMNDIR)/inline.hpp $(CMNDIR)/seq.hpp $(CMNDIR)/statistics.hpp $(CMNDIR)/util.hpp

$(OBJDIR)/ParseMapfile.o: $(SRCDIR)/ShiftProfile.hpp
$(OBJDIR)/ShiftProfile.o: $(SRCDIR)/ShiftProfile_p.hpp $(SRCDIR)/ShiftProfile.hpp
$(OBJDIR)/FragmentCluterScore.o: $(SRCDIR)/FragmentCluterScore_p.hpp
$(OBJS): Makefile $(HEADS)
