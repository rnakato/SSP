CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS = -lz -lgsl -lgslcblas -lboost_thread
LIBS += -lboost_program_options -lboost_system -lboost_filesystem 

SRCDIR = ./src
OBJDIR = ./obj
ALGLIBDIR = ./src/alglib
BINDIR = ./bin

PROGRAMS = ssp 
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(OBJDIR)/readdata.o $(OBJDIR)/util.o $(OBJDIR)/BoostOptions.o
OBJS = $(OBJDIR)/ssp_main.o $(OBJDIR)/Mapfile.o $(OBJDIR)/ParseMapfile.o $(OBJDIR)/LibraryComplexity.o $(OBJDIR)/ssp_shiftprofile.o $(OBJDIR)/statistics.o $(ALGLIBDIR)/libalglib.a

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/ssp: $(OBJS) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LDFLAGS)

$(ALGLIBDIR)/libalglib.a:
	make -C $(ALGLIBDIR)

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin lib obj
	make clean -C $(ALGLIBDIR)

HEADS_UTIL = $(SRCDIR)/BoostOptions.hpp $(SRCDIR)/util.h $(SRCDIR)/readdata.h $(SRCDIR)/macro.h $(SRCDIR)/seq.h $(SRCDIR)/mthread.h $(SRCDIR)/SeqStats.hpp $(SRCDIR)/bpstatus.h

$(OBJDIR)/ParseMapfile.o: $(SRCDIR)/ssp_shiftprofile.h
$(OBJDIR)/ssp_shiftprofile.o: Makefile $(SRCDIR)/ssp_shiftprofile_p.h $(SRCDIR)/ssp_shiftprofile.h
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS): Makefile $(SRCDIR)/ssp_gv.h $(SRCDIR)/Mapfile.hpp $(SRCDIR)/ParseMapfile.hpp $(SRCDIR)/LibraryComplexity.hpp $(SRCDIR)/statistics.h $(HEADS_UTIL)
