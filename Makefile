CC = g++
CFLAGS += -std=c++11 -O2 -Wall -W
LDFLAGS =
LIBS += -lboost_program_options -lboost_system -lboost_filesystem -lboost_iostreams
LIBS_DP += -lz -lgsl -lgslcblas -lboost_thread

SRCDIR = ./src
OBJDIR = ./obj
LIBDIR = ./lib
BINDIR = ./bin
ALGLBDIR = $(SRCDIR)/alglib-3.10.0/src

PROGRAMS = ssp 
TARGET = $(addprefix $(BINDIR)/,$(PROGRAMS))
#$(warning $(TARGET))

ifdef DEBUG
CFLAGS += -DDEBUG
endif

OBJS_UTIL = $(OBJDIR)/readdata.o $(OBJDIR)/util.o
OBJS_PW = $(OBJDIR)/ssp_main.o $(OBJDIR)/ssp_estFlen.o $(OBJDIR)/pw_readmapfile.o $(OBJDIR)/pw_makefile.o $(OBJDIR)/pw_gc.o $(OBJDIR)/pw_shiftprofile.o $(OBJDIR)/statistics.o $(LIBDIR)/libalglib.a

.PHONY: all clean

all: $(TARGET)

$(BINDIR)/ssp: $(OBJS_PW) $(OBJS_UTIL)
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ $^ $(LIBS) $(LIBS_DP)

$(LIBDIR)/libalglib.a: $(SRCDIR)/alglib.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -c $< $(ALGLBDIR)/specialfunctions.cpp $(ALGLBDIR)/ap.cpp $(ALGLBDIR)/alglibinternal.cpp $(CFLAGS)
	ar r $@ alglib.o alglibinternal.o ap.o specialfunctions.o
	rm alglib.o alglibinternal.o ap.o specialfunctions.o

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@if [ ! -e `dirname $@` ]; then mkdir -p `dirname $@`; fi
	$(CC) -o $@ -c $< $(CFLAGS) $(WFLAGS)

clean:
	rm -rf bin lib obj

HEADS_UTIL = $(SRCDIR)/util.h $(SRCDIR)/readdata.h $(SRCDIR)/macro.h $(SRCDIR)/seq.h

$(LIBDIR)/libalglib.a: Makefile $(SRCDIR)/alglib.h
$(OBJDIR)/ssp_main.o: $(SRCDIR)/pw_makefile.h $(SRCDIR)/pw_gc.h
$(OBJDIR)/pw_readmapfile.o: $(SRCDIR)/pw_shiftprofile.h
$(OBJDIR)/pw_makefile.o: $(SRCDIR)/pw_makefile.h
$(OBJDIR)/pw_gc.o: $(SRCDIR)/pw_gc.h
$(OBJDIR)/pw_shiftprofile.o: Makefile $(SRCDIR)/pw_shiftprofile_p.h $(SRCDIR)/pw_shiftprofile.h
$(OBJS_UTIL): Makefile $(HEADS_UTIL)
$(OBJS_PW): Makefile $(SRCDIR)/pw_gv.h $(SRCDIR)/pw_readmapfile.h $(SRCDIR)/statistics.h $(HEADS_UTIL) $(SRCDIR)/ssp_estFlen.h
