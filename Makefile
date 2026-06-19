.PHONY: all configure build clean

HTSLIBDIR = src/htslib-1.10.2
BUILDDIR  = build
BINDIR    = bin

ifdef DEBUG
CMAKEFLAGS += -DENABLE_DEBUG=ON
endif

all: build

configure: $(HTSLIBDIR)/libhts.a
	cmake -S . -B $(BUILDDIR) $(CMAKEFLAGS)

bin/ssp: $(HTSLIBDIR)/libhts.a
	mkdir -p build
	cd build && cmake $(CMAKEFLAGS) ..
	$(MAKE) -C build
	mkdir -p bin
	cp build/test/ssp bin

build: configure
	cmake --build $(BUILDDIR)
	mkdir -p $(BINDIR)
	cp $(BUILDDIR)/test/ssp $(BINDIR)/ssp

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf $(BUILDDIR) $(BINDIR)
	$(MAKE) -C $(HTSLIBDIR) clean
