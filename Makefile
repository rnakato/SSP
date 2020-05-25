.PHONY: all clean

HTSLIBDIR=src/htslib-1.10.2/

all: bin/ssp $(HTSLIBDIR)/libhts.a

bin/ssp: $(HTSLIBDIR)/libhts.a
	mkdir -p build && cd build && cmake .. && make
	mkdir -p bin
	cp build/test/ssp bin

$(HTSLIBDIR)/libhts.a:
	$(MAKE) -C $(HTSLIBDIR)

clean:
	rm -rf build bin
	make -C $(HTSLIBDIR) clean
