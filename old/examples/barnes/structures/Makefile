FLAG_BIGKEYS = 
CXXFLAGS += -O3  -I.. $(FLAG_BIGKEYS)
CFLAGS += -O3  -I.. $(FLAG_BIGKEYS)
CXX = g++ -m64  -m64 -fPIC  
CC = gcc -m64  -m64 -fPIC  
AR = ar 

NAME = structures
VERSION = 1.1
INSTALL_LOCATION = $(HOME)/my_stuff/bin/
LIBTARGETS = TipsyReader.o TipsyFile.o SS.o NChilReader.o hilbert.o
EXECTARGETS = browse see_header see_treexdr random_sample tipsy2ifrit ss2ifrit
LDFLAGS += -L../libs
LDLIBS += 

.PHONY: all docs dist clean

all: libTipsy.a($(LIBTARGETS)) $(EXECTARGETS)  SFC.o

libTipsy.a: $(LIBTARGETS)
	-rm -f libTipsy.a
	$(AR) libTipsy.a $(LIBTARGETS)
	ranlib libTipsy.a

TipsyFile.o: TipsyFile.h TipsyParticles.h TipsyReader.h TipsyFile.cpp \
             Vector3D.h OrientedBox.h

TipsyReader.o: TipsyReader.h TipsyReader.cpp TipsyParticles.h Vector3D.h xdr_template.h
NChilReader.o: NChilReader.h NChilReader.cpp tree_xdr.h Vector3D.h xdr_template.h

SS.o: xdr_template.h Vector3D.h

# SFC.o is only used by TreeDataFormat/tipsy2tree:
SFC.o: SFC.cpp SFC.h
	$(CXX) -O -ffloat-store -c SFC.cpp

$(EXECTARGETS): libTipsy.a

docs:
	doxygen structures.doxygen

dist:
	mkdir $(NAME)-$(VERSION)
	cp Makefile *.h *.cpp structures.doxygen $(NAME)-$(VERSION)/
	tar zcf $(NAME)-$(VERSION).tar.gz $(NAME)-$(VERSION)
	rm -Rf $(NAME)-$(VERSION)

clean:
	rm -f core* *.o *~ $(TARGETS) libTipsy.a

install: all
	cp $(EXECTARGETS) $(INSTALL_LOCATION)
