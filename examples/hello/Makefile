# This is for CHARM_PATH
DOTS=../..
include $(DOTS)/Makefile.in


MDTPATH = $(DOTS)/src
OPTS = -g -O3
CXXFLAGS += $(OPTS)
LDFLAGS += $(OPTS) -L$(MDTPATH) -module mdt -module CkCache -language charm++  #-memory paranoid -tracemode projections 
CPPFLAGS = -I$(CHARM_PATH)/include -I$(MDTPATH)

CHARMC = $(CHARM_PATH)/bin/charmc $(INCPATH)

CXX = $(CHARMC) 
CC = $(CXX)
CXX_DEPEND = $(CXX) -M -MM -MG $(CXXFLAGS)

SRCS = Main.cpp 

TARGET = pgm
all: $(TARGET)

$(TARGET): pgm.def.h $(SRCS)
	$(CHARMC) $(CPPFLAGS) $(SRCS) $(LDFLAGS)  -o $(TARGET) 

%.decl.h %.def.h : %.ci
	$(CHARMC) -E $(OPTS) $(DEFINE_FLAGS) $<

%.o: Makefile

clean:
	rm -f core* $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 

projclean:
	rm -f *.gz *.projrc *.sts *.topo

depends:
	$(CXX_DEPEND) $(SRCS) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep

test: $(TARGET)
	./charmrun +p2 ./$(TARGET)

.PHONY: all docs dist clean depend test

include Makefile.dep
