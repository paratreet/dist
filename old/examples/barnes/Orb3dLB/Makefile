DOTS=../../..
include ../../../Makefile.in

VPATH = $(STRUCTURES_PATH)

INCPATH = -I$(STRUCTURES_PATH) 
OPTS = -g -O3
CXXFLAGS += $(OPTS) $(DEFINE_FLAGS) -I$(STRUCTURES_PATH)
LDFLAGS += $(OPTS) -L$(STRUCTURES_PATH) -L. -L$(CHARM_PATH)/lib -language charm++ 

CHARMC = $(CHARM_PATH)/bin/charmc $(INCPATH)

CXX = $(CHARMC) 
CC = $(CXX)
AR = ar q 
CXX_DEPEND = $(CXX) -M -MM -MG $(CXXFLAGS)

OBJECTS = Orb3dLB.o

SRCS = Orb3dLB.cpp

TARGET = libmoduleOrb3dLB.a
all: $(TARGET) 

$(TARGET): $(OBJECTS) 
	$(CHARMC) $(LDFLAGS) -o $(TARGET) $(OBJECTS) 

%.decl.h %.def.h : %.ci
	$(CHARMC) -E $(OPTS) $(DEFINE_FLAGS) $<

%.o: Makefile

clean:
	rm -f core* $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 

projclean:
	rm -f *.gz *.projrc *.sts *.topo

depends:
	$(CXX_DEPEND) $(SRCS) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep

.PHONY: all docs dist clean depend test

include Makefile.dep
