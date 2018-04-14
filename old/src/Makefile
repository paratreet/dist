DOTS=..
include ../Makefile.in

DEFINE_FLAGS = #-DUSE_MESH_STREAMER

INCPATH = -I.. #-I$(CHARM_PATH)/tmp/libs/ck-libs/NDMeshStreamer 
OPTS = -g -O3
CXXFLAGS += $(OPTS) $(DEFINE_FLAGS) 
LDFLAGS += $(OPTS) -L$(STRUCTURES_PATH) -L. -L$(CHARM_PATH)/lib -language charm++ 

CHARMC = $(CHARM_PATH)/bin/charmc $(INCPATH)

CXX = $(CHARMC) 
CC = $(CXX)
AR = ar q 
CXX_DEPEND = $(CXX) -M -MM -MG $(CXXFLAGS)

OBJECTS = mdt.o

SRCS = mdt.cpp

TARGET = libmodulemdt.a
all: $(TARGET) 

$(TARGET): $(OBJECTS) 
	$(CHARMC) $(LDFLAGS) -o $(TARGET) $(OBJECTS) 

%.decl.h %.def.h : %.ci
	$(CHARMC) -E $(OPTS) $(DEFINE_FLAGS) $<

%.o: Makefile

clean:
	rm -f core* $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h 

depends:
	$(CXX_DEPEND) $(SRCS) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep

test: $(TARGET)
	cd $(BARNES_PATH); $(MAKE) test

.PHONY: all docs dist clean depend test

include Makefile.dep
