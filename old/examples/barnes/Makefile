DOTS=../..
include ../../Makefile.in

DEFINE_FLAGS = -DSTATIC #-DDEBUG_TRAVERSAL #-DLOW_OVERHEAD_DEBUG_TRAVERSAL
ORB3DLB_PATH = ./Orb3dLB
INPUTFILE_PATH= ./inputfiles
VPATH = $(STRUCTURES_PATH)

INCPATH = -I$(DOTS)/src -I$(ORB3DLB_PATH) #-I$(CHARM_PATH)/tmp/libs/ck-libs/NDMeshStreamer
OPTS = -g -O3
CXXFLAGS += $(OPTS) $(DEFINE_FLAGS) -I$(STRUCTURES_PATH)
LDFLAGS += $(OPTS) -L$(STRUCTURES_PATH) -L. -L$(CHARM_PATH)/lib -L$(DOTS)/src -L$(ORB3DLB_PATH) -lTipsy -module mdt -module CkCache -language charm++ -module completion -balancer RotateLB -balancer Orb3dLB #-memory paranoid -module NDMeshStreamer -tracemode projections 
CPPFLAGS = -I$(CHARM_PATH)/include -I$(STRUCTURES_PATH)

CHARMC = $(CHARM_PATH)/bin/charmc $(INCPATH)

CXX = $(CHARMC) 
CC = $(CXX)
AR = ar qrv 
CXX_DEPEND = $(CXX) -M -MM -MG $(CXXFLAGS)

OBJECTS = Main.o \
          TreePiece.o \
          Visitor.o \
          Decomposer.o \
          BoundingBox.o \
          Particle.o \
          ParticleMsg.o \
          NodePayload.o \
          TraversalData.o \
          Splitter.o \
          SplitterGroup.o \
          LbOnOffGroup.o \
          Environment.o \
          ParamStorage.o

SRCS = Main.cpp \
       TreePiece.cpp \
       Visitor.cpp \
       Decomposer.cpp \
       BoundingBox.cpp \
       Particle.cpp \
       ParticleMsg.cpp \
       NodePayload.cpp \
       TraversalData.cpp \
       Splitter.cpp \
       SplitterGroup.cpp \
       LbOnOffGroup.cpp \
       Environment.cpp \
       ParamStorage.cpp \

TARGET = bh
all: $(ORB3DLB_PATH)/libmoduleOrb3dLB.a $(STRUCTURES_PATH)/libTipsy.a $(INPUTFILE_PATH)/tipsyPlummer $(TARGET)

$(TARGET): $(ORB3DLB_PATH)/libmoduleOrb3dLB.a $(STRUCTURES_PATH)/libTipsy.a $(OBJECTS)
	$(CHARMC) $(LDFLAGS) -o $(TARGET) $(OBJECTS)

$(INPUTFILE_PATH)/tipsyPlummer:
	cd $(INPUTFILE_PATH); $(MAKE)

$(STRUCTURES_PATH)/libTipsy.a: $(STRUCTURES_PATH)/Makefile
	cd $(STRUCTURES_PATH); $(MAKE)

$(STRUCTURES_PATH)/Makefile:
	cd $(STRUCTURES_PATH); ./configure

$(ORB3DLB_PATH)/libmoduleOrb3dLB.a:
	cd $(ORB3DLB_PATH); $(MAKE) libmoduleOrb3dLB.a

%.decl.h %.def.h : %.ci
	$(CHARMC) -E $(OPTS) $(DEFINE_FLAGS) $<

%.o: Makefile

clean:
	rm -f core* $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 

projclean:
	rm -f *.gz *.projrc *.sts *.topo

depends:
	$(CXX_DEPEND) $(SRCS) | while read i;do echo $$i| awk -F' ' '{for (i=1;i<NF;++i) print $$i" \\"}';echo;done|grep -v "$(CHARM_PATH)/bin" > Makefile.dep

test: bh $(INPUTFILE_PATH)/tipsyPlummer
	./charmrun +p2 ./bh --nPieces 33 --nIterations 20 --maxppc 500 --inFile ./inputfiles/1k.tipsy +balancer Orb3dLB ++ppn 1

.PHONY: all docs dist clean depend test

include Makefile.dep
