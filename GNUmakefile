name := SPL_SUPERBEAM
G4TARGET := $(name)
G4EXLIB := true

.PHONY: all
all: lib bin

ROOTCFLAGS    = $(shell root-config --cflags)
CPPFLAGS += -I$(ROOTSYS)/
CPPFLAGS += -I$/usr/include
CPPFLAGS += -I/usr/include/root
CPPFLAGS += $(ROOTCFLAGS)

CPPFLAGS += -I$(ROOTSYS)/include
EXTRALIBS = $(shell root-config --glibs)

include $(G4INSTALL)/config/analysis.gmk
include $(G4INSTALL)/config/architecture.gmk
include $(G4INSTALL)/config/binmake.gmk

#ANALYSISLIBS += `aida-config --libs`
#CPPFLAGS += `aida-config --include`
#LIBS += `aida-config --lib`
#LIBS += `root-config --libs`
LIBS += `root-config --cflags --glibs`

SBDetectorConstructionCint.cc: $(CINT)
	@echo "Generating dictionary $@..."
	/home/irfulx114/mnt/alonghin/work/root/bin/rootcint -f $@ -c -p $(INCFLAGS) ./include/SBDetectorConstruction.hh $^
#SBDetectorConstructionCintLinkDef.hh $^

MyClassCint.cc: $(CINT)
	@echo "Generating dictionary $@..."
	rootcint -f ./src/$@ -c -p -I/home/irfulx114/mnt/alonghin/work/EUROnu/sim/g4/project_Lyon/project/trunk/include MyClass.hh $^

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

histclean:
	rm ${G4WORKDIR}/tmp/${G4SYSTEM}/${G4TARGET}/HistoManager.o
#include TargetsDef.mk
#rootcint -f SBDetectorConstructionCint.C -c -p -I/home/irfulx114/mnt/alonghin/work/EUROnu/sim/g4/v2.2/include -I/home/irfulx114/mnt/alonghin/work/GEANT/geant4.9.2/source/run/include -I/home/irfulx114/mnt/alonghin/work/GEANT/geant4.9.2/source/global/management/include SBDetectorConstruction.hh SBDetectorConstructionLinkDef.h
#g++ -shared SBDetectorConstruction.o  SBDetectorConstructionCint.o -o SBDetectorContruction.so
