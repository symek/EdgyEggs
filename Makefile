# Minimal Makefile
INSTDIR = $(HIH)
SOURCES =  ./src/converters.cpp\
   ./src/SOP_IGLUVproject.cpp\
   ./src/SOP_IGLDiscreteGeo.cpp\
   ./src/SOP_ShapeOp.cpp\
   ./src/SOP_Eltopo.cpp\
   ./src/SOP_IGLMain.cpp

INCDIRS = -I$(HOME)/work/eigen\
	-I$(HOME)/work/libigl/include\
	-I$(HOME)/work/ShapeOp\
	-I$(HOME)/work/eltopo

LIBS = $(HOME)/work/eltopo/eltopo3d/libeltopo_release.a -llapack -lblas

OPTIMIZER = -O3
DSONAME = SOP_EdgyEggs.so
# Include HDK Makefile.
include $(HT)/makefiles/Makefile.gnu
