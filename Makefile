# Minimal Makefile
INSTDIR = $(HIH)
SOURCES =  ./src/converters.cpp\
   ./src/SOP_IGLUVproject.cpp\
   ./src/SOP_IGLDiscreteGeo.cpp\
   ./src/SOP_ShapeOp.cpp\
   ./src/SOP_IGLMain.cpp

INCDIRS = -I$(HOME)/work/eigen\
	-I$(HOME)/work/libigl/include\
	-I$(HOME)/work/ShapeOp\
	-I$(HOME)/work/ShapeOp/libShapeOp/src

OPTIMIZER = -O3
DSONAME = IglSOP.so
# Include HDK Makefile.
include $(HT)/makefiles/Makefile.gnu
