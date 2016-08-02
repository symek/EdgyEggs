# Minimal Automattes Makefile
INSTDIR = $(HIH)
SOURCES =  ./src/SOP_IGLDiscreteGeo.cpp
INCDIRS = -I$(HOME)/work/eigen -I$(HOME)/work/libigl/include
OPTIMIZER = -O3
DSONAME = IglSOP.so
# Include HDK Makefile.
include $(HT)/makefiles/Makefile.gnu	