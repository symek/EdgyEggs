# Minimal Automattes Makefile
INSTDIR = $(HIH)
SOURCES =  ./src/SOP_IGLDiscreteGeo.cpp
INCDIRS = $(HOME)/work/eigen $(HOME)/work/libigl/include
OPTIMIZER = -O3
DSONAME = libIglSOP.so
# Include HDK Makefile.
include $(HT)/makefiles/Makefile.gnu	