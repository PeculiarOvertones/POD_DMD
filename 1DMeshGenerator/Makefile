# I am a comment, and I want to say that the variable CC will be
# the compiler to use.
#CC=icpc
#CC=CC 
CC=g++
ifneq (,$(findstring DDebug,$(COMPFLAGS)))
CFLAGS=-c 
#CFLAGS=-c -std=gnu++0x -Wno-deprecated
else
CFLAGS=-c 
endif 
GFLAGS= 
PLATFORM=BW
DEBUG=none
OPTIMIZATION=none


CFLAGS=-c -std=gnu++0x 

# Merge all the flags

INCLINK=  $(VTK_INCLINK)

all: maestro.x Makefile
	rm -rf *.o

ifeq ($(PLATFORM),lionx)
	echo Lionx
else
	echo "Somewhere else"
endif

maestro.x: outputMesh.o amrMesh.o Makefile
	$(CC) $(LIBLINK) $(GFLAGS) outputMesh.o amrMesh.o -o maestro.x $(LIBLINK2)
	
amrMesh.o: amrMesh.cpp Makefile
	$(CC) $(COMPFLAGS) $(CFLAGS) amrMesh.cpp

outputMesh.o: outputMesh.cpp Makefile
	$(CC) $(COMPFLAGS) $(CFLAGS) outputMesh.cpp

clean:
	rm -rf *.o

veryclean:
	rm -rf *.o *.x
