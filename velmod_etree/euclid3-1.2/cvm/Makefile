#
# Etree directory and library
#
ETREE_DIR = $(CURDIR)/../libsrc
ETREE_LIB = $(ETREE_DIR)/libetree.so


#
# Set the following variable to the absolute pathname of the
# destination directory in which to to store the CVM etree
#
DESTDIR = /usr1/tutk/scratch

#
# Uncomment this if compiling on an Alpha/DUX40D platform
#
#ALPHAFLAGS = -DDUX40D  
#


#
# Uncommnet this if you want to measure performance 
#
PERFFLAGS = -DPERF

##################################################
# You shouldn't need to change anything below here
##################################################

#
# We use GCC to compile this program
#
CC = gcc
LD = gcc


# Compiler flags Choose one Compiler flags
#CFLAGS = -O2 -Wall -I $(ETREE_DIR) -D_LARGEFILE64_SOURCE \
	-D_FILE_OFFSET_BITS=64 \
	$(ALPHAFLAGS) $(PERFFLAGS) \
	-DENABLE_OUTPUT -DOUTPUT_BINARY -DREPORT_NAN
CFLAGS = -g -Wall -I $(ETREE_DIR) $(ALPHAFLAGS) $(PERFFLAGS)

# Linking flags
LDFLAGS = -lm 


TARGET = definedomain createcvm packcvm querycvm 


all: $(TARGET)

definedomain: definedomain.o domain.o
createcvm: createcvm.o domain.o $(ETREE_LIB)
packcvm: packcvm.o $(ETREE_LIB)
querycvm: querycvm.o $(ETREE_LIB) 

clean: 
	rm -f $(TARGET) *.o *~ 

cleanquery:
	rm -f *.txt

# 
# create CVM etree for LA basin
#
MODELDIR = $(DESTDIR)/model

# Parameters to createcvm 
DOMAIN = labase.domain
LAYERS = labase.2 #labase.layers
ETREENAME_UNPACKED = lasmall_unpacked.e
ETREENAME_PACKED = lasmall.e

cvmetree: 
	rm -rf $(MODELDIR)
	mkdir $(MODELDIR)
	tar zxvf Version3.0.tar.gz -C $(MODELDIR)
	g77 $(MODELDIR)/version3.0.f -o $(MODELDIR)/version3.0
	clear
	./createcvm $(DOMAIN) $(LAYERS) $(DESTDIR) $(ETREENAME_UNPACKED)
	./packcvm $(DESTDIR)/$(ETREENAME_UNPACKED) $(DESTDIR)/$(ETREENAME_PACKED)
	rm -rf $(MODELDIR)
	rm -f $(DESTDIR)/$(ETREENAME_UNPACKED)


#
# Examples of querying cvm etree
#
pointquery:
	./querycvm $(DESTDIR)/labase.e -118.5 33.6 0 -118.5 33.6 0 0 0 Vs point.txt

linequery:
	./querycvm $(DESTDIR)/labase.e -118.5 33.6 0 -118.25 33.8 1000 1 100 Vp line.txt

slicequery:
	./querycvm $(DESTDIR)/labase.e -118.5 33.6 0 -118.25 33.8 1000 2 100 all slice.txt

volumequery:
	./querycvm $(DESTDIR)/labase.e -118.5 33.6 0 -118.25 33.8 1000 3 100 density volume.txt


# 
# Performance measurement 
#
performance:
# 1D lines 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 1000 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 500 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 200 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 150 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 100 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 50 density $(DESTDIR)/line.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 1 20 density $(DESTDIR)/line.txt 


#2D slices 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 1000 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 500 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 200 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 150 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 100 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 50 density $(DESTDIR)/slice.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 2 20 density $(DESTDIR)/slice.txt 

#3D slices 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 3 1000 density $(DESTDIR)/volume.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 3 500 density $(DESTDIR)/volume.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 3 200 density $(DESTDIR)/volume.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 3 150 density $(DESTDIR)/volume.txt 
	./querycvm $(DESTDIR)/labase.e -118.68 33.6 0 -118.18 34.2 30000 3 100 density $(DESTDIR)/volume.txt 
