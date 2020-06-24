# -*- Makefile -*-

# 
# Uncomment this if compiling on an Alpha/DUX40D platform
#
#ALPHAFLAGS = -DALPHA -DALIGNMENT -D_STACK_T


##################################################
# You shouldn't need to change anything below here
##################################################

CC = gcc

ETREE_DIR = $(CURDIR)/../libsrc
ETREE_LIB = $(ETREE_DIR)/libetree.so


# Compiler flags
CFLAGS = -g -Wall -I $(ETREE_DIR) $(ALPHAFLAGS)

TARGET = tiny query txt2etree etree2txt refine 

all: $(TARGET)

tiny: tiny.o $(ETREE_LIB)
query: query.o $(ETREE_LIB)
txt2etree: txt2etree.o $(ETREE_LIB)
etree2txt: etree2txt.o traverse.o $(ETREE_LIB)
refine: refine.o traverse.o $(ETREE_LIB)


test: all
	./tiny tree.e
	./txt2etree tree.e < tree.txt
	./etree2txt tree.e
	./refine tree.e
	./query tree.e <tree.txt

clean:
	rm -f *.o *~

cleanall:
	rm -f $(TARGET) *.o *~
