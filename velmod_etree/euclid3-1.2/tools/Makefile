#
# Etree directory and library
#
ETREE_DIR = $(CURDIR)/../libsrc
ETREE_LIB = $(ETREE_DIR)/libetree.so

#
# We use GCC to compile this program
#
CC = gcc
LD = gcc

# Compiler flags
DEBUGFLAGS = -g

CFLAGS = $(DEBUGFLAGS) -Wall -gstabs+3\
	-I$(ETREE_DIR) -D_LARGEFILE64_SOURCE -D_FILE_OFFSET_BITS=64 \

LDFLAGS = -lm

-include user.make

TARGET = findpixel etreeinfo

all: $(TARGET)

findpixel: findpixel.o $(ETREE_LIB)
etreeinfo: etreeinfo.o $(ETREE_LIB)

clean: 
	rm -f core *.o *~

cleanall:
	rm -f $(TARGET) *.o core *~


