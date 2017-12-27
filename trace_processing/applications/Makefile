
# Build environment can be configured the following
# environment variables:
#   CC : Specify the C compiler to use
#   CFLAGS : Specify compiler options to use

#  IMPORTANT: change here to set the directory for binary executable files
ifdef MYBIN
BINDIR=${MYBIN}
else
# with the following, binary executables will be placed in bin subdirectory of your home directory
#BINDIR=~/bin/
# if in doubt, use the following - binary executables will be placed in the current directory
BINDIR=.
endif

# Options specific for GCC
ifndef CC
CC = gcc
endif
#
ifndef CCFLAGS
#
CCFLAGS_BASIC =  -Wall -fPIC -I..
#
# optimized
CCFLAGS = -O3 $(CCFLAGS_BASIC)
#
# profile
#CCFLAGS=-O3 -pg $(CCFLAGS_BASIC)
#
# debug - gdb, valgrind, ...
#CCFLAGS = $(CCFLAGS_BASIC) -g
# valgrind --leak-check=yes  exe_name <args>
# valgrind --leak-check=full --show-reachable=yes exe_name exe_name <args>
endif

LDFLAGS = $(LDFLAGS_RMQ) -L..
LDLIBS = $(LDLIBS_RMQ) -lm -lmseed -ltdproc -lm `xml2-config --cflags` `xml2-config --libs`

OBJ= $(LIB_OBJS) app_lib.o



all: $(BINDIR)/seedlink_monitor${EXE_NAME} $(BINDIR)/miniseed_process${EXE_NAME} $(BINDIR)/pick_process${EXE_NAME}

$(BINDIR)/pick_process${EXE_NAME}: pick_process.o $(OBJ)
	$(CC) $(CCFLAGS) -o $(BINDIR)/pick_process${EXE_NAME} pick_process.o $(OBJ) $(LDFLAGS) $(LDLIBS)


$(BINDIR)/miniseed_process${EXE_NAME}: miniseed_process.o mseedrtstream_lib.o $(OBJ)
	$(CC) $(CCFLAGS) -o $(BINDIR)/miniseed_process${EXE_NAME} miniseed_process.o mseedrtstream_lib.o $(OBJ) $(LDFLAGS) $(LDLIBS)


$(BINDIR)/seedlink_monitor${EXE_NAME}: seedlink_monitor.o $(OBJ)
	$(CC) $(CCFLAGS) -o $(BINDIR)/seedlink_monitor${EXE_NAME} seedlink_monitor.o $(OBJ) $(LDFLAGS) $(LDLIBS) -lslink

clean:
	rm -f $(BINDIR)/pick_process${EXE_NAME} $(BINDIR)/miniseed_process${EXE_NAME} $(BINDIR)/seedlink_monitor${EXE_NAME}
	rm -f *.o
	rm -f *.a


# Implicit rule for building object files
%.o: %.c
	$(CC) -c $< $(CCFLAGS)

install:
	@echo
	@echo "No install target, copy the executable(s) yourself"
	@echo
