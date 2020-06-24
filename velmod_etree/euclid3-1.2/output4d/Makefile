CC = gcc
ETREE_DIR = ../libsrc
ETREE_LIB = $(ETREE_DIR)/libetree.a
CFLAGS = -g -Wall -I $(ETREE_DIR)

FILES = etree4d2txt bin2etree4d query4d

all: $(FILES)

query4d: query4d.c $(ETREE_LIB)

bin2etree4d: bin2etree4d.c $(ETREE_LIB)

etree4d2txt: etree4d2txt.c $(ETREE_LIB)

%.o: %.c 
	$(CC) -I$(INCDIR) -c $(CFLAGS) $< -o $@

clean:
	rm -f $(FILES) *.o *~

cleanall:
	rm -f $(FILES) *.e *.em *.o *~
