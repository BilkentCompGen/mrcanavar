CC=gcc
CFLAGS = -c -O3 -g -Wall
LDFLAGS = -lz -lm
SOURCES = mrcanavar.c utils.c prep.c sam.c callcnv.c gcnorm.c
OBJECTS = $(SOURCES:.c=.o)
EXECUTABLE = mrcanavar
BIN = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -f $(EXECUTABLE) *.o *~ \#* 

install:	
	cp mrcanavar $(BIN)
