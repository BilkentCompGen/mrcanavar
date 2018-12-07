CC=gcc
CFLAGS = -c -O3 -g -Wall
LDFLAGS = -lz -lm
SOURCES = mrcanavar.c utils.c prep.c sam.c callcnv.c gcnorm.c
OBJECTS = $(SOURCES:.c=.o)
AUTOSOURCES = mrcanavar-auto.c utils.c
AUTOOBJ = $(AUTOSOURCES:.c=.o)
EXECUTABLE = mrcanavar
AUTOEXEC = mrcanavar-auto
BIN = /usr/local/bin/

all: $(SOURCES) $(EXECUTABLE) $(AUTOEXEC)
	rm -rf *.o

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

$(AUTOEXEC): $(AUTOOBJ) 
	$(CC) $(AUTOOBJ) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@

auto: utils.c mrcanavar-auto.c
	$(CC) $(CFLAGS) utils.c
	$(CC) mrcanavar-auto.c utils.o -o mrcanavar-auto $(LDFLAGS)

clean: 
	rm -f $(EXECUTABLE) $(AUTOEXEC) *.o *~ \#* 

install:	
	cp $(EXECUTABLE) $(AUTOEXEC) $(BIN)
