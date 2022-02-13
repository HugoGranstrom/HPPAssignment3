CC = gcc
LD = gcc
CFLAGS = -Wall
LDFLAGS = 
LDLIBS = -lm
RM = /bin/rm -f
OBJS = galsim.o graphics/graphics.o 
SORT = galsim 
all:$(SORT)

$(SORT): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(SORT)

graphics.o: graphics/graphics.c graphics/graphics.h 
	$(MAKE) -C ./graphics

galsim.o: galsim.c
	$(CC) $(CFLAGS) -c galsim.c
clean:
	$(RM) $(SORT) $(OBJS)
