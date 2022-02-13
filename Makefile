CC = gcc
LD = gcc
CFLAGS = -Wall -O3
LDFLAGS = 
LDLIBS = -lm -L/opt/X11/lib -lX11 -lm
RM = /bin/rm -f
OBJS = galsim.o graphics/graphics.o 
SORT = galsim 
all:$(SORT)

$(SORT): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LDLIBS) -o $(SORT)

graphics.o: graphics/graphics.c graphics/graphics.h 
	$(MAKE) -C ./graphics -f Makefile

galsim.o: galsim.c
	$(CC) $(CFLAGS) -c galsim.c
clean:
	$(RM) $(SORT) $(OBJS)
