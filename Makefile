# Reset the terminal
# Reset=`tput sgr0`

RED='\e[1;31m'
GREEN='\e[1;32m'
YELLOW='\e[1;33m'
WHITE='\e[0m'

OK_STRING=$(GREEN)[OK]$(WHITE)
ERROR_STRING=$(RED)[ERRORS]$(WHITE)
WARN_STRING=$(YELLOW)[WARNINGS]$(WHITE)

CC = g++
CFLAGS = -Wall
LDFLAGS = 

LSDIR = $(shell ls)

OBJS = particles.o 1DPlasmaModel.o

Plasma: $(OBJS)
	@$(CC) -o Run $^ $(LDFLAGS)

%.o: %.cpp 
	@$(CC) $(CFLAGS) -c $<

clean: 
	rm -f *.o
	rm -f *.exe
	rm -f *.path
	rm -f Run