TARGETS = kc2bc qp2bc

CC = gcc
CFLAGS = -std=c99

default : $(TARGETS)

% : %.c
	$(CC) $(CFLAGS) -o $@ $<

clean :
	rm $(TARGETS) *.bc *output*

