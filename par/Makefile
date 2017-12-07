include Makefile.in

.PHONY: default
default: cgp3d.x

# === Code

cgp3d.x: cgp3d.o pcg.o params.o timing.o
	$(CC) -o $@ $^ $(LFLAGS) $(LIBS) $(OMPFLAGS)

.c.o:
	$(CC) $(CFLAGS) $(OMPFLAGS) $(INCLUDE) -c $<

.PHONY: clean
clean:
	rm -f *.o *~ cgp3d.x
