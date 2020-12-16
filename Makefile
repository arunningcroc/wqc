all: scquantumoptics

WARNINGS = -Wall
DEBUG =-std=c99 -ggdb -fno-omit-frame-pointer
OPTIMIZE = -O2 -lm

scquantumoptics: Makefile scquantumoptics.c operators.c pool.c ode.c
	$(CC) -o $@ $(WARNINGS) $(DEBUG) $(OPTIMIZE) scquantumoptics.c \
	operators.c pool.c ode.c

debug: Makefile scquantumoptics.c operators.c pool.c ode.c
	$(CC) -o $@ $(WARNINGS) $(DEBUG) $(OPTIMIZE) -pg scquantumoptics.c operators.c pool.c ode.c
clean:
	rm -f scquantumoptics

# Builder will call this to install the application before running.
install:
	echo "Installing is not supported"

# Builder uses this target to run your application.
run:
	./scquantumoptics
test:
	./scquantumoptics
