EXEC = main.exe
OBJS = main.o src/xi_mm/xi_mm_at_r.o src/xi_mm/xi_mm.o
CC = gcc
INCL = -I/home/tom/code/gsl/include
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) $(EXEC)
	rm -f *~
