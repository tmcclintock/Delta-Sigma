EXEC = main.exe
OBJS = main.o src/xi_mm/xi_mm_at_r.o src/xi_mm/xi_mm.o src/tinker_bias/tinker_bias.o src/tinker_bias/tinker_bias_at_M.o src/xi_nfw/xi_nfw.o
CC = gcc
INCL = -I/home/tom/code/gsl/include -fopenmp
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm -fopenmp
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
