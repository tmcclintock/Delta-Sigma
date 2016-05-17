EXEC = main.exe
OBJS = main.o src/xi_mm/xi_mm_at_r.o src/xi_mm/xi_mm.o src/tinker_bias/tinker_bias.o src/tinker_bias/tinker_bias_at_M.o src/xi_nfw/xi_nfw.o src/xi_2halo/xi_2halo.o src/xi_hm/xi_hm.o src/sigma_r/sigma_r.o src/sigma_r/sigma_r_at_r.o src/delta_sigma/delta_sigma.o src/delta_sigma/delta_sigma_at_r.o src/miscentered_sigma_r/miscentered_sigma_r.o src/miscentered_sigma_r/miscentered_sigma_r_at_r.o src/miscentered_delta_sigma/miscentered_delta_sigma.o src/miscentered_delta_sigma/miscentered_delta_sigma_at_r.o src/ave_delta_sigma/ave_delta_sigma.o src/ave_delta_sigma/ave_delta_sigma_in_bin.o src/wrapper/wrapper.o
CC = gcc
INCL = -I/home/tom/code/gsl/include -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/home/tom/code/gsl/lib -lm -fopenmp -O3
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
