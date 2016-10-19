OBJS = main.o src/xi_mm/xi_mm.o src/tinker_bias/tinker_bias.o src/tinker_bias/tinker_bias_at_M.o src/xi_nfw/xi_nfw.o src/xi_2halo/xi_2halo.o src/xi_hm/xi_hm.o src/sigma_r/sigma_r.o src/sigma_r/sigma_r_at_r.o src/delta_sigma/delta_sigma.o src/delta_sigma/delta_sigma_at_r.o src/miscentered_sigma_r/miscentered_sigma_r.o src/miscentered_sigma_r/miscentered_sigma_r_at_r.o src/miscentered_delta_sigma/miscentered_delta_sigma.o src/miscentered_delta_sigma/miscentered_delta_sigma_at_r.o src/ave_delta_sigma/ave_delta_sigma.o src/ave_delta_sigma/ave_delta_sigma_in_bin.o src/miscentered_ave_delta_sigma/miscentered_ave_delta_sigma.o src/miscentered_ave_delta_sigma/miscentered_ave_delta_sigma_in_bin.o src/mc_relation/mc_relation.o src/wrapper/wrapper.o

CC = gcc
ifdef SHARED
ifeq ($(SHARED),yes)
$(info Building shared library)
EXEC = Delta_Sigma.so
CFLAGS = -fPIC
OFLAGS = -shared 
#-W1,-soname=$(EXEC)
endif
else
$(info Building executable)
EXEC = main.exe
CFLAGS = 
OFLAGS = 
endif

INCL = -I/${GSLI} -fopenmp -O3
LIBS = -lgsl -lgslcblas -L/${GSLL} -lm -fopenmp -O3
.SUFFIXES : .c .o
%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

#$(OBJS): $(INCL)

.PHONY : clean

clean:
	rm -f $(OBJS) main.exe Delta_Sigma.so
	rm -f *~
