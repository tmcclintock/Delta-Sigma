OBJS = src/xi_mm/xi_mm.o src/tinker_bias/tinker_bias.o src/tinker_bias/tinker_bias_at_M.o src/xi_nfw/xi_nfw.o src/xi_2halo/xi_2halo.o src/xi_hm/xi_hm.o src/sigma/sigma.o src/sigma/sigma_at_r.o src/delta_sigma/delta_sigma.o src/delta_sigma/delta_sigma_at_r.o src/sigma_mis/sigma_mis_at_r.o src/sigma_mis/sigma_mis.o src/delta_sigma_mis/delta_sigma_mis_at_r.o src/delta_sigma_mis/delta_sigma_mis.o src/miscentered_sigma/miscentered_sigma.o src/miscentered_sigma/miscentered_sigma_at_r.o src/miscentered_delta_sigma/miscentered_delta_sigma.o src/miscentered_delta_sigma/miscentered_delta_sigma_at_r.o src/ave_delta_sigma/ave_delta_sigma.o src/ave_delta_sigma/ave_delta_sigma_in_bin.o src/ave_miscentered_delta_sigma/ave_miscentered_delta_sigma.o src/ave_miscentered_delta_sigma/ave_miscentered_delta_sigma_in_bin.o src/mc_relation/mc_relation.o src/wrapper/wrapper.o

CC = gcc
EXEC =
CEXEC = ./main.exe
PYEXEC = ./src/wrapper/Delta_Sigma.so

ifdef ALONE
ifeq ($(ALONE),yes)
EXEC += $(CEXEC)
CFLAGS = -g
OFLAGS = 
OBJS += main.o
endif
else
EXEC += $(PYEXEC)
CFLAGS = -fPIC
OFLAGS = -shared 
endif

INCL = -I${GSLI} -fopenmp -O2
LIBS = -lgsl -lgslcblas -L${GSLL} -lm -fopenmp  -O2
DFLAGS =

all : $(EXEC)
	@echo "Compilation complete."

timing : DFLAGS += -DTIMING
timing : $(EXEC)

%.o: %.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCL) -c $< -o $@

$(EXEC) : $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

.PHONY : clean all

clean:
	@rm -f $(OBJS) main.o $(CEXEC) $(PYEXEC)
	@echo "Cleanup complete."