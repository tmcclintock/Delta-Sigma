OBJS = src/xi_mm/xi_mm.o src/tinker_bias/tinker_bias.o src/tinker_bias/tinker_bias_at_M.o src/xi_nfw/xi_nfw.o src/xi_2halo/xi_2halo.o src/xi_hm/xi_hm.o src/sigma/sigma.o src/sigma/sigma_at_r.o src/delta_sigma/delta_sigma.o src/delta_sigma/delta_sigma_at_r.o src/ave_delta_sigma/ave_delta_sigma.o src/ave_delta_sigma/ave_delta_sigma_in_bin.o src/wrapper/wrapper.o

CC = gcc
EXEC =
CEXEC = ./main.exe
PYEXEC = ./src/wrapper/Delta_Sigma.so

ifdef ALONE
ifeq ($(ALONE),yes)
EXEC += $(CEXEC)
CFLAGS = -pg
OFLAGS = -pg
OBJS += main.o
endif
else
EXEC += $(PYEXEC)
CFLAGS = -fPIC
OFLAGS = -shared 
endif

INCL = -I${GSLI} -O2
LIBS = -lgsl -lgslcblas -L${GSLL} -lm -O2
DFLAGS =

all : $(EXEC)
	@echo "Compilation complete."

%.o: %.c
	$(CC) $(CFLAGS) $(DFLAGS) $(INCL) -c $^ -o $@

$(EXEC) : $(OBJS)
	$(CC) $(OFLAGS) $^ $(LIBS) -o $(EXEC)

.PHONY : clean all

clean:
	@rm -f $(OBJS) main.o $(CEXEC) $(PYEXEC)
	@echo "Cleanup complete."
