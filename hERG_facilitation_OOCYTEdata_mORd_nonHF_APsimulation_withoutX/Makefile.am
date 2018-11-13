#INCLUDES = syspara.h

bin_PROGRAMS = ORd 
ORd_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg

ORd_INCLUDE = -I./include
ORd_LDADD = -L/opt/intel/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core ./lib/xhplot.o

ORd_SOURCES = ./src/ord.c ./src/val_constants.c ./src/mem.c\
	./src/win_out.c ./src/runge.c ./src/eular.c ./src/input_para.c ./src/system.c ./src/currents.c\
	./src/eventloop.c ./src/exptable.c ./src/data_out.c ./src/syspara.h

