#INCLUDES = syspara.h

bin_PROGRAMS = ORd_without 
#ORd_CFLAGS = -O2 -Wall -axsse4.2 -mkl -parallel -par-report2 -Winline -inline-level=2 -ip -vec-report2 -opt-mem-layout-trans=3 -pg
ORd_without_CFLAGS = -O2 -axsse4.2 -mkl -parallel -Winline -inline-level=2 -ip -opt-mem-layout-trans=3 -pg

ORd_without_INCLUDE = -I./include
ORd_without_LDADD = -L/opt/intel/composer_xe_2015.1.133/mkl/lib/intel64 -liomp5 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core

ORd_without_SOURCES = ./src/ord.c ./src/val_constants.c ./src/mem.c\
	./src/eular.c ./src/input_para.c ./src/system.c ./src/currents.c\
	./src/exptable.c ./src/syspara.h

