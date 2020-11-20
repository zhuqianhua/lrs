GXX = g++
LDP = -L/hwfswh2/BC_PUB/Software/07.User-defined/01.Development/c_plus_plus_lib/htslib_v1.10.2/lib
LDL = -lpthread -lhts
INC = /hwfswh2/BC_PUB/Software/07.User-defined/01.Development/c_plus_plus_lib/htslib_v1.10.2/include
OBJ = lrs_opts.o lrs_subs.o lrs_core.o lrs.o
VPATH = src

lrs : ${OBJ}
	${GXX} -o lrs src/sw/ssw_cpp.cpp src/sw/ssw.c ${OBJ} ${LDL} ${LDP}
.PHONY : clean
clean :
	rm -f ${OBJ}

lrs_opts.o : lrs_opts.h
lrs_subs.o : lrs_subs.h
lrs_core.o : lrs_core.h
lrs.o : lrs.cpp lrs.h
	${GXX} -c $< -I${INC}
