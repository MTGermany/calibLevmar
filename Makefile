#
# Unix/Linux GCC Makefile for Levenberg - Marquardt minimization
# Under windows, use Makefile.vc for MSVC
#

CC=gcc -fPIC
CPP=g++ -fPIC
CONFIGFLAGS=#-ULINSOLVERS_RETAIN_MEMORY
#ARCHFLAGS=-march=pentium4 # YOU MIGHT WANT TO UNCOMMENT THIS FOR P4
CFLAGS=$(CONFIGFLAGS) $(ARCHFLAGS) -O3 -funroll-loops -Wall #-ffast-math #-pg
LAPACKLIBS_PATH=/usr/local/lib # WHEN USING LAPACK, CHANGE THIS TO WHERE YOUR COMPILED LIBS ARE!
LDFLAGS=-L$(LAPACKLIBS_PATH) -L.
LIB_LEVMAR=levmar-2.5
LIBOBJS=${LIB_LEVMAR}/lm.o ${LIB_LEVMAR}/Axb.o ${LIB_LEVMAR}/misc.o ${LIB_LEVMAR}/lmlec.o ${LIB_LEVMAR}/lmbc.o ${LIB_LEVMAR}/lmblec.o ${LIB_LEVMAR}/lmbleic.o
LIBSRCS=${LIB_LEVMAR}/lm.c ${LIB_LEVMAR}/Axb.c ${LIB_LEVMAR}/misc.c ${LIB_LEVMAR}/lmlec.c ${LIB_LEVMAR}/lmbc.c ${LIB_LEVMAR}/lmblec.c ${LIB_LEVMAR}/lmbleic.c
DEMOBJS=lmdemo.o
DEMOSRCS=lmdemo.c
DEMOBJS2=expfit.o
DEMOSRCS2=expfit.cpp
CALIBTRAJOBJ=calibTraj.o
CALIBTRAJSRC=calibTraj.cpp
EFCOBJ=makeEFCfromFCdata.o
EFCSRC=makeEFCfromFCdata.cpp

AR=ar
RANLIB=ranlib

#LAPACKLIBS=-llapack -lblas -lf2c # comment this line if you are not using LAPACK.
                                 # On systems with a FORTRAN (not f2c'ed) version of LAPACK, -lf2c is
                                 # not necessary; on others, -lf2c is equivalent to -lF77 -lI77

#LAPACKLIBS=-L/usr/local/atlas/lib -llapack -lcblas -lf77blas -latlas -lf2c # This works with the ATLAS updated lapack and Linux_P4SSE2
                                                                            # from http://www.netlib.org/atlas/archives/linux/

#LAPACKLIBS=-llapack -lgoto2 -lpthread -lf2c # This works with GotoBLAS
                                             # from http://www.tacc.utexas.edu/research-development/tacc-projects/

#LAPACKLIBS=-L/opt/intel/mkl/8.0.1/lib/32/ -lmkl_lapack -lmkl_ia32 -lguide -lf2c # This works with MKL 8.0.1 from
                                            # http://www.intel.com/cd/software/products/asmo-na/eng/perflib/mkl/index.htm

LIBS=$(LAPACKLIBS)

#own
MY_LIBOBJECTS=InOut.o Statistics.o Math.o RandomUtils.o general.o
MY_LIBDIR=my-lib
BINDIR=~/bin
###########################################################
# suffix regel: mache aus *.cpp ein *.o, und zwar fuer eingabedatei $<
.cpp.o:
	${CC} -I ${MY_LIBDIR}  -c $<

############################################################

general.o: ${MY_LIBDIR}/general.cpp
	${CC} -c ${MY_LIBDIR}/general.cpp -o general.o

InOut.o: ${MY_LIBDIR}/InOut.cpp
	${CC} -c ${MY_LIBDIR}/InOut.cpp -o InOut.o

Statistics.o: ${MY_LIBDIR}/Statistics.cpp
	${CC} -c ${MY_LIBDIR}/Statistics.cpp -o Statistics.o

Math.o: ${MY_LIBDIR}/Math.cpp
	${CC} -c ${MY_LIBDIR}/Math.cpp -o Math.o

RandomUtils.o: ${MY_LIBDIR}/RandomUtils.cpp
	${CC} -c ${MY_LIBDIR}/RandomUtils.cpp -o RandomUtils.o
#end own

all: liblevmar.a lmdemo expfit calibTraj calibTrajTest calibDiscrChoice

liblevmar.a: $(LIBOBJS)
	$(AR) crv liblevmar.a $(LIBOBJS)
	$(RANLIB) liblevmar.a

lmdemo: $(DEMOBJS) liblevmar.a
	$(CC) $(LDFLAGS) $(DEMOBJS) -o lmdemo -llevmar $(LIBS) -lm
expfit: $(DEMOBJS2) liblevmar.a
	$(CPP) $(LDFLAGS) $(DEMOBJS2) -o expfit -llevmar $(LIBS) -lm
#own
calibTraj: $(CALIBTRAJOBJ) $(MY_LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) $(CALIBTRAJOBJ) $(MY_LIBOBJECTS) -o $(BINDIR)/calibTraj -llevmar $(LIBS) -lm
calibTrajTest: calibTrajTest.o $(MY_LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) calibTrajTest.o $(MY_LIBOBJECTS) -o $(BINDIR)/calibTrajTest -llevmar $(LIBS) -lm

makeEFCfromFCdata: $(EFCOBJ) $(MY_LIBOBJECTS)
	$(CPP) $(LDFLAGS) $(EFCOBJ) $(MY_LIBOBJECTS) -o $(BINDIR)/makeEFCfromFCdata -lm

determineDrivingRegimes: determineDrivingRegimes.o $(MY_LIBOBJECTS)
	$(CPP) $(LDFLAGS) determineDrivingRegimes.o $(MY_LIBOBJECTS) -o $(BINDIR)/determineDrivingRegimes -lm

calibDiscrChoice: calibDiscrChoice.o $(MY_LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) calibDiscrChoice.o  $(MY_LIBOBJECTS) -o $(BINDIR)/calibDiscrChoice -llevmar $(LIBS) -lm

#end own

lm.o: ${LIB_LEVMAR}/lm.c ${LIB_LEVMAR}/lm_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h ${LIB_LEVMAR}/compiler.h
Axb.o: ${LIB_LEVMAR}/Axb.c ${LIB_LEVMAR}/Axb_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h
misc.o:${LIB_LEVMAR}/misc.c ${LIB_LEVMAR}/misc_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h
lmlec.o: ${LIB_LEVMAR}/lmlec.c ${LIB_LEVMAR}/lmlec_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h
lmbc.o: ${LIB_LEVMAR}/lmbc.c ${LIB_LEVMAR}/lmbc_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h ${LIB_LEVMAR}/compiler.h
lmblec.o: ${LIB_LEVMAR}/lmblec.c ${LIB_LEVMAR}/lmblec_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h
lmbleic.o: ${LIB_LEVMAR}/lmbleic.c ${LIB_LEVMAR}/lmbleic_core.c ${LIB_LEVMAR}/levmar.h ${LIB_LEVMAR}/misc.h

lmdemo.o: ${LIB_LEVMAR}/levmar.h
expfit.o: ${LIB_LEVMAR}/levmar.h
calibTraj.o: ${LIB_LEVMAR}/levmar.h
calibTrajTest.o: ${LIB_LEVMAR}/levmar.h

clean:
	@rm -f $(LIBOBJS) $(DEMOBJS) $(MY_LIBOBJECTS) *.o expfit

cleanall: clean
	@rm -f lmdemo
	@rm -f expfit
	@rm -f calibTraj
	@rm -f calibTrajTest
	@rm -f liblevmar.a

