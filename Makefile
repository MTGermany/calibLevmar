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
LIBOBJS=lm.o Axb.o misc.o lmlec.o lmbc.o lmblec.o lmbleic.o
LIBSRCS=lm.c Axb.c misc.c lmlec.c lmbc.c lmblec.c lmbleic.c
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
LIBOBJECTS=InOut.o Statistics.o Math.o RandomUtils.o general.o
LIBDIR=~/versionedProjects/lib/trunk
BINDIR=~/bin
###########################################################
# suffix regel: mache aus *.cpp ein *.o, und zwar fuer eingabedatei $<
.cpp.o:
	${CC} -I ${LIBDIR}  -c $<

############################################################

general.o: ${LIBDIR}/general.cpp
	${CC} -c ${LIBDIR}/general.cpp -o general.o

InOut.o: ${LIBDIR}/InOut.cpp
	${CC} -c ${LIBDIR}/InOut.cpp -o InOut.o

Statistics.o: ${LIBDIR}/Statistics.cpp
	${CC} -c ${LIBDIR}/Statistics.cpp -o Statistics.o

Math.o: ${LIBDIR}/Math.cpp
	${CC} -c ${LIBDIR}/Math.cpp -o Math.o

RandomUtils.o: ${LIBDIR}/RandomUtils.cpp
	${CC} -c ${LIBDIR}/RandomUtils.cpp -o RandomUtils.o
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
calibTraj: $(CALIBTRAJOBJ) $(LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) $(CALIBTRAJOBJ) $(LIBOBJECTS) -o $(BINDIR)/calibTraj -llevmar $(LIBS) -lm
calibTrajTest: calibTrajTest.o $(LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) calibTrajTest.o $(LIBOBJECTS) -o $(BINDIR)/calibTrajTest -llevmar $(LIBS) -lm

makeEFCfromFCdata: $(EFCOBJ) $(LIBOBJECTS)
	$(CPP) $(LDFLAGS) $(EFCOBJ) $(LIBOBJECTS) -o $(BINDIR)/makeEFCfromFCdata -lm

determineDrivingRegimes: determineDrivingRegimes.o $(LIBOBJECTS)
	$(CPP) $(LDFLAGS) determineDrivingRegimes.o $(LIBOBJECTS) -o $(BINDIR)/determineDrivingRegimes -lm

calibDiscrChoice: calibDiscrChoice.o $(LIBOBJECTS) liblevmar.a
	$(CPP) $(LDFLAGS) calibDiscrChoice.o  $(LIBOBJECTS) -o $(BINDIR)/calibDiscrChoice -llevmar $(LIBS) -lm

#end own

lm.o: lm.c lm_core.c levmar.h misc.h compiler.h
Axb.o: Axb.c Axb_core.c levmar.h misc.h
misc.o: misc.c misc_core.c levmar.h misc.h
lmlec.o: lmlec.c lmlec_core.c levmar.h misc.h
lmbc.o: lmbc.c lmbc_core.c levmar.h misc.h compiler.h
lmblec.o: lmblec.c lmblec_core.c levmar.h misc.h
lmbleic.o: lmbleic.c lmbleic_core.c levmar.h misc.h

lmdemo.o: levmar.h
expfit.o: levmar.h
calibTraj.o: levmar.h
calibTrajTest.o: levmar.h

clean:
	@rm -f $(LIBOBJS) $(DEMOBJS)

cleanall: clean
	@rm -f lmdemo
	@rm -f expfit
	@rm -f calibTraj
	@rm -f calibTrajTest
	@rm -f liblevmar.a

