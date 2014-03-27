#FC=gfortran
FC = ifort

FFLAGS=-O3 
SOURCE=yukawa_beltrami.f
TARGET=yukawa_beltrami

OBJS=\
auxil.o\
conic.o\
cons.o\
dasum.o\
daxpy.o\
dcfft.o\
dcopy.o\
ddot.o\
dgeco.o\
dgefa.o\
dgesl.o\
dgmres.o\
dhels.o\
dheqr.o\
dnrm2.o\
dorth.o\
dpigmr.o\
drlcal.o\
dscal.o\
dxlcal.o\
hyp_2F1.o\
idamax.o\
isdgmr.o\
matplot.o\
prini.o\
quad2.o

yukawa_beltrami: ${OBJS} yukawa_beltrami.f
	${FC} ${FFLAGS} yukawa_beltrami.f -o ${TARGET} \
	  ${OBJS}

auxil.o: auxil.f90
	${FC} ${FFLAGS} -c auxil.f90

conic.o: conic.f90
	${FC} ${FFLAGS} -c conic.f90

cons.o: cons.f90
	${FC} ${FFLAGS} -c cons.f90

dasum.o: dasum.f
	${FC} ${FFLAGS} -c dasum.f

daxpy.o: daxpy.f
	${FC} ${FFLAGS} -c daxpy.f

dcfft.o: dcfft.f
	${FC} ${FFLAGS} -c dcfft.f

dcopy.o: dcopy.f
	${FC} ${FFLAGS} -c dcopy.f

ddot.o: ddot.f
	${FC} ${FFLAGS} -c ddot.f

dgeco.o: dgeco.f
	${FC} ${FFLAGS} -c dgeco.f

dgefa.o: dgefa.f
	${FC} ${FFLAGS} -c dgefa.f

dgesl.o: dgesl.f
	${FC} ${FFLAGS} -c dgesl.f

dgmres.o: dgmres.f
	${FC} ${FFLAGS} -c dgmres.f

dhels.o: dhels.f
	${FC} ${FFLAGS} -c dhels.f

dheqr.o: dheqr.f
	${FC} ${FFLAGS} -c dheqr.f

dnrm2.o: dnrm2.f
	${FC} ${FFLAGS} -c dnrm2.f

dorth.o: dorth.f
	${FC} ${FFLAGS} -c dorth.f

dpigmr.o: dpigmr.f
	${FC} ${FFLAGS} -c dpigmr.f

drlcal.o: drlcal.f
	${FC} ${FFLAGS} -c drlcal.f

dscal.o: dscal.f
	${FC} ${FFLAGS} -c dscal.f

dxlcal.o: dxlcal.f
	${FC} ${FFLAGS} -c dxlcal.f

hyp_2F1.o: hyp_2F1.f90
	${FC} ${FFLAGS} -c hyp_2F1.f90

idamax.o: idamax.f
	${FC} ${FFLAGS} -c idamax.f

isdgmr.o: isdgmr.f
	${FC} ${FFLAGS} -c isdgmr.f

matplot.o: matplot.f
	${FC} ${FFLAGS} -c matplot.f

prini.o: prini.f
	${FC} ${FFLAGS} -c prini.f

quad2.o: quad2.f
	${FC} ${FFLAGS} -c quad2.f


clean:
	rm *.o
	rm yukawa_beltrami
