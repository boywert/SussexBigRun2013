simuPdens:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuCroCo:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuPk:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuXi:			DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuVstat:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
simuSigmaR:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
MergerTree:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindHalo:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindHaloPairs:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}
ahfFindCrossPairs:		DEFINEFLAGS	=	-DVERBOSE -DHEADER ${MASTER_DEFINEFLAGS}

CFLAGS  =       ${OPTIMIZE} ${CCFLAGS} ${DEFINEFLAGS}
LFLAGS  =       ${LNFLAGS}

export DEFINEFLAGS

LIBS = 	../src/libio_serial/libio_serial.a \
	../src/libio/libio.a \
	../src/libahf/libahf.a \
	../src/libgravity/libgravity.a \
	../src/libutility/libutility.a \
	../src/libamr_serial/libamr_serial.a


# HOW TO MAKE THE ACTUAL PROGRAMMES
simuPdens:	simuPdens.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuPdens simuPdens.c ${LIBS} ${LIBS} -lm

simuCroCo:	simuCroCo.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuCroCo simuCroCo.c ${LIBS} ${LIBS} -lm

simuPk:	simuPk.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuPk simuPk.c ${LIBS} ${LIBS} -lm

simuXi:	simuXi.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuXi simuXi.c ${LIBS} ${LIBS} -lm

simuVstat:	simuVstat.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuVstat simuVstat.c ${LIBS} ${LIBS} -lm

simuSigmaR:	simuSigmaR.c ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o simuSigmaR simuSigmaR.c ${LIBS} ${LIBS} -lm

MergerTree:	MergerTree.c  ${LIBS}
	$(CC) $(CFLAGS) $(LFLAGS) -o MergerTree  ${LIBS} MergerTree.c  ${LIBS} -lm

ahfFindHalo:	ahfFindHalo.c
	$(CC) $(CFLAGS) $(LFLAGS) -o ahfFindHalo ahfFindHalo.c -lm

ahfFindHaloPairs:	ahfFindHaloPairs.c
	$(CC) $(CFLAGS) $(LFLAGS) -o ahfFindHaloPairs ahfFindHaloPairs.c -lm

ahfFindCrossPairs:	ahfFindCrossPairs.c
	$(CC) $(CFLAGS) $(LFLAGS) -o ahfFindCrossPairs ahfFindCrossPairs.c -lm

ahfHaloHistory:		ahfHaloHistory.f
	$(FC) -O2 -o ahfHaloHistory ahfHaloHistory.f

ahfSubCheck:		ahfSubCheck.f
	$(FC) -O2 -o ahfSubCheck ahfSubCheck.f

ahfXi:		ahfXi.f
	$(FC) -O2 -o ahfXi ahfXi.f


# HOW TO MAKE THE LIBRARIES
../src/libio_serial/libio_serial.a:
	cd ../src/libio_serial ; ${MAKE} libio_serial.a

../src/libio/libio.a:
	cd ../src/libio ; ${MAKE} libio.a

../src/libahf/libahf.a:
	cd ../src/libahf ; ${MAKE} libahf.a

../src/libgravity/libgravity.a:
	cd ../src/libgravity ; ${MAKE} libgravity.a

../src/libamr_serial/libamr_serial.a:
	cd ../src/libamr_serial ; ${MAKE} libamr_serial.a

../src/libutility/libutility.a:
	cd ../src/libutility ; ${MAKE} libutility.a



clean:
	rm -f  *.o

veryclean:
	rm -f  *.o *~ *~.*
