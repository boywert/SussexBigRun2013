SYSTEM =	"Sisu+Intel"

CC=mpicc
FC=mpif90

ifeq ($(SYSTEM), "Sisu+Intel")
	CC         	=	cc -std=c99 -Wall -W -fopenmp
	FC         	=	ftn -fpp -O3 -DBINARY

endif

all: chunk

chunk: 
	$(FC) chunk_cubep3m_v1.5.f90 -o chunk 	

clean:
	rm -f chunk