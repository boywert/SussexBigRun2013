SYSTEM =	"Sisu+Intel"

CC=mpicc
FC=mpif90

ifeq ($(SYSTEM), "Sisu+Intel")
	CC         	=	cc -std=c99 -Wall -W -fopenmp
	FC         	=	ftn

endif

all: chunk

chunk: 
	$(FC) -cpp -O3 chunk_cubep3m_v1.5.f90 -o chunk 	

clean:
	rm -f chunk