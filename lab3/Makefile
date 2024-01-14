#default build: gameoflife
#all: gameoflife
all: gameoflife_mpi

measure_time_mpi: measure_time_mpi.c
	mpicc -O3 -std=gnu99 measure_time_mpi.c -o measure_time_mpi

gameoflife_mpi: gameoflife_mpi.c
	mpicc -O3 -std=gnu99 gameoflife_mpi.c -o gameoflife_mpi

heat_equation_mpi: heat_equation_mpi.c
	mpicc -O3 -lpng -std=gnu99 heat_equation_mpi.c -o heat_equation_mpi

clean:
	rm -r gol
	rm -r heq
	rm measure_time_mpi
	rm gameoflife_mpi
	rm heat_equation_mpi
