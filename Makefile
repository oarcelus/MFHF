FC     = ifort

default: hf.x 

hf.x: 
	$(FC) -mkl constant.f90 utility.f90 parameters.f90 density.f90 hamiltonian.f90 energy.f90 plot.f90 main.f90

clean: 
	rm -f  *.mod 

