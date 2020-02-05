
CXX = g++
CXXFLAGS = -O3 -Wall -fopenmp -I/n/app/gsl/2.3/include

ld_mat_calc: main.o bgen.o geno.o snps.o
	$(CXX) $(CXXFLAGS) -L/n/app/gsl/2.3/lib -L/n/app/intel/2016/mkl/lib/intel64 -o ld_mat_calc main.o bgen.o geno.o snps.o -lgsl -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -lz -Wl,-rpath,/n/app/intel/2016/mkl/lib/intel64:/n/app/intel/2016/mkl/lib/intel64 -Wl,-rpath,/n/app/gsl/2.3/lib:/n/app/gsl/2.3/lib

main.o: main.cpp bgen.h geno.h snps.h
	$(CXX) $(CXXFLAGS) -c main.cpp

bgen.o: bgen.cpp bgen.h
	$(CXX) $(CXXFLAGS) -c bgen.cpp

geno.o: geno.cpp geno.h bgen.h
	$(CXX) $(CXXFLAGS) -c geno.cpp

snps.o: snps.cpp snps.h
	$(CXX) $(CXXFLAGS) -c snps.cpp


approxld: approxld.o bgen.o geno.o snps.o
	$(CXX) $(CXXFLAGS) -L/n/app/gsl/2.3/lib -L/n/app/intel/2016/mkl/lib/intel64 -o approxld approxld.o bgen.o geno.o snps.o -lgsl -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -lpthread -lm -lz -Wl,-rpath,/n/app/intel/2016/mkl/lib/intel64:/n/app/intel/2016/mkl/lib/intel64 -Wl,-rpath,/n/app/gsl/2.3/lib:/n/app/gsl/2.3/lib

approxld.o: approxld.cpp bgen.h geno.h snps.h
	$(CXX) $(CXXFLAGS) -c approxld.cpp
