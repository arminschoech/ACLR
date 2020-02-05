#ifndef SNPS_H_
#define SNPS_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <gsl/gsl_vector.h>


class Snps 
{

public:
    Snps(const char* impSnpFilePath, const char* regSnpFilePath, int inBatchSize);
    gsl_vector_int* getRegSnps();
    gsl_vector_int* getImpSnps();
    gsl_vector_int* getNextImpSnpBatch();
    bool anyImpSnpsLeft;
    int snpsInBatch;
    int getNextImpSnpNum();

private:
    int batchSize;
    int nextImpSnpNum;
	gsl_vector_int* impSnps;
	gsl_vector_int* regSnps;
    
};


inline gsl_vector_int* Snps::getImpSnps()
{
    return impSnps;
}

inline gsl_vector_int* Snps::getRegSnps()
{
    return regSnps;
}

inline int Snps::getNextImpSnpNum() {
    return nextImpSnpNum;
}

#endif
