#include "geno.h"


Geno::Geno(int inN, int inM) 
{
	n = inN;
	mBuffer = inM;
	mNextToFill = 0;

	pGenotypes = gsl_matrix_float_calloc(n, mBuffer);
    // TODO: use alloc instead for speed?
	pTempVec = gsl_vector_float_alloc(n);

}


void Geno::readInGenVec(Bgen* bgenLink)
{
	if(mNextToFill >= mBuffer) {
      std::cerr << "ERROR: Cannot write to genotype matrix. Genotype matrix is full. " << std::endl;
      exit(1);
    }

	float freq = bgenLink->getFrequency();
    float maf = freq;
	if(freq > 0.5) maf = 1.0 - maf;

	// calculate SNP frequency p and normalization constant	
	if(maf <= 0.0) 
	{
		gsl_vector_float_set_all(pTempVec, 0.0);
		gsl_matrix_float_set_col(pGenotypes, mNextToFill, pTempVec);
		// if all genotypes are the same (AA or BB), fill in zeros in the 
	} 
	else 
	{
		if(freq <= 0.5) {
			gsl_blas_scopy(bgenLink->genVec, pTempVec);
			gsl_vector_float_add_constant(pTempVec, -2.0 * maf);
			// gsl_blas_sscal(normalizationConst, pTempVec);
		} else {
			gsl_blas_scopy(bgenLink->genVec, pTempVec);
			gsl_blas_sscal(-1.0, pTempVec);
			gsl_vector_float_add_constant(pTempVec, 2.0);
			gsl_vector_float_add_constant(pTempVec, -2.0 * maf);
			// gsl_blas_sscal(normalizationConst, pTempVec);
		}

		// write normalized vector to column #mNextToFill in the genotype matrix
		gsl_matrix_float_set_col(pGenotypes, mNextToFill, pTempVec);
	}
	mNextToFill++;
}


void Geno::readInGenMat(Bgen* bgenLink, gsl_vector_int* snpIndices)
{
    int snpIndicesNum = (int) snpIndices->size;
    if(snpIndicesNum > mBuffer) {
        std::cerr << "ERROR: SNP index vector to be read in is larger than buffer." << std::endl;
      exit(1);
    }

    resetGenotypeMatrix();
    bgenLink->rewindBgenFile();
    int lastSnp;
    int nextSnp;

    lastSnp = -1;
    for(int i = 0; i < snpIndicesNum; i++) {
        nextSnp = gsl_vector_int_get(snpIndices, i);
        bgenLink->skipSnps(nextSnp - lastSnp - 1);
        bgenLink->calcNextSnpDosage();
        readInGenVec(bgenLink);
        lastSnp = nextSnp;
    }
}


void Geno::printGenotypes()
{
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < mBuffer; j++) {
			std::cout << gsl_matrix_float_get(pGenotypes, i, j) << " ";
		}
		std::cout << std::endl;
	}
}


void Geno::printGenotypes(int maxN, int maxM)
{
	for(int i = 0; i < maxN; i++) {
		for(int j = 0; j < maxM; j++) {
			std::cout << gsl_matrix_float_get(pGenotypes, i, j) << " ";
		}
		std::cout << std::endl;
	}
}



