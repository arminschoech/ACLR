#include "snps.h"


Snps::Snps(const char* impSnpFilePath, const char* regSnpFilePath, int inBatchSize)
{
    batchSize = inBatchSize;
    nextImpSnpNum = 0;

    // read in imputation SNP indices
    FILE* impSnpFile = fopen(impSnpFilePath, "rb");
    if (impSnpFile == NULL) {
        std::cerr << "ERROR: 'imp_snp.bin' file can't be opened." << std::endl;
        exit(1);
    }

    fseek(impSnpFile, 0, SEEK_END);
    int impSnpNum = ftell(impSnpFile) / 4;
    rewind(impSnpFile);
    anyImpSnpsLeft = (impSnpNum > 0);

    impSnps = gsl_vector_int_alloc(impSnpNum);
    int currentIndex;
    for(int i = 0; i < impSnpNum; i++) {
        fread(&currentIndex, 4, 1, impSnpFile);
        gsl_vector_int_set(impSnps, i, currentIndex);
        // std::cout << currentIndex << std::endl;
    }

    FILE* regSnpFile = fopen(regSnpFilePath, "rb");
    if (regSnpFile == NULL) {
        std::cerr << "ERROR: 'reg_snp.bin' file can't be opened." << std::endl;
        exit(1);
    }
    fseek(regSnpFile, 0, SEEK_END);
    int regSnpNum = ftell(regSnpFile) / 4;
    rewind(regSnpFile);

    regSnps = gsl_vector_int_alloc(regSnpNum);
    for(int i = 0; i < regSnpNum; i++) {
        fread(&currentIndex, 4, 1, regSnpFile);
        gsl_vector_int_set(regSnps, i, currentIndex);
        // std::cout << currentIndex << std::endl;
    }
    
    // std::cout << "There are " << impSnpNum << " imp SNPs and " << regSnpNum << " reg SNPs." << std::endl;

}

gsl_vector_int* Snps::getNextImpSnpBatch() {
    int fromSnp = nextImpSnpNum; //included
    int toSnp; //not included

    if(fromSnp + batchSize < (int) impSnps->size) {
        toSnp = fromSnp + batchSize;
    } else {
        toSnp = impSnps->size;
        anyImpSnpsLeft = false;
    }
    snpsInBatch = toSnp - fromSnp;
    
    gsl_vector_int* impSnpBatch = gsl_vector_int_alloc(snpsInBatch);
    for(int i = 0; i < snpsInBatch; i++) {
        gsl_vector_int_set(impSnpBatch, i, gsl_vector_int_get(impSnps, fromSnp + i));
    }

    nextImpSnpNum = toSnp;
    return impSnpBatch;
}

