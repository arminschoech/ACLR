#include "bgen.h"
#include "geno.h"
#include "snps.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <string>

using std::cout;
using std::endl;
using std::cerr;


int main(int argc, char* argv[])
{
    int impBatchSize = 100;
    const char* impSnpBinFilePath = "imp_snps.bin";
    const char* regSnpBinFilePath = "reg_snps.bin";
    const char* outFilePath = "ld_mat.bin";
    const char* bgenPathFilePath = "bgen_path.txt";
    const char* indiPathFilePath = "indi_path.txt";

    // read in BGEN file path
    std::ifstream bgenPathFile(bgenPathFilePath);
    std::string bgenPathString;
    if(bgenPathFile.is_open()) {
        getline(bgenPathFile, bgenPathString);
    } else {
        cerr << "ERROR: unable to open '" << bgenPathFilePath << "'." << endl;
        exit(1);
    }
    const char* bgenFilePath = bgenPathString.c_str();

    // read in individual file path
    std::ifstream indiPathFile(indiPathFilePath);
    std::string indiPathString;
    if(indiPathFile.is_open()) {
        getline(indiPathFile, indiPathString);
    } else {
        cerr << "ERROR: unable to open '" << indiPathFilePath << "'." << endl;
        exit(1);
    }
    const char* indiFilePath = indiPathString.c_str();

	Bgen bgenLink(bgenFilePath, indiFilePath);
    int indiNum = bgenLink.getN();
    Snps snpIndices(impSnpBinFilePath, regSnpBinFilePath, impBatchSize);
    int regSnpNum = snpIndices.getRegSnps()->size;

    Geno regGenotypes(indiNum, regSnpNum); 
    Geno impGenotypes(indiNum, impBatchSize);
    gsl_matrix_float* ldMatrix = gsl_matrix_float_alloc(impBatchSize, regSnpNum);
    FILE* outFile = fopen(outFilePath, "wb");
    
	// cout << endl << "Reading in regression SNPs ..." << endl; 
    regGenotypes.readInGenMat(&bgenLink, snpIndices.getRegSnps());

    while(snpIndices.anyImpSnpsLeft) {
        // cout << "Processing imputed SNPs starting at SNP number " << snpIndices.getNextImpSnpNum() << " ..." << endl;
        impGenotypes.readInGenMat(&bgenLink, snpIndices.getNextImpSnpBatch());

        // matrix multiplication
        gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1.0 / indiNum, impGenotypes.getGenMatAddress(), regGenotypes.getGenMatAddress(), 0.0, ldMatrix);

        // write to file
        gsl_matrix_float_fwrite(outFile, ldMatrix);
    }

    fclose(outFile);
	return 0;
}


