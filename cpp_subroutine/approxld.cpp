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
    int batchSize = 100;
    const char* impSnpBinFilePath = "imp_snps.bin";
    const char* outFilePath = "approx_ld_mat.bin";
    const char* bgenPathFilePath = "bgen_path.txt";
    const char* indiFilePath = "indi_bgen_pos_subset.bin";

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

	Bgen bgenLink(bgenFilePath, indiFilePath);
    int indiNum = bgenLink.getN();
    Snps snpIndices(impSnpBinFilePath, impSnpBinFilePath, batchSize);
    int snpNum = snpIndices.getImpSnps()->size;

    Geno allGenotypes(indiNum, snpNum); 
    Geno batchGenotypes(indiNum, batchSize);
    gsl_matrix_float* ldMatrix = gsl_matrix_float_alloc(batchSize, snpNum);
    FILE* outFile = fopen(outFilePath, "wb");
    
	// cout << endl << "Reading in regression SNPs ..." << endl; 
    allGenotypes.readInGenMat(&bgenLink, snpIndices.getImpSnps());

    while(snpIndices.anyImpSnpsLeft) {
        // cout << "Processing imputed SNPs starting at SNP number " << snpIndices.getNextImpSnpNum() << " ..." << endl;
        batchGenotypes.readInGenMat(&bgenLink, snpIndices.getNextImpSnpBatch());

        // matrix multiplication
        gsl_blas_sgemm(CblasTrans, CblasNoTrans, 1.0 / indiNum, batchGenotypes.getGenMatAddress(), allGenotypes.getGenMatAddress(), 0.0, ldMatrix);

        // write to file
        gsl_matrix_float_fwrite(outFile, ldMatrix);
    }

    fclose(outFile);
	return 0;
}


