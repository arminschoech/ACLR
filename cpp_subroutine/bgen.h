#ifndef BGEN_H_
#define BGEN_H_


#include <iostream>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <cstdlib>
#include <string>
#include <gsl/gsl_vector.h>
#include <zlib.h>

using std::cout;
using std::endl;
using std::cerr;


class Bgen 
{
public:
	
	Bgen(const char* fileName, const char* indiFilePath);
	// takes individuals from individuals file

	int getN();
	float getFrequency();
	void calcNextSnpDosage();
	void skipSnps(int skipM);
    void rewindBgenFile();
	
	gsl_vector_float* genVec;
	
	void printGenProb();
	void printGenProb(int maxN);
	void printGenVec(int maxN);


// private:
	FILE* bgenFile;
	int nBgen; // number of individuals in the bgen file
	int mBgen; // the total number of SNPs in the BGEN file
    int skipBytesUzBuf; 

	int n; // number of individuals used
    gsl_vector_int* indiBgenPos;

	int currM; // current SNP

	unsigned char* zBuf; // buffer for compressed data
  	unsigned char* uzBuf; // buffer for unpacked data

  	float frequency; // frequency of B alleles in genVec

    void readFileHeader();
    void readSnpHeader(int &zipLen, int &unZipLen); // reads in SNP header and returns size of compressed and uncompressed SNP data; used by loadNextSnp()
	void dosageSnpFromBuffer();
};


inline int Bgen::getN()
{
	return n;
}

inline float Bgen::getFrequency()
{
	return frequency;
}

#endif
