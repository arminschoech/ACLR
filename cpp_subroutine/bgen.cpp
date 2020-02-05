#include "bgen.h"

Bgen::Bgen(const char* fileName, const char* indiFilePath)
{
	bgenFile = fopen(fileName, "rb");
	if (bgenFile == NULL) {
    	std::cerr << "ERROR: BGEN file cannot be opened." << std::endl;
    	exit(1);
  	}

    currM = 0;

    readFileHeader();
    skipBytesUzBuf = nBgen + 10;

    // read in individuals used
    FILE* indiFile = fopen(indiFilePath, "rb");
    if (indiFile == NULL) {
        std::cerr << "ERROR: " << indiFilePath << " can't be opened." << 
            std::endl;
        exit(1);
    }
    fseek(indiFile, 0, SEEK_END);
    n = ftell(indiFile) / 4;
    rewind(indiFile);
    indiBgenPos = gsl_vector_int_alloc(n);
    int currentNum;
    for(int i = 0; i < n; i++) {
        fread(&currentNum, 4, 1, indiFile);
        if(currentNum >= nBgen) {
            std::cerr << "ERROR: individual file includes index that exceeds BGEN file size." << std::endl;
            exit(1);
        }
        gsl_vector_int_set(indiBgenPos, i, currentNum);
    }
    // cout << "BGEN N: " << nBgen << endl;
    // cout << "Indi N: " << n << endl;

  	// declare variables and buffer
  	zBuf = (unsigned char *) malloc(3*nBgen + 10);
  	uzBuf = (unsigned char *) malloc(3*nBgen + 10);

  	// initialize genotype vector
  	genVec = gsl_vector_float_alloc(n);
}


void Bgen::readFileHeader() {
  	unsigned int offset; fread(&offset, 4, 1, bgenFile); 
  	unsigned int H; fread(&H, 4, 1, bgenFile);
  	unsigned int uIntM; fread(&uIntM, 4, 1, bgenFile); mBgen = (int) uIntM;
  	unsigned int uIntN; fread(&uIntN, 4, 1, bgenFile); nBgen = (int) uIntN;
    char magicBytes[5]; fread(magicBytes, 1, 4, bgenFile); magicBytes[4] = '\0'; // cout << "magic bytes: " << std::string(magicBytes) << endl;
    if (magicBytes != std::string("bgen")) {
        std::cerr << "ERROR: BGEN file format identification bytes are incorrect" << endl;
        exit(1);
    }
    unsigned int flags; fread(&flags, 4, 1, bgenFile);
     if (flags != 9 && flags != 2147483657) { // 2147483657 if sample IDs are included
        std::cerr << "ERROR: BGEN file does not have the expected format (either SNP blocks are not compressed or it is not in 'Layout 2' format)" << endl;
        exit(1);
    }
  	fseek(bgenFile, offset+4, SEEK_SET);
}


void Bgen::calcNextSnpDosage()
{
    if(currM >= mBgen) {
        std::cerr << "ERROR: accessing SNP number larger than BGEN file size" << std::endl;
        exit(1);
    }
	// read in genotype data
	int zLen, uzLen; 
    readSnpHeader(zLen, uzLen); 
		
	fread(zBuf, 1, zLen, bgenFile); //read in compressed data
	uLongf destLen = 3*nBgen + 10;
   	if (uncompress((Bytef *) uzBuf, &destLen, zBuf, zLen) != Z_OK || destLen != 3 * (unsigned int) nBgen + 10)
    {
        std::cerr << "ERROR: uncompress() failed" << std::endl;
      		exit(1);
   	}

    unsigned char nOfSnpChars[4];
    for (int i = 0; i < 4; i++) {
        nOfSnpChars[i] = uzBuf[i];
    }
    int nOfSnp = *(int *)nOfSnpChars;
    if (nOfSnp != nBgen) {
        cerr << "ERROR: individual number at SNP " << currM << " does not match the number in BGEN file header" << endl;
        exit(1);
    }

    if ((int) uzBuf[4] != 2) {
        cerr << "ERROR: SNP " << currM << " is not bi-allelic" << endl;
        exit(1);
    }

    if ((int) uzBuf[6] != 2 || (int) uzBuf[7] != 2) {
        cerr << "ERROR: SNP " << currM << " maximum or minimum ploidy differs from 2" << endl;
        exit(1);
    }

	dosageSnpFromBuffer();
   	currM++;
}


void Bgen::readSnpHeader(int &zipLen, int &unZipLen)
{
    unsigned short LS; fread(&LS, 2, 1, bgenFile); 
    char snpID[1000]; fread(snpID, 1, LS, bgenFile); snpID[LS] = '\0'; // cout << "snpID: " << snpID << endl;
    unsigned short LR; fread(&LR, 2, 1, bgenFile); 
    char rsID[1000]; fread(rsID, 1, LR, bgenFile); rsID[LR] = '\0'; // cout << "rsID: " << rsID << endl;
    unsigned short LC; fread(&LC, 2, 1, bgenFile); 
    char chrStr[1000]; fread(chrStr, 1, LC, bgenFile); chrStr[LC] = '\0';
    int chrom;
    if (sscanf(chrStr, "%d", &chrom)!=1 || !(1<=chrom&&chrom<=22)) {
      std::cerr << "ERROR: Invalid chrom (expecting integer 1-22): " << std::string(chrStr) << endl;
      exit(1);
    }
    unsigned int physpos; fread(&physpos, 4, 1, bgenFile); // cout << "physpos: " << physpos << endl;
    unsigned short alleleNum; fread(&alleleNum, 2, 1, bgenFile); // cout << "alleleNum: " << alleleNum << endl;
    if (alleleNum != 2U) {
        std::cerr << "ERROR: " << snpID << " has " << alleleNum << " alleles (non-bi-allelic)" << endl;    
        exit(1);
    }
    unsigned int LA; fread(&LA, 4, 1, bgenFile);
    char allele1[1000]; fread(allele1, 1, LA, bgenFile); allele1[LA] = '\0';
    // cout << "Allele 1: " << allele1 << endl;
    unsigned int LB; fread(&LB, 4, 1, bgenFile);
    char allele0[1000]; fread(allele0, 1, LB, bgenFile); allele0[LB] = '\0';
    // cout << "Allele 0: " << allele0 << endl;
    unsigned int uZipLen; fread(&uZipLen, 4, 1, bgenFile);
    zipLen = (int) uZipLen - 4; // 4 bytes for unZipLen subtracted
    // cout << "zipLen: " << zipLen << endl;
    unsigned int uUnZipLen; fread(&uUnZipLen, 4, 1, bgenFile);
    unZipLen = (int) uUnZipLen; // cout << "unZipLen: " << unZipLen << endl;
}


void Bgen::skipSnps(int skipM)
{
    if(currM + skipM > mBgen) {
        std::cerr << "ERROR: accessing SNP number larger than BGEN file size" << std::endl;
        exit(1);
    }

	for(int i = 0; i < skipM; i++)
	{
		// read in genotype data
		int zLen, uzLen; 
        readSnpHeader(zLen, uzLen); 
		fseek(bgenFile, zLen, SEEK_CUR);
		currM++;
	}
}


void Bgen::rewindBgenFile()
{
    rewind(bgenFile);
    readFileHeader();
    currM = 0;
}


void Bgen::dosageSnpFromBuffer()
{
	float AA;
	float AB;
	float BB;
	float genoSum = 255.0; // since probabilities are stored in single bytes

	double dosageSum = 0;
	float dosage; 

	for(int i = 0; i < n; i++)
	{
		// look up the BGEN file position of the ith used individual
		// and get the genotype probabilities there
		AA = float(uzBuf[2*gsl_vector_int_get(indiBgenPos, i) + skipBytesUzBuf]);
		AB = float(uzBuf[2*gsl_vector_int_get(indiBgenPos, i) + skipBytesUzBuf + 1]);
		BB = genoSum - AA - AB;
		
		dosage = (AB + 2.0*BB) / genoSum;
		gsl_vector_float_set(genVec, i, dosage); 
		dosageSum += dosage; 
	}
	frequency = dosageSum / 2.0 / (float) n;
}



void Bgen::printGenProb(int maxN)
{
	std::cout << "Bgen probability matrix:" << std::endl;
	for(int i = 0; i < maxN; i++)
	{
		std::cout << uzBuf[3*i]/32768.0 << '\t' << uzBuf[3*i+1]/32768.0 << '\t' << uzBuf[3*i+2]/32768.0 << std::endl;
	}
}

void Bgen::printGenProb()
{
	std::cout << "Bgen probability matrix:" << std::endl;
	for(int i = 0; i < n; i++)
	{
		std::cout << uzBuf[3*i]/32768.0 << '\t' << uzBuf[3*i+1]/32768.0 << '\t' << uzBuf[3*i+2]/32768.0 << std::endl;
	}
}

void Bgen::printGenVec(int maxN)
{
	std::cout << "Genotype vector:" << std::endl;
	for(int i = 0; i < maxN; i++) {
		std::cout << gsl_vector_float_get(genVec, i) << std::endl;
	}
	std::cout << std::endl;
}



