#include <stdio.h>
#include <string.h>
#ifdef BN254
#include <mcl/bn_c256.h>
#endif
#ifdef BLS12_381
#include <mcl/bn_c384_256.h>
#endif

typedef struct PCsrs {
	mclBnG1 G1;
	mclBnG2 G2;
	mclBnG1* G1PK;
	mclBnG2* G2PK;
	int srs_len;
} PCsrs;

int PCsrs_init(PCsrs* srs, const char* G1sk, const char* G2sk, const char* Ask,
               int srs_len);
