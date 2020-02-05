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
	int srs_len;	// support polynomial with degree up to srs_len-1
} PCsrs;

typedef struct PCprecompute {
	int eval_len;	// support evaluation with x=0..=eval_len-1
	mclBnGT eG1G2;	// e(G1, G2)
	uint64_t* mG2;	// for millerloop(P, G2)
	uint64_t** mG2AG2I;	// millerloop(P, G2^A/G2^I)
} PCprecompute;

int PCsrs_init(PCsrs* srs, const char* G1sk, const char* G2sk, const char* Ask,
               int srs_len);

int PCprecompute_init(PCprecompute* pc, const PCsrs* srs, int eval_len);
