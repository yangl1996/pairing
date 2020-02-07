#ifdef BLS
#include "polycommit_bls.h"
#define RSIZE 32	// size of Fr in bytes
#endif
#ifdef BN
#include "polycommit_bn.h"
#define RSIZE 32
#endif
#include "string.h"
#include "stdio.h"

int PCsrs_init(PCsrs* srs, const char* G1sk, const char* G2sk, const char* Ask,
               int srs_len) {
	// init G1, G2, alpha from input strings
	if (mclBnG1_hashAndMapTo(&srs->G1, G1sk, strlen(G1sk))) {
		return 1;
	}
	if (mclBnG2_hashAndMapTo(&srs->G2, G2sk, strlen(G2sk))) {
		return 1;
	}
	mclBnFr A;
	if (mclBnFr_setHashOf(&A, Ask, strlen(Ask))) {
		return 1;
	}
	srs->srs_len = srs_len;

	srs->G1PK = malloc(srs_len * sizeof(mclBnG1));
	srs->G2PK = malloc(srs_len * sizeof(mclBnG2));
	srs->G1PK[0] = srs->G1;
	srs->G2PK[0] = srs->G2;
	for (int i = 1; i < srs_len; i++) {
		// compute G1^A^i ad G2^A^i
		mclBnG1_mul(srs->G1PK + i, srs->G1PK + i - 1, &A);
		mclBnG2_mul(srs->G2PK + i, srs->G2PK + i - 1, &A);
	}
	return 0;
}

int PCprecompute_init(PCprecompute* pc, const PCsrs* srs, int eval_len) {
	mclBn_pairing(&pc->eG1G2, &srs->G1, &srs->G2);			// e(G1, G2)
	pc->mG2AG2I = (uint64_t**)malloc(eval_len * sizeof(uint64_t*));
	pc->mG2 = (uint64_t*)malloc(mclBn_getUint64NumToPrecompute() * sizeof(uint64_t));
	mclBn_precomputeG2(pc->mG2, &srs->G2);
	mclBnG2 mulres;
	mclBnG2 G2AdivG2I;
	mclBnFr I;
	for (int i = 0; i < eval_len; i++) {
		mclBnFr_setInt(&I, i);
		mclBnG2_mul(&mulres, &srs->G2, &I);			// G2^I
		mclBnG2_sub(&G2AdivG2I, srs->G2PK + 1, &mulres);	// G2^A/G2^I
		pc->mG2AG2I[i] = (uint64_t*)malloc(mclBn_getUint64NumToPrecompute() * sizeof(uint64_t));
		mclBn_precomputeG2(pc->mG2AG2I[i], &G2AdivG2I);
	}
	int srs_len = srs->srs_len;
	pc->expG1 = (mclBnG1***)malloc(srs_len * sizeof(mclBnG1**));
	for (int i = 0; i < srs_len; i++) {
		// precompute for G1^A^i
		pc->expG1[i] = (mclBnG1**)malloc(RSIZE * sizeof(mclBnG1*));
		// set the base for k = 2^(j*8), we start with j=0
		mclBnFr expBase, expStep;
		mclBnFr_setInt(&expBase, 1);
		mclBnFr_setInt(&expStep, 1 << 8);
		for (int j = 0; j < RSIZE; j++) {
			// precompute for the j-th least significant byte
			mclBnG1 base;
			mclBnG1_mul(&base, srs->G1PK + i, &expBase);
			pc->expG1[i][j] = (mclBnG1*)malloc((1<<8) * sizeof(mclBnG1));
			mclBnG1_clear(pc->expG1[i][j] + 0);
			for (int k = 1; k < (1<<8); k++) {
				mclBnG1_add(pc->expG1[i][j] + k, pc->expG1[i][j] + k - 1, &base);
			}
			// bump up the exp base
			mclBnFr_mul(&expBase, &expBase, &expStep);
		}
	}
	return 0;
}

inline int PCmulVecG1Precomputed(mclBnG1* c, const PCprecompute* pc, const mclBnFr* poly, int len) {
	// C = Prod(G1PK[i] ^ poly[i])
	mclBnG1_clear(c);
	uint8_t ef[RSIZE];
	for (int i = 0; i < len; i++) {
		mclBnG1 cthis;
		mclBnG1_clear(&cthis);
		// first decompose the coefficient into bit string
		int bufsize = mclBnFr_getLittleEndian(ef, RSIZE, poly + i);
		// then add them up. notice that ef is little endian
		for (int j = 0; j < bufsize; j++) {
			mclBnG1_add(&cthis, &cthis, pc->expG1[i][j] + ef[j]);
		}
		mclBnG1_add(c, c, &cthis);
	}
	return 0;
}

int PCwitness(mclBnG1* w, mclBnFr* evalRes, int evalPoint, const mclBnFr* poly, int poly_len, const PCsrs* srs, const PCprecompute* pc) {
	// first, evaluate the polynomial
	// i = evalPoint
	mclBnFr I;
	mclBnFr_setInt(&I, evalPoint);
	mclBn_FrEvaluatePolynomial(evalRes, poly, poly_len, &I);
	// then, calculate (Poly(x)-Poly(i))/(x-i)
	mclBnFr* newPoly;
	newPoly = (mclBnFr*)malloc(poly_len * sizeof(mclBnFr));
	for (int k = 0; k < poly_len; k++) {
		// coefficient for x^k is Poly_t i^t-k-1 + Poly_t-1 i^t-k-2 + ... + Poly_k+1 i^0
		// where t is degree of the polynomial
		mclBnFr_clear(newPoly + k);
		mclBn_FrEvaluatePolynomial(newPoly + k, poly + k + 1, poly_len - 1 - k, &I);
	}
	// finally, calculate the witness
	PCmulVecG1Precomputed(w, pc, newPoly, poly_len);

	free(newPoly);
	return 0;
}

int PCcommit(mclBnG1* c, const mclBnFr* poly, int len, const PCsrs* srs, const PCprecompute* pc) {
	PCmulVecG1Precomputed(c, pc, poly, len);
	return 0;
}
/*
int PCbatchWitness(mclBnG1* w, mclBnFr* r, mclBnFr* evalRes, const int* evalPoint, int evalLen, const mclBnFr* poly,
                   int poly_len, const PCsrs* srs) {
	// first, evaluate at all the requested points
	mclBnFr* I;
	I = (mclBnFr*)malloc(evalLen * sizeof(mclBnFr));
	for (int i = 0; i < evalLen; i++) {
		mclBnFr_setInt(I + i, evalPoint[i]);
		mclBn_FrEvaluatePolynomial(evalRes + i, poly, poly_len, I + i);
	}
	// run interpolation to get r(x), b/c r(x)=poly(x) for all x in evalPoint
	// r(x) has degree evalLen - 1. b/c if r(x) has degree evalLen, r(x) can be devided by Prod((x-poly(i))
	// TODO: potential improvement http://fourier.eng.hmc.edu/e176/lectures/ch7/node3.html
	//
	// the library function only gives y(0), so we use the following method:
	// first, get r_0=y(0). then r(x) = x (r_t x^t-1 + ... + r_1) + r_0
	// so (r(x) - r_0) / x = r_t x^t-1 + ... + r_1, so the same interpolation for r_1, r_2...r_t
	const mclBnFr* x = I;
	mclBnFr* y;
	y = (mclBnFr*)malloc(evalLen * sizeof(mclBnFr));
	memcpy(y, evalRes, sizeof(mclBnFr) * evalLen);
	for (int i = 0; i < evalLen; i++) {
		int remaining_uk = evalLen - i;
		mclBn_FrLagrangeInterpolation(r + i, x, y, remaining_uk);
		for (int j = 0; j < evalLen; j++) {
			// y[j] <- (y[j] - r_i) / x[j])
			mclBnFr_sub(y + j, y + j, r + i);
			mclBnFr_div(y + j, y + j, x + j);
		}
	}

	free(I);
	free(y);
	return 0;
}
*/
