#include <stdio.h>
#include <string.h>

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
	mclBnG1*** expG1; // expG1[i][j][k] = G1^Ai^(k * 2^(8*j))
} PCprecompute;

int PCsrs_init(PCsrs* srs, const char* G1sk, const char* G2sk, const char* Ask,
               int srs_len);

int PCprecompute_init(PCprecompute* pc, const PCsrs* srs, int eval_len);

int PCcommit(mclBnG1* c, const mclBnFr* poly, int len, const PCsrs* srs, const PCprecompute* pc);

// poly is defined as the coefficients, from low degree to high degree
int PCwitness(mclBnG1* w, mclBnFr* evalRes, int evalPoint, const mclBnFr* poly, int poly_len, const PCsrs* srs, const PCprecompute* pc);

inline int PCverifyEval_computeCG2(mclBnGT* CG2, const mclBnG1* c, const PCprecompute* pc) {
	// CG2 = e(C, G2);
	mclBn_precomputedMillerLoop(CG2, c, pc->mG2);
	mclBn_finalExp(CG2, CG2);
	return 0;
}

inline int PCverifyEval_precomputed(int evalPoint, const mclBnFr* evalRes,
		                    const mclBnG1* w, const mclBnGT* CG2, const PCprecompute* pc) {
	mclBnGT e1, e2;
	// I = evalpoint
	// e1 = e(w, G2^A/G2^I)
	mclBn_precomputedMillerLoop(&e1, w, pc->mG2AG2I[evalPoint]);
	mclBn_finalExp(&e1, &e1);
	// e2 = e(G1, G2)^Poly(I)
	mclBnGT_pow(&e2, &pc->eG1G2, evalRes);
	// e = e1 e2
	mclBnGT e;
	mclBnGT_mul(&e, &e1, &e2);
	return mclBnGT_isEqual(CG2, &e);
}

inline int PCverifyEval(int evalPoint, const mclBnFr* evalRes,
		        const mclBnG1* c, const mclBnG1* w,
                        const PCprecompute* pc) {
	mclBnGT CG2;
	PCverifyEval_computeCG2(&CG2, c, pc);
	return PCverifyEval_precomputed(evalPoint, evalRes, w, &CG2, pc);
}

