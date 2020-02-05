#include <stdio.h>
#ifdef BN254
#include <mcl/bn_c256.h>
#endif
#ifdef BLS12_381
#include <mcl/bn_c384_256.h>
#endif
#include <time.h>
#include "polycommit.h"

int main()
{
	char buf[1600];
#ifdef BN254
	int ret = mclBn_init(MCL_BN254, MCLBN_COMPILED_TIME_VAR);
#endif
#ifdef BLS12_381
	int ret = mclBn_init(MCL_BLS12_381, MCLBN_COMPILED_TIME_VAR);
#endif
	if (ret != 0) {
		printf("error initialize crypto library, ret=%d\n", ret);
		return 1;
	}

	// init the SRS from hashes
	PCsrs srs;
	PCsrs_init(&srs, "g1sk", "g2sk", "alphask", 2048);
	mclBnG1 G1;
	mclBnG2 G2;
	G1 = srs.G1;
	G2 = srs.G2;

	mclBnG1_getStr(buf, sizeof(buf), &G1, 16);
	printf("G1=%s\n", buf);
	mclBnG2_getStr(buf, sizeof(buf), &G2, 16);
	printf("G2=%s\n", buf);

	mclBnG1* G1PK;
	mclBnG2* G2PK;
	G1PK = srs.G1PK;
	G2PK = srs.G2PK;
	printf("Public keys generated\n");

	// precompute stuff
	PCprecompute pc;
	PCprecompute_init(&pc, &srs, 4096);
	uint64_t** Qbuf;
	Qbuf = pc.mG2AG2I;
	mclBnGT eG1G2;
	eG1G2 = pc.eG1G2;
	uint64_t* Pbuf;
	Pbuf = pc.mG2;
	printf("millerloop(P, G2), millerloop(P, G2^A / G2^I), and e(G1, G2) precomputed\n");

	// generate the polynomial to be encoded
	// the coefficients of the polynomial are the message chunks
	mclBnFr* data;
	data = (mclBnFr*)malloc(2048 * sizeof(mclBnFr));
	for (int i = 0; i < 2048; i++) {
		mclBnFr_setByCSPRNG(data + i);
	}
	printf("Input polynomial generated\n");

	// generate commitment
	mclBnG1 C;
	mclBnG1_clear(&C);
	mclBnG1* CItems;
	clock_t start_commitment, end_commitment;
	start_commitment = clock();
	mclBnG1_mulVec(&C, G1PK, data, 2048);	// C = G1PK[] . data[]
	end_commitment = clock();
	double time_used;
	time_used = ((double) (end_commitment - start_commitment)) / CLOCKS_PER_SEC;
	mclBnG1_getStr(buf, sizeof(buf), &C, 16);
	printf("Commitment generated, time=%.2lfms, throughput=%.2lfMbps\n", time_used * 1000, 64 / 1024.0 * 8 / time_used);

	// generate witnesses for each evaluation point (suppose there are 4096)
	mclBnG1* W;
	W = (mclBnG1*)malloc(4096 * sizeof(mclBnG1));
	for (int i = 0; i < 4096; i++) {
		printf("Witness for i=%d...", i);
		fflush(stdout);
		mclBnFr I;
		mclBnFr_setInt(&I, i);
		mclBnFr* IExp;	// IExp[m] = I^m
		IExp = (mclBnFr*)malloc(2048 * sizeof(mclBnFr));
		mclBnFr_setInt(IExp + 0, 1);
		for (int m = 1; m < 2048; m++) {
			mclBnFr_mul(IExp + m, IExp + m - 1, &I);
		}

		// first, we need to calculate (Phi(x)-Phi(i))/(x-i) for the given i
		mclBnFr* PDiv;	// coefficients
		PDiv = (mclBnFr*)malloc(2048 * sizeof(mclBnFr));
		for (int k = 0; k < 2048; k++) {
			mclBnFr_clear(PDiv + k);
			mclBnFr PTemp;
			// the coefficient for x^k is Phi_t i^t-k-1 + Phi_t-1 i^t-k-2 + ... + Phi_k+1 i^0
			// t = 2047 here (degree of the polynomial)
			for (int j = 0; j < 2047 - k; j++) {
				mclBnFr_mul(&PTemp, data + (k+1+j), IExp + j);	// calculate Phi_(k+1+j) I^j
				mclBnFr_add(PDiv + k, PDiv + k, &PTemp);	// accumulate to PDiv[k]
			}
		}

		// verify the calculation is good: evaluate (Phi(x)-Phi(i))/(x-i) at random point ep
		mclBnFr ep, er, ei, ex, epExp;
		mclBnFr_clear(&er);	// er: Phi(ep)-Phi(i)/ep-i
		mclBnFr_clear(&ei);	// ei: Phi(i)
		mclBnFr_clear(&ex);	// ex: Phi(ep)
		mclBnFr_setByCSPRNG(&ep);
		mclBnFr_setInt(&epExp, 1);
		for (int m = 0; m < 2048; m++) {
			mclBnFr temp;
			mclBnFr_mul(&temp, PDiv + m, &epExp);
			mclBnFr_add(&er, &er, &temp);
			
			mclBnFr_mul(&temp, data + m, &epExp);
			mclBnFr_add(&ex, &ex, &temp);

			mclBnFr_mul(&temp, data + m, IExp + m);
			mclBnFr_add(&ei, &ei, &temp);

			mclBnFr_mul(&epExp, &epExp, &ep);
		}
		mclBnFr epsubi, exsubei, er2;
		mclBnFr_sub(&epsubi, &ep, &I);	// epsubi = ep - i
		mclBnFr_sub(&exsubei, &ex, &ei);// exsubei = Phi(ep) - Phi(i)
		mclBnFr_div(&er2, &exsubei, &epsubi);
		int eval_eq = mclBnFr_isEqual(&er, &er2);
		if (eval_eq) {
			printf("P(x)-P(i)/x-i ✔...");
		} else {
			printf("P(x)-P(i)/x-i ❌...");
		}
		fflush(stdout);

		// now we can calculate the witness
		mclBnG1_mulVec(W + i, G1PK, PDiv, 2048);	// W[i] = G1PK[] . PDiv[]

		// verify witness
		// first calculate the evaluation at I
		mclBnFr EvalI, EvalItem;
		mclBnFr_clear(&EvalI);
		for (int j = 0; j < 2048; j++) {
			mclBnFr_mul(&EvalItem, IExp + j, data + j);
			mclBnFr_add(&EvalI, &EvalI, &EvalItem);
		}
		clock_t start_witness, mid_witness, end_witness;
		mclBnGT e1, e2_1, e2_2, e2;
		start_witness = clock();
		mclBn_precomputedMillerLoop(&e1, &C, Pbuf);
		mclBn_finalExp(&e1, &e1);			// e1 = e(C, G2)
		mid_witness = clock();
		mclBn_precomputedMillerLoop(&e2_1, W + i, Qbuf[i]);
		mclBn_finalExp(&e2_1, &e2_1);			// e2_1 = e(W_i, G2^A / G2^I)
		mclBnGT_pow(&e2_2, &eG1G2, &EvalI);		// e2_2 = e(G1, G2)^Phi(I)
		mclBnGT_mul(&e2, &e2_1, &e2_2);
		int eq = mclBnGT_isEqual(&e1, &e2);
		end_witness = clock();
		double time_witness, time_unique;
		time_witness = ((double) (end_witness - start_witness)) / CLOCKS_PER_SEC;
		time_unique = ((double) (end_witness - mid_witness)) / CLOCKS_PER_SEC;
		if (eq) {
			printf("verification ✔ in %.2lfms (%.2lfms nonshared)\n", time_witness * 1000, time_unique * 1000);
		} else {
			printf("verification ❌\n");
		}
	}

	return 0;
}

