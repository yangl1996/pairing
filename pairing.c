#include <stdio.h>
#include <string.h>
#include <mcl/bn_c384_256.h>

int main()
{
	char buf[1600];
	int ret = mclBn_init(MCL_BLS12_381, MCLBN_COMPILED_TIME_VAR);
	if (ret != 0) {
		printf("error initialize crypto library, ret=%d\n", ret);
		return 1;
	}

	// define generators of G1 and G2 and init them from random hashes
	mclBnG1 G1;
	mclBnG2 G2;
	if (mclBnG1_hashAndMapTo(&G1, "g1blah", 6)) {
		printf("error init generator of G1\n");
		return 1;
	}
	if (mclBnG2_hashAndMapTo(&G2, "g2hey", 5)) {
		printf("error init generator of G2\n");
		return 1;
	}
	mclBnG1_getStr(buf, sizeof(buf), &G1, 16);
	printf("G1=%s\n", buf);
	mclBnG2_getStr(buf, sizeof(buf), &G2, 16);
	printf("G2=%s\n", buf);

	// define the secret key (alpha) and init from random hash
	mclBnFr A;
	if (mclBnFr_setHashOf(&A, "alphask", 7)) {
		printf("error init secret key\n");
		return 1;
	}
	mclBnFr_getStr(buf, sizeof(buf), &A, 16);
	printf("Secret key Alpha=%s\n", buf);

	// generate the public key from G1 and G2
	// we need to handle polynomials of degree t=2048
	mclBnG1* G1PK;
	mclBnG2* G2PK;
	G1PK = (mclBnG1*)malloc(2048 * sizeof(mclBnG1));
	G2PK = (mclBnG2*)malloc(2048 * sizeof(mclBnG2));
	mclBnG1 G1Mul;	// used to store G1^t, will be updated in the loop
	mclBnG2 G2Mul;
	G1Mul = G1;
	G2Mul = G2;
	for (int i = 0; i < 2048; i++) {
		G1PK[i] = G1Mul;
		G2PK[i] = G2Mul;
		mclBnG1_mul(&G1Mul, &G1Mul, &A);
		mclBnG2_mul(&G2Mul, &G2Mul, &A);
	}
	printf("Public keys generated\n");

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
	mclBnG1 CItem;
	for (int i = 0; i < 2048; i++) {
		mclBnG1_mul(&CItem, G1PK + i, data + i);
		mclBnG1_add(&C, &C, &CItem);
	}
	mclBnG1_getStr(buf, sizeof(buf), &C, 16);
	printf("Commitment C=%s\n", buf);

	// generate witnesses for each evaluation point (suppose there are 4096)
	mclBnG1* W;
	W = (mclBnG1*)malloc(4096 * sizeof(mclBnG1));
	for (int i = 0; i < 4096; i++) {
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
			mclBnFr_clear(PDiv + i);
			mclBnFr PTemp;
			// the coefficient for x^k is Phi_t i^t-k-1 + Phi_t-1 i^t-k-2 + ... + Phi_k+1 i^0
			// t = 2047 here (degree of the polynomial)
			for (int j = 0; j < 2047 - k; j++) {
				mclBnFr_mul(&PTemp, data + (k+1+j), IExp + j);	// calculate Phi_(k+1+j) I^j
				mclBnFr_add(PDiv + k, PDiv + k, &PTemp);	// accumulate to PDiv[k]
			}
		}

		// now we can calculate the witness
		mclBnG1_clear(W + i);
		mclBnG1 WItem;
		for (int k = 0; k < 2048; k++) {
			mclBnG1_mul(&WItem, G1PK + k, PDiv + k);
			mclBnG1_add(W + i, W + i, &WItem);
		}
		mclBnG1_getStr(buf, sizeof(buf), W + i, 16);
		printf("Witness W[%d]=%s\n", i, buf);

		// verify witness
		// first calculate the evaluation at I
		mclBnFr EvalI, EvalItem;
		for (int j = 0; j < 2048; j++) {
			mclBnFr_mul(&EvalItem, IExp + j, data + j);
			mclBnFr_add(&EvalI, &EvalI, &EvalItem);
		}
		mclBnGT e1, e2_1, e2_2, e2;
		mclBn_pairing(&e1, &C, G2PK + 0);
		mclBnG2 GI, GAdivGI;
		mclBnG2_mul(&GI, G2PK + 0, &I);
		mclBnG2_sub(&GAdivGI, G2PK + 1, &GI);
		mclBn_pairing(&e2_1, W + i, &GAdivGI);
		mclBn_pairing(&e2_2, G1PK + 0, G2PK + 0);
		mclBnGT_pow(&e2_2, &e2_2, &EvalI);
		mclBnGT_mul(&e2, &e2_1, &e2_2);
		int eq = mclBnGT_isEqual(&e1, &e2);
		if (eq) {
			printf("... verification passed\n");
		} else {
			printf("!!! verification failed\n");
			mclBnGT_getStr(buf, sizeof(buf), &e1, 16);
			printf("E1=%s\n", buf);
			mclBnGT_getStr(buf, sizeof(buf), &e2, 16);
			printf("E2=%s\n", buf);
		}
	}

	return 0;
}

