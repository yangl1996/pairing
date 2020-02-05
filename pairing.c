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
	clock_t start_commitment, end_commitment;
	start_commitment = clock();
	PCcommit(&C, &srs, data, 2048);
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
		mclBnFr EvalI;
		// generate the witness
		PCwitness(W + i, &EvalI, i, data, 2048, &srs);

		// verify witness
		clock_t start_witness, mid_witness, end_witness;
		mclBnGT e1;
		start_witness = clock();
		PCverifyEval_computeCG2(&e1, &C, &pc);
		mid_witness = clock();
		int eq = PCverifyEval_precomputed(i, &EvalI, W + i, &e1, &pc);
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

