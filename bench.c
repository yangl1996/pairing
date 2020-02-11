#include <stdio.h>
#include <time.h>
#ifdef BN
#include "polycommit_bn.h"
#endif
#ifdef BLS
#include "polycommit_bls.h"
#endif

#define DEG 2048
#define NUM_WITNESS 20

int main()
{
#ifdef BN
	int ret = mclBn_init(MCL_BN254, MCLBN_COMPILED_TIME_VAR);
#endif
#ifdef BLS
	int ret = mclBn_init(MCL_BLS12_381, MCLBN_COMPILED_TIME_VAR);
#endif
	if (ret != 0) {
		printf("error initialize crypto library, ret=%d\n", ret);
		return 1;
	}

	// init the SRS from hashes
	PCsrs srs;
	PCsrs_init(&srs, "g1sk", "g2sk", "alphask", DEG);

	// precompute stuff
	PCprecompute pc;
	PCprecompute_init(&pc, &srs, DEG * 2);

	printf("block size %.2lf kilobytes, %d evaluation points, %d checks per node\n", DEG * 31.0 / 1024.0, DEG * 2, NUM_WITNESS);

	// generate the polynomial to be encoded
	// the coefficients of the polynomial are the message chunks
	for (int test = 0; test < 10; test++) {
		mclBnFr* data;
		data = (mclBnFr*)malloc(DEG * sizeof(mclBnFr));
		for (int i = 0; i < DEG; i++) {
			mclBnFr_setByCSPRNG(data + i);
		}

		// generate commitment
		mclBnG1 C;
		clock_t start_commitment, end_commitment;
		start_commitment = clock();
		PCcommit(&C, data, DEG, &srs, &pc);
		end_commitment = clock();
		double time_commitment;
		time_commitment = ((double) (end_commitment - start_commitment)) / CLOCKS_PER_SEC;

		// generate witnesses
		mclBnG1* W;
		mclBnFr* eval;
		W = (mclBnG1*)malloc(NUM_WITNESS* sizeof(mclBnG1));
		eval = (mclBnFr*)malloc(NUM_WITNESS * sizeof(mclBnFr));
		clock_t start_witness, end_witness;
		start_witness = clock();
		for (int wtns = 0; wtns < NUM_WITNESS; wtns++) {
			PCwitness(W + wtns, eval + wtns, wtns, data, DEG, &srs, &pc);
		}
		end_witness = clock();
		double time_witness;
		time_witness = ((double) (end_witness - start_witness)) / CLOCKS_PER_SEC / NUM_WITNESS;

		// verify witnesses
		int res = 0;
		clock_t start_verify, end_verify;
		start_verify = clock();
		mclBnGT e1;
		PCverifyEval_computeCG2(&e1, &C, &pc);
		for (int wtns = 0; wtns < NUM_WITNESS; wtns++) {
			res += !PCverifyEval_precomputed(wtns, eval + wtns, W + wtns, &e1, &pc);
		}
		end_verify = clock();
		double time_verify;
		time_verify = ((double) (end_verify - start_verify)) / CLOCKS_PER_SEC / NUM_WITNESS;
		if (res != 0) {
			printf("Verification Error\n");
		}
		printf("Commit=%.2lf ms (%.2lf Mbps); Eval=%.2lf ms (%.2lf Mbps); Check=%.2lf ms (%.2lf Mbps)\n",
                        time_commitment * 1000.0,
			DEG * 31.0 / 1024.0 * 8 / 1024.0 / time_commitment,
			time_witness * 1000.0,
			DEG * 31.0 / 1024.0 * 8 / 1024.0 / time_witness,
			time_verify * 1000.0,
			DEG * 31.0 / 1024.0 * 8 / 1024.0 / time_verify);
	}
	return 0;
}

