binary: demo_bn demo_bls bench_bn bench_bls

static: libpolycommitbn.a libpolycommitbls.a

libpolycommitbn.a: polycommit.c
	clang -DBN -Imcl/include -o polycommit_bn.o -Ofast -Wall -c polycommit.c
	ar rcs libpolycommitbn.a polycommit_bn.o

libpolycommitbls.a: polycommit.c
	clang -DBLS -Imcl/include -o polycommit_bls.o -Ofast -Wall -c polycommit.c
	ar rcs libpolycommitbls.a polycommit_bls.o

demo_bn: demo.c libpolycommitbn.a
	clang -DBN -Imcl/include -Lmcl/lib -L. -o demo_bn -Ofast -Wall demo.c -lpolycommitbn -lmcl -lmclbn256 

demo_bls: demo.c libpolycommitbls.a
	clang -DBLS -Imcl/include -Lmcl/lib -L. -o demo_bls -Ofast -Wall demo.c -lpolycommitbls -lmcl -lmclbn384_256 

bench_bn: bench.c libpolycommitbn.a
	clang -DBN -Imcl/include -Lmcl/lib -L. -o bench_bn -Ofast -Wall bench.c -lpolycommitbn -lmcl -lmclbn256 

bench_bls: bench.c libpolycommitbls.a
	clang -DBLS -Imcl/include -Lmcl/lib -L. -o bench_bls -Ofast -Wall bench.c -lpolycommitbls -lmcl -lmclbn384_256 

.PHONY: clean
clean:
	rm -f *.o *.a *.so demo_bn demo_bls bench_bn bench_bls
