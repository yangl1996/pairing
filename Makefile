bn: pairing.c include/mcl/bn_c256.h lib/libmclbn256.dylib lib/libmcl.dylib
	clang -DBN254 -Iinclude -Llib -o pairing -lmclbn256 -lmcl -Ofast pairing.c polycommit.c

bls: pairing.c include/mcl/bn_c384_256.h lib/libmclbn384_256.dylib lib/libmcl.dylib
	clang -DBLS12_381 -Iinclude -Llib -o pairing -lmclbn384_256 -lmcl -Ofast pairing.c

lib/libmclbn384_256.dylib lib/libmcl.dylib lib/libmclbn256.dylib:
	cd mcl && make
	cp -r mcl/lib .
	cp -r mcl/include .

.PHONY: clean
clean:
	rm -f pairing polycommit.o
