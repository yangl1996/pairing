pairing: pairing.c include/mcl/bn_c384_256.h lib/libmclbn384_256.dylib lib/libmcl.dylib
	clang -Iinclude -Llib -o pairing -lmclbn384_256 -lmcl -Ofast pairing.c

lib/libmclbn384_256.dylib lib/libmcl.dylib include/mcl/bn_c384_256.h:
	cd mcl && make
	cp -r mcl/lib .
	cp -r mcl/include .

.PHONY: clean
clean:
	rm -f pairing
