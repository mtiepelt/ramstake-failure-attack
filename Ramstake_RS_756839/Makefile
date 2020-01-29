CC=gcc
OPTFLAGS=-O3
CFLAGS=-I XKCP/bin/generic64/ # point this and the one below to wherever the keccak code package is located
LFLAGS=-L XKCP/bin/generic64/ -lgmp -lkeccak -lcrypto
TEST_HEADERS = csprng.h  gf256x.h  ramstake.h  reedsolomon.h 
FAILURES_SOURCE = csprng.c  gf256x.c  ramstake.c  reedsolomon.c  getfailures.c codec_rs.c 
SAMPLES_SOURCE = csprng.c  gf256x.c  ramstake.c  reedsolomon.c  samples.c codec_rs.c 
PERFORM_SOURCE = csprng.c  gf256x.c  ramstake.c  reedsolomon.c  perform.c codec_rs.c 
KAT_HEADERS = csprng.h  gf256x.h  ramstake.h  reedsolomon.h  rng.h 
KAT_SOURCES = csprng.c  gf256x.c  ramstake.c  reedsolomon.c  PQCgenKAT_kem.c  rng.c  kem.c codec_rs.c 

getfailures: $(FAILURES_SOURCE) $(TEST_HEADERS)
	gcc -o $@ $(FAILURES_SOURCE) $(CFLAGS) $(LFLAGS) $(OPTFLAGS)

samples: $(SAMPLES_SOURCE) $(TEST_HEADERS)
	gcc -o $@ $(SAMPLES_SOURCE) $(CFLAGS) $(LFLAGS) $(OPTFLAGS)

perform: $(PERFORM_SOURCE) $(TEST_HEADERS)
	gcc -o $@ $(PERFORM_SOURCE) $(CFLAGS) $(LFLAGS) $(OPTFLAGS) -Wall -lrt

kat: $(KAT_HEADERS) $(KAT_SOURCES)
	gcc -o PQCgenKAT_kem $(KAT_SOURCES) $(CFLAGS) $(LFLAGS) $(OPTFLAGS) -std=c99 -DKAT

api: ramstake.h genapi.c
	gcc -o genapi genapi.c $(OPTFLAGS)
	./genapi > api.h
	rm -f genapi

clean:
	rm -f *.o
	rm -f getfailures
	rm -f samples
	rm -f perform
	rm -f kat
	rm -f genapi

