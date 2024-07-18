# Demonstrator for exploiting decryption failures in the Ramstake KEM

This code package contains implementations in python to demonstrate an attack on the Ramtake KEM using decryption failures.

## Dependencies

* eXtended Keccak Code Package [XKCP](https://github.com/XKCP/XKCP)
  * Tested with commit `3c43b02f7e120bdb942810fd3a29ba65c79ce84f` from Nov 16 2023  
* [Ramstake] (https://github.com/0x64616E69656C/ramstake)
  * Uses local copy, **single line in ramstake.h was changed**
* Python 3.10.12, with libraries 
  * [peakdetect](https://github.com/avhn/peakdetect
) Version 1.2 
  * [pathos](https://github.com/uqfoundation/pathos) Version 0.3.2 
  * [tqdm](https://tqdm.github.io) Version 4.65.0 

## Files

### Run
|-- run.Unix.sh -- main file to be executed.

Note that execution may take several days.

Adjusting *numfailures* and *numsamples* may decrease execution time, but also significantly decreases quality of results.

> numfailures=1024

> numsamples=128

Seeds for pseudorandom generation of parameters may be adjusted arbitrarily:
> seedsamples="c0ffee"

> seedfailures="d118b7"


### Samples 

samples/

|-- estimates_params.npy

|-- SEED/estimate-a.npy

|-- SEED/estimate-b.npy

|-- SEED/secret-a.npy

|-- SEED/secret-b.npy

Contains the estimated parameters from the precomputation phase and precomputed estimates from decryption failures generated from the respective seeds.

### Other Files

| File | Description | Comment |
| ----------- |:-----------:|--------:|
| Ramstake_RS_756839/ramstake.h | Modified header file using only a single codeword. | |
| Ramstake_RS_756839/samples.c | Collects random failing and succeeding ciphertexts. | called by run.Unix.sh |
| Ramstake_RS_756839/getfailures.c | Collects decryption failures for target public key. | called by run.Unix.sh |
| estimate.py | Estimation of secrets. | Called by run.Unix.sh |
| stagedPartitioning.py | Partitioning using the estimation and evaluation of partition quality. | called by run.Unix.sh |
| checkPartition.py | Evaluates quality of intervals resulting from partitioning. | called by stagedPartitioning.py |
| plotGraphs.py | Constructs the graphs of estimates provided in the paper. | called by stagedPartitioning.py |
| plotResults.py | Constructs the results graph in the paper.  |  |
