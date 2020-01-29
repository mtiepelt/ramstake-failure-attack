#!/usr/bin/env bash

numfailures=2048        # should be at least 2048
numsamples=64          # should be at least 64
seedsamples="c0ffee"
seedfailures="d118b7"
# 0d86a4
# 170784
# 195bc2
# c14532

echo "--------------------------------------------------------------------"
echo " Demonstator for exploiting decryption failures in the Ramstake KEM "
echo "--------------------------------------------------------------------"
echo "--------------------------------------------------------------------"
echo " "
echo " Precomputation phase"
echo "---------"
echo " Estimating parameters with $numsamples samples..."
cd Ramstake_RS_756839/XKCP
make generic64/libkeccak.a
cd ..
make samples
cd ..
Ramstake_RS_756839/samples $numsamples $seedsamples
echo "----------------------------------"
echo " "
echo " Attack phase: querying decryption oracle"
echo "---------"
echo " Generating failures with seed: $seedfailures"
cd Ramstake_RS_756839
make getfailures
cd ..
Ramstake_RS_756839/getfailures $numfailures $seedfailures
echo "---------"
echo " Estimating secrets..."
python estimate.py -s failures.txt -p samplesfail.txt samplessucces.txt
echo "---------"
echo " Partitioning with secrets as additional input to verify partitions."
python stagedPartitioning.py -e estimate-a.npy estimate-b.npy -s secret-a.npy secret-b.npy -m 
echo "---------"

