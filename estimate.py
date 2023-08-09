import argparse
import numpy as np
import os
from tqdm import tqdm

"""
    --------------------------------
    HELPER FUNCITONS 
    --------------------------------
"""

# get binary representation with LSB first
def tobin(n, modulus_bitsize):
    n = bin(long(n))[2:]
    addzero = modulus_bitsize - len(n)
    binary = '0' * addzero + n
    binary = [int(i) for i in binary]
    return binary[::-1]

def mult(a, pos, modulus_bitsize):
    tmp = a+a
    pos = modulus_bitsize-pos
    return tmp[pos:pos+modulus_bitsize]


### data loading functions
def getestimatorsamples(location):
    with open(location) as f:
        # get rid of first line and comments
        f.next()
        for line in f:
            tmp = line.split(' ')
            yield int(tmp[0]), int(tmp[1]), int(tmp[2]), int(tmp[3])


def getsamples(location, modulus_bitsize):
    with open(location) as f:
        # get rid of first line
        f.next()
        secret = []
        for i in range(0, 3):
            tmp = int(f.next())
            secret.append(tobin(tmp, modulus_bitsize))
        yield secret

        # get rid of errors line
        f.next()
        for line in f:
            tmp = line.split(' ')
            yield int(tmp[0]), int(tmp[1])


### interesting functions :)
def geterrors(n, rs_codeword_length, modulus_bitsize, pos=0):
    # get the full error bytes
    # stra = tobin(n, modulus_bitsize)

    num_errors = 0
    for j in range(0, rs_codeword_length / 8):
        if sum(n[pos + (8 * j): pos + (8 * j + 8)]) > 0:
            num_errors += 1

    return num_errors


def calculateestimatorInner(a, d, ai, modulus_bitsize, q, rs_codeword_length):
    """
        Compute propability that a certain number of errors occurs for each bit position in secret.
    """
    abin = tobin(a, modulus_bitsize)
    dbin = tobin(d, modulus_bitsize)
    
    # find all ones in a
    apos = [i for i, x in enumerate(abin) if x == 1]

    # find how likely a certain amount of errors is in d for a one in a
    for i in apos:
        tmp = mult(dbin, i, modulus_bitsize)
        idx = geterrors(tmp, rs_codeword_length, modulus_bitsize)
        ai[idx] = ai.get(idx, 0) + 1

    return ai

def calculateestimator(location, locationsucces, modulus_bitsize, q, rs_codeword_length):
    """
        Construct estimator for byte errors of failing ciphertext.
    """
    samplesfail = getestimatorsamples(location)
    samplessucces = getestimatorsamples(locationsucces)
    
    # calculate the error distribution of a one on failure samples
    aione = {}
    aisucces = {}
    for a, b, c, d in tqdm(samplesfail):
        # calculate two times, one time for a-d, one time for c-b
        aione = calculateestimatorInner(a, d, aione, modulus_bitsize, q, rs_codeword_length)
        aione = calculateestimatorInner(c, b, aione, modulus_bitsize, q, rs_codeword_length)


        # calculate the error distribution of a one on success samples
        a, b, c, d = samplessucces.next()

        # calculate two times, one time for a-d, one time for c-b
        aisucces = calculateestimatorInner(a, d, aisucces, modulus_bitsize, q, rs_codeword_length)
        aisucces = calculateestimatorInner(c, b, aisucces, modulus_bitsize, q, rs_codeword_length)

        # determine the keys in both estimates
        aikeys = aione.keys()
        keys = sorted([i for i in aisucces.keys() if i in aikeys])

        # calculate estimator
        norm = float(sum(aione.values())) / sum(aisucces.values())
        estimator = [float(aione[i]) / aisucces[i] / norm for i in keys]
        print
        print [(i, aione[i]) for i in keys]
        print [(i, aisucces[i]) for i in keys]
        print estimator

    return estimator


### use the estimator
def getestimate(sample, est, modulus_bitsize, q, rs_codeword_length):
    """
        
    """
    adapt = np.zeros(modulus_bitsize)
    samplebin = tobin(sample, modulus_bitsize)
    samplebin = samplebin + samplebin

    # for every bit in the secret
    for i in range(0, modulus_bitsize):
        # get the amount of errors of ai * d
        pos = modulus_bitsize - i
        errors = geterrors(samplebin, rs_codeword_length, modulus_bitsize, pos)
        # the amount of errors determines the estimate
        # adapt[i] = est.get(errors,1)
        adapt[i] = np.interp(errors, range(0, len(est)), est, right=1)

        # shift by one
        # samplebin = [samplebin[-1]] + samplebin[0:-1]

    return adapt


def main():
    """
        Hardcoded parameters for Ramstake
    """
    rs_codeword_length = 255 * 8
    modulus_bitsize = 756839
    q = 2**modulus_bitsize - 1

    # Parse arguments
    parser = argparse.ArgumentParser("Secret Estimator. Compute either estimator parameters or estimate secrets based on a set of failures (requires precomputed parameters).")
    parser.add_argument("-p", "--params", nargs=2, dest="params", metavar=("samplesfail", "samplessucces"), help="Compute estimator parameter using random decryption failures (filename) and successes (filename).")
    parser.add_argument("-s", "--secrets", nargs=1, dest="secrets", metavar=("estimator"), required=True, help="Compute estimates of the secrets using a set of decryption failures (filename).")
    args = parser.parse_args() 

    print "---------------------------------"
    print " modulus_bitsize: " + str(modulus_bitsize)
    print "---------------------------------"

    if args.params:
        print " Calculating estimations from samples..."
        est = calculateestimator(args.params[0], args.params[1], modulus_bitsize, q, rs_codeword_length)
        np.save('estimated_params', est)

    """
        If --params not given as argument, file 'estimates_params.npy' should exist is same directory 
    """
    if not os.path.exists('estimated_params.npy'):
        print('Parameters from sampling not yet estimates. Please specify -p sample_fail samples_success')
        return 0

    est = np.load('estimated_params.npy')

    """
        Estimate secrets.
    """
    if args.secrets:
        #estimate a and b
        # b=1; a=0
        for aorb in [1, 0]:
            samples = getsamples(args.secrets[0], modulus_bitsize)

            # skip the secret, save and then delete
            secret = samples.next()
            np.save('secret-a',secret[1])
            np.save('secret-b',secret[2])
            del secret

            # get a standard probability estimate of the zero's and ones
            prob = np.ones(modulus_bitsize)

            # loop over all the samples, and estimate the required change
            from pathos.multiprocessing import ProcessingPool as Pool
            pool = Pool(nodes=4)
            
            f = lambda x: getestimate(x[aorb], est, modulus_bitsize, q, rs_codeword_length)

            print 'start'
            for i in tqdm(pool.uimap(f, samples)):
                prob = prob * i
                np.save('estimate-tmp',prob)

            if aorb==1:
                tmp = 'b'
            else:
                tmp = 'a'
            np.save('estimate-'+tmp,prob)

if __name__ == '__main__':
    main()
