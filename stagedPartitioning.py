# -*- coding: utf-8 -*-
import numpy as np
from peakdetect import peakdetect
import scipy.signal as sc

from operator import itemgetter

import argparse
import math

import random

# custom modules
import plotGraphs as pG
import checkPartition as cP
 
class Partitioning:
    """
        Naming convention:
            estimates   -   lst mapping each bit position to Pr[b_1i = 1/ b_i = 0]
            peaks       -   all peaks from find_peaks
            bit_range   -   all bit ranges derived from peaks
            bit_range_good  -   bit ranges with larger trailing empty space
            bit_range_bad   -   bit range with smaller trailing empty space
    """

    def __init__(self, modulus_bitlength, hw_omega):
        self.modulus_bitlength = int(modulus_bitlength)
        self.hw_omega = int(hw_omega)

        self.dummy = True

        self.plot_dir = 'results'

    """
        ---------------------------------
        EXTRACTING PARTS
        ---------------------------------
    """
    def classifyBitRangeP(self, bit_range):
        """
            Takes a range of peaks and outputs a set of correct parts and a set of sampling ranges
            @input: [START, END, e]
            @return: [START, END, e, e_p]
        """
        good = []
        bad = []

        for i in range(len(bit_range)):
            triple = bit_range[i][:]
            s = bit_range[i][0]
            e = bit_range[i][1]
            empty = bit_range[i][2]

            if (e-s) <= empty:
                good.append(triple)
            else:
                triple.extend([0])

                if i > 0:
                    prev_s = bit_range[i-1][0]
                    prev_e = bit_range[i-1][1]
                    prev_empty = bit_range[i-1][2]

                    if prev_e - prev_s <= prev_empty:
                        new_prev = [prev_s, prev_e, prev_e - prev_s]
                        good[-1] = new_prev

                        pro_empty = min(s - prev_e - new_prev[2], empty)
                        bad.append([s, e, empty, pro_empty])

                    else:
                        bad.append(triple) 
                else:
                    bad.append(triple)

        return good, bad

    """
        STAGE I: find_peaks
    """
    def completePartsFromFindPeaks(self, peaks, peak_look_ahead):
        """
            return [START, END, e]
        """
        lst_peaks = list(peaks)
        lst_peaks.sort()
        bit_range = []
        
        """

        """
        i = 0
        while i < len(lst_peaks):

            peak = lst_peaks[i]
            first_peak = lst_peaks[i]

            """
                Identify close peaks
            """
            i += 1
            while i < len(lst_peaks) and (peak + peak_look_ahead) > lst_peaks[i]:
                peak = lst_peaks[i]
                i += 1
            # lst_peak[i] is not close

            # [START, END, EMPTY DISTANCE]
            if i >= len(lst_peaks):
                # last peak. take distance to first part
                bit_range.append([first_peak, peak, modulus_bitlength - peak + lst_peaks[0]])
            else:
                bit_range.append([first_peak, peak, lst_peaks[i] - peak])

        return bit_range

    def recomputeEmptySpace(self, bit_range):
        """
            ASSUME: bit_range is sorted

            Recomputes the empty space for each bit range
        """
        new_bit_range = []

        for i in range(len(bit_range)):

            s = bit_range[i][0]
            e = bit_range[i][1]

            if i == len(bit_range) - 1:
                s_next = bit_range[0][0]
                new_part = [s, e, self.modulus_bitlength - e + s_next]
            else:
                s_next = bit_range[i+1][0]
                new_part = [s, e, s_next - e]
                
            new_bit_range.append(new_part)

        return new_bit_range

    def stageI_findInitialRange(self, estimates):
        """
            Uses SciPy to identify peaks over a threshold and merges peaks into bit ranges.
            @input: estimates as mapping from bit position to Pr[b_i=1]/Pr[b_i=0]
            @output: [START, END, e]
        """

        def getAvgHeight(estimates):
            if len(estimates) == 0:
                return 0

            avg = 0
            for e in estimates:
                if e > 0:
                    avg += math.log(e, 10)
                else:
                    avg += 0
            return int(avg/len(estimates))

        p_wlen = 512    # area evaluated for prominence
        p_distance = 50 # distance to next peak, bigger values make result imprecise

        heuThres = getAvgHeight(estimates)

        findPeaksLookAhead = modulus_bitlength / (10 * self.hw_omega)

        peaks,_ = sc.find_peaks(estimates, height=10**heuThres, distance=p_distance, width = 1, wlen=p_wlen)
        bit_range = self.completePartsFromFindPeaks(peaks, findPeaksLookAhead)

        # build graph
        if self.dummy:
            demoGraph = pG.demonstratorPlots(self.modulus_bitlength, self.hw_omega)
            demoGraph.xstart = 96000
            demoGraph.xend = 164000
            demoGraph.ystart = 1.e-70
            demoGraph.yend = 1.e93
            demoGraph.pltEstimatesPeaksThres(estimates, peaks, heuThres, "a", show=False)
            self.dummy = False

        return bit_range

    def stageII_getTighterRanges(self, estimates, bit_range):
        """
            Tighten bit ranges relative to their height.
            @input: estimates, [START, END, e]
            @output: tighter [START, END, e]
        """
        tight_bitRange = []

        for i in range(len(bit_range)):
            start = bit_range[i][0]
            end = bit_range[i][1]
            empty = bit_range[i][2]

            """
            Do not tighten, if range is 
                good AND LEFT range is good
                --> no added value
            """
            if end - start <= empty:
                if bit_range[i-1][2] >= (bit_range[i-1][1] - bit_range[i-1][0]):
                    tight_bitRange.append(bit_range[i])
                    continue

            """
                Super small part, just ignore
            """
            if end - start < 64:
                tight_bitRange.append(bit_range[i])
                continue

            """
                Slightly tighten relative to height, increase threshold by 1/32 seems to work well.
                Height threshold is calculated seperately for least and most significant bit,
                relative height = local max / local min
            """

            # get local maxima in bit range: [RELATIVE pos, height]
            lookahead = 256 # arbitrary power of 2
            local_max, _ = peakdetect(estimates[start:end], lookahead=lookahead)

            # if no local maxima detected, take maxima of bit range
            if len(local_max) == 0:
                pos_m = 0
                h_m = 0
                for i in range(end-start):
                    if estimates[start + i] > h_m:
                        pos_m = i
                        h_m = estimates[start + i]

                # same format peakdetect
                local_max.append([i, h_m])
            

            # Heuristic value. Smaller values give better results
            lsb_ratio = 32
            msb_ratio = 32

            # LEFT
            # calculate new threshold
            h_lsb_lmax = math.log(local_max[0][1], 10)
            h_lsb = math.log(estimates[start], 10)
            h_lsb_thres = h_lsb + (h_lsb_lmax - h_lsb) * (1.0 / lsb_ratio)

            # move start until threshold reached
            pos_newStart = start
            for i in range(end-start):
                pos = start + i
                if math.log(estimates[pos], 10) >= h_lsb_thres:
                    pos_newStart = pos
                    break

            # RIGHT
            # calculate new threshold
            h_msb_lmax = math.log(local_max[-1][1], 10)
            h_msb = math.log(estimates[end], 10)
            h_msb_thres = h_msb + (h_msb_lmax - h_msb) * (1.0 / msb_ratio)

            # move end until threshold reached
            pos_newEnd = end
            for i in range(end-start, 0, -1):
                pos = start + i
                #print "est: " + str(math.log(estimates[pos], 10))
                if math.log(estimates[pos], 10) >= h_msb_thres:
                    pos_newEnd = pos
                    break

            # define new bit range
            bit_range_tight = [int(pos_newStart), int(pos_newEnd), empty + (end - pos_newEnd)]
            tight_bitRange.append(bit_range_tight)
        
        # bit_range should be naturally sorted
        tight_bitRange = self.recomputeEmptySpace(tight_bitRange)
        return tight_bitRange

    def stageIII_advantageOfEmptySpace(self, bit_range):
        """
            Merge bad parts with (succeeding) good part, if trailing empty space is large enough to form a new good part.
            @input: [START, END, e]
            @output: [START, END, e] (merged)
        """

        """
            Check for bit ranges that can be merged.
            Remember index of first and last bit range, and number of parts that are enclosed.
        """

        # [index first, index last, number of merged bad]
        lst_merging = []

        for i in range(len(bit_range)):
            start = bit_range[i][0]
            end = bit_range[i][1]
            empty = bit_range[i][2]

            # skip good parts
            if end - start <= empty:
                continue

            num_merged_bad_parts = 1
            num_merged_parts = 1
            # iterate over all other bit ranges
            for j in range(len(bit_range) - 1):
                last = (i + j) % len(bit_range)

                # good part (chance of large empty space)
                if bit_range[last][1] - bit_range[last][0] <= bit_range[last][2]:

                    # good part has large empty space

                    if last > i:
                        size_merged = bit_range[last][1] - start
                    else:
                        size_merged = self.modulus_bitlength - bit_range[last][1] + start

                    if size_merged <= bit_range[last][2]:
                        lst_merging.append([i, last, num_merged_bad_parts])
                else:
                    num_merged_parts += 1
                    if bit_range[last][1] - bit_range[last][0] > bit_range[last][2]:
                        num_merged_bad_parts += 1

        """
            Avoid subsummed merging.
        """

        # Sort for largest number of subsumed parts.
        lst_merging.sort(key=lambda x: x[2])

        # Removes subsummed merged part sets
        lst_merging_new = []
        for m in lst_merging:
            isSubsumed = False

            # if not contained in any other, take it
            for n in lst_merging:
                if m[0] > n[0] and m[0] < n[1]:
                    isSubsumed = True

            if not isSubsumed:
                lst_merging_new.append(m)

        # Creates set of non-merged prats
        bit_range_new = []
        for r in bit_range:

            isSubmerged = False
            for m in lst_merging_new:

                i_first = m[0]
                i_last = m[1]
                if r[0] >= bit_range[i_first][0] and r[1] <= bit_range[i_last][1]:
                    isSubmerged = True
                    break

            if not isSubmerged:
                bit_range_new.append(r)

        # Add merged parts
        for m in lst_merging_new:

            i_start = m[0]
            i_end = m[1]

            r_new = [bit_range[i_start][0], bit_range[i_end][1], bit_range[i_end][2]]
            bit_range_new.append(r_new)

        return bit_range_new

    def getPartitioning(self, lst_est_a, lst_est_b):
        bit_range_a = self.stageI_findInitialRange(lst_est_a)
        tight_bit_range_a = self.stageII_getTighterRanges(lst_est_a, bit_range_a)
        merged_bit_range_a = self.stageIII_advantageOfEmptySpace(tight_bit_range_a)

        # -----------------------

        bit_range_b = self.stageI_findInitialRange(lst_est_b)
        tight_bit_range_b = self.stageII_getTighterRanges(lst_est_b, bit_range_b)
        merged_bit_range_b = self.stageIII_advantageOfEmptySpace(tight_bit_range_b)

        return merged_bit_range_a, merged_bit_range_b

  
"""
    NONE PARTITIONING
"""

def readFailures(fn_failures):
    """
     Read G,a,b and reconstruct public key.
    """
    with open(fn_failures, 'r') as f:
        line = f.readline() # skip first line
        G = int(f.readline())
        a = int(f.readline())
        b = int(f.readline())
        #H = (a * G + b) % (2**self.modulus_bitlength - 1)
        return a,b,G

def checkSecrets(lst_sec_a, lst_sec_b, a, b, hw_omega):
    """
        Compare secrets with numpy secrets file
    """
    a_hw = 0
    a_rec = 0
    for i in range(len(lst_sec_a)):
        if lst_sec_a[i] == 1:
            a_rec += 2**i 
            a_hw += 1

    b_hw = 0
    b_rec = 0
    for i in range(len(lst_sec_b)):
        if lst_sec_b[i] == 1:
            b_rec += 2**i 
            b_hw += 1
    
    if a_hw != b_hw or a_hw != hw_omega:
        print ("INPUT FILES INCONSISTENT")
        print ("Found HW a: " + str(a_hw) + ", HW b: " + str(b_hw) + ", but should be " + str(hw_omega))
        exit(0)

    if a_rec != a or b_rec != b:
        print ("INPUT FILES INCONSISTENT")
        print ("Could not recreate secrets from files.")
        exit(0)

def successProbP(B, e, ep):
    """
        Return success probability for a single part.
    """

    opt_ep = min(B + B, ep)

    exp_correct = 0.0
    for i in range(1,B):
        exp_correct += min(e + B - i, i + opt_ep)

    # multiply by probability of event
    exp_correct = exp_correct * 1.0 / B

    return exp_correct / (B + opt_ep)

def successProbR(good, bad):
    """
        Return success probability accumulated over all bit ranges
    """
    sum_pr_success = 0
    sum_w = 0
    sum_w_e = 0
    sum_w_ep = 0

    for b in bad:
        pr_success = successProbP(b[1]-b[0], b[2], b[3])
        sum_w += b[1] - b[0]
        sum_w_e += b[2]
        sum_w_ep += min(b[2], b[3])
        sum_pr_success += pr_success

    return sum_pr_success/len(bad), sum_w/len(bad), sum_w_e/len(bad), sum_w_ep/len(bad)


"""
    Setup argument parser}
"""
parser = argparse.ArgumentParser("Extract partitions from estimates.")
parser.add_argument("-e", "--estimates", nargs=2, dest="fn_ests", help="Filenames estimates: estimate-a.npy estimate-b.npy.", required=True)
parser.add_argument("-s", "--secrets", nargs=2, dest="fn_sec", help="Filenames secrets: secret-a.npy secret-b.npy.")
parser.add_argument("-f", "--failures", dest="fn_fails", help="Filename failures (contains a,b,G): failures.txt")
parser.add_argument("-p", "--params", nargs=2, dest="params", help="Modulus bit-length and Hamming weight, optional default is n=756839, w = 128")
parser.add_argument("-m", "--missing", action='store_true', default=False, help="Display positions of unidentified bits.")

args = parser.parse_args()

"""
    Initialize paramaters
"""
modulus_bitlength = 756839
hw_omega = 128

if args.params:
    modulus_bitlength = args.params[0]
    hw_omega = args.params[1]

print ("--------------------------------------------")
print (" Demonstator for Failure Attack on Ramstake")
print ("--------------------------------------------")
print (" Parameters")
print (" n=" + str(modulus_bitlength) + ", w=" + str(hw_omega))

"""
    Read estimates
"""
# NOTE, SWAPPED
lst_est_b = list(np.load(args.fn_ests[0]))
lst_est_a = list(np.load(args.fn_ests[1]))

"""
    Read secrets from file generated by estimates.py and compare
"""
lst_sec_a = []
lst_sec_b = []
if args.fn_sec:
    lst_sec_a = list(np.load(args.fn_sec[0]))
    lst_sec_b = list(np.load(args.fn_sec[1]))
else:
    print ("No secrets provided.")

"""
    Read a,b,G from failures file
"""
if args.fn_fails and args.fn_sec:
    a, b, G = readFailures(args.fn_fails)
    checkSecrets(lst_sec_a, lst_sec_b, a, b, hw_omega);
print ("--------------------------------------------")


def main():
    
    print ("Partitioning...")
    """
        Partitioning
    """
    Part = Partitioning(modulus_bitlength, hw_omega)
    
    # [START, END, ENMPTY DISTANCE]
    bit_range_a, bit_range_b = Part.getPartitioning(lst_est_a, lst_est_b)

    # [START, END, EMPTY DISTANCE, PRECEDING EMPTY]
    good_a, bad_a = Part.classifyBitRangeP(bit_range_a)
    good_b, bad_b = Part.classifyBitRangeP(bit_range_b)

    """
        TABLE
    """
    def formatedOutput(bit_range_good, bit_range_bad, bit_range, lst_sec, sid, offs = "   ", showFalse = False):
        num_bits_inRange, _ = demoPartPCheck.getCoveredBits(bit_range, lst_sec, showFalse = showFalse, id_range = sid)
        num_bits_inGood, num_good_affected = demoPartPCheck.getCoveredBits(bit_range_good, lst_sec, showFalse = False, id_range = sid)
        num_bits_sampleSpace = demoPartPCheck.getNumSampleSpace(bit_range)
        
        print (offs + "# secret positions missing: " + str(hw_omega - num_bits_inRange))
        print (offs + "|I_correct|: " + str(len(bit_range_good)) + ", # used intervals: " + str(num_good_affected))
        print (offs + "|I_sample|: " + str(len(bit_range_bad)))
        print (offs + "# secret positions in I_correct: " + str(num_bits_inGood) + " of " + str(hw_omega))

    offs = "   "
    showNonidentifiedPositions = args.missing

    demoPartPCheck = cP.demonstratorPartitionChecker()
    num_bits_inRange_a, _ = demoPartPCheck.getCoveredBits(bit_range_a, lst_sec_a)
    num_bits_inRange_b, _ = demoPartPCheck.getCoveredBits(bit_range_b, lst_sec_b)

    num_bits_inGood_a, num_good_affected = demoPartPCheck.getCoveredBits(good_a, lst_sec_a, showFalse = False, id_range = "a")
    num_bits_inGood_b, num_good_affected = demoPartPCheck.getCoveredBits(good_b, lst_sec_b, showFalse = False, id_range = "b")

    num_secretBitsMissing = (2 * hw_omega) - (num_bits_inRange_a + num_bits_inRange_b)

    if args.fn_sec:
        if num_secretBitsMissing == 0:
            print (offs + "SUCCESS")
            print (offs + "Estimation allows successful attack!")
        else:
            print (offs + "Estimation was not able to identify " + str(num_secretBitsMissing) + " bit positions of the secret!")

        print ("--------------------------------------------")
        print (offs + "|I_correct|: " + str(len(good_a) + len(good_b)))
        print (offs + "|I_sample|: " + str(len(bad_a) + len(bad_b)))
        print (offs + "# secret positions in I_correct: " + str(num_bits_inGood_a + num_bits_inGood_b))
        print (offs + "# secret positions missing: " + str(2 * hw_omega - (num_bits_inGood_a + num_bits_inGood_b)))
        #print offs + "Expected number of quantum steps: exp(" + str((2 * hw_omega - (num_bits_inGood_a + num_bits_inGood_b)) / 2) + ")"
        print ("--------------------------------------------")
        #print "  A"
        #formatedOutput(good_a, bad_a, bit_range_a, lst_sec_a, sid="a", showFalse=num_secretBitsMissing)
        #print "==="
        #print "  B"
        #formatedOutput(good_b, bad_b, bit_range_b, lst_sec_b, sid="b", showFalse=num_secretBitsMissing)
        #print "--------------------------------------------"

    """
        Average ratio of empty space of bad parts
    """
    print ("Partitioning Results ")

    pr_success_a, avg_w_a, avg_e_a, avg_ep_a = successProbR(good_a, bad_a)
    pr_success_b, avg_w_b, avg_e_b, avg_ep_b = successProbR(good_b, bad_b)

    print (offs + "Avg size of B: " + str((avg_w_a + avg_w_b)/2))
    print (offs + "Avg size of e: " + str((avg_e_a + avg_e_b)/2))
    print (offs + "Avg size of e_p: " + str((avg_ep_a + avg_ep_b)/2))
    print (offs + "Avg Pr[success]: " + str((pr_success_a + pr_success_b)/2))

    all_missing_bits = (2 * hw_omega - (num_bits_inGood_a + num_bits_inGood_b))
    all_pr_success = (pr_success_a + pr_success_b)/2.0

    q_steps = math.log(1/all_pr_success**all_missing_bits, 2)
    print (offs + "E[# quantum steps] (exp): " + str(q_steps/2))
    print ("--------------------------------------------")

    print ("Plotting graphs...")
    demoGraph = pG.demonstratorPlots(modulus_bitlength, hw_omega)

    # Figure from paper
    demoGraph.ystart = 1.e-70
    demoGraph.yend = 1.e93

    demoGraph.xstart = 96000
    demoGraph.xend = 164000
    demoGraph.pltEstimates(lst_est_a, lst_sec_a, "a", show=False)
    demoGraph.pltEstimates(lst_est_b, lst_sec_b, "b", show=False)

    demoGraph.xstart = 96000
    demoGraph.xend = 200000

    demoGraph.ystart = 1.e-70
    demoGraph.yend = 1.e106
    demoGraph.pltPartitioning(lst_est_a, lst_sec_a, good_a, bad_a, "a", "", show=False, showtxt=False)
    demoGraph.pltPartitioning(lst_est_b, lst_sec_b, good_b, bad_b, "b", "", show=False, showtxt=False)

    print ("Done.")

if __name__ == "__main__":
    main()

 