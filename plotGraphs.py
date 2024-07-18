# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import os

class demonstratorPlots:
    """
        Construct plots.
    """
    def __init__(self, modulus_bitlength, hw_omega, plot_dir='results/'):

        self.plot_dir = plot_dir

        os.system('mkdir -p ' + plot_dir)
        os.system('touch ' + plot_dir + '/PLACEHOLDER')

        self.modulus_bitlength = int(modulus_bitlength)
        self.hw_omega = int(hw_omega)

        self.xstart = 0 #0 #250000
        self.xend = self.modulus_bitlength #756839 #50000 #756839 #270000 #756839 #300000

        # legend
        self.fontsize = 10

        # axis
        self.font = {
            'size': 10
        }

        self.ystart = 1.e-70
        self.yend = 1.e130

    def getPartLowerHalfOfGood(self, good, bad):
        """
            Computes the space spanning the lower half of each part
            @return: dict_bitRange_lowerHalf    - bit range --> size lower half
        """

        bit_range = good + bad
        bit_range.sort(key=lambda tup: tup[0])

        dict_pStart_lowerHalf = {}

        for i in range(len(bit_range)):

            if bit_range[i] in bad:
                continue

            pStart = bit_range[i][0]
            pEnd = bit_range[i][1]

            lower_half = 0

            # Last part has warp around offset
            if i == (len(bit_range) - 1):
                lower_half = int((self.modulus_bitlength - pStart + bit_range[0][0])/2)
            else:
                lower_half = int(pEnd - pStart)
            
            dict_pStart_lowerHalf[pStart] = lower_half

        return dict_pStart_lowerHalf

    def pltSetup(self, estimates, sid):
        plt.clf()
        plt.xlabel(u'bit position of the secret', fontdict=self.font)

        ylabel = "P[" + str(sid) + "$_1$ = 1]/P[" + str(sid) + "$_1$ = 0]"
        plt.ylabel(ylabel, fontdict=self.font)

        plt.xlim([self.xstart,self.xend])
        plt.ylim([self.ystart, self.yend])

        plt.semilogy(estimates, alpha=1, label=r'estimates', color='#87cefa', linewidth=.8, linestyle='-')

    def pltFitting(self, loc='upper right', bbox_anc=(1,1)):
        plt.tight_layout()
        plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)
        plt.legend(loc=loc, bbox_to_anchor=bbox_anc, prop={'size': self.fontsize})

    def pltPrepareBitRange(self, estimates, lst_range):
        """
            Prepares a list of values for plotting with the estimates
        """
        plt_range = [0] * self.modulus_bitlength

        for r in lst_range:
            height = max(estimates[r[0]], estimates[r[1]])
            plt_range[r[0] - r[3]] = height
            plt_range[r[1]] = height

        return plt_range

    def pltPrepareGood(self, estimates, good):
        """
            Prepares a list of values for plotting with the estimates
        """
        plt_range = [0] * self.modulus_bitlength

        for r in good:
            plt_range[r[0]] = estimates[r[0]]

        return plt_range

    # ===================================================================
    # ===================================================================

    def pltEstimates(self, estimates, secrets, sid, show=True):
        """
            Plots estimate and a threshold line.
            @input: estimates   estimates
                    thres       exponent of threshold
        """
        self.pltSetup(estimates, sid)

        # SECRETS
        for i in range(len(secrets)):
            if secrets[i] == 1:
                secrets[i] = estimates[i]

        plt.semilogy(secrets, alpha=1, label=r'secrets ' + str(sid), color='black', linewidth=1.2, linestyle=':')

        self.pltFitting()
        filename= self.plot_dir + 'estimates-' + str(sid) + '.png'
        plt.savefig(filename, bbox_inches='tight', dpi=500)
        print(f" (plotGraphs) pltEstimates wrote file {filename}")

        # vertical threshold line
        if show:
            plt.show()

       
    def pltEstimatesThres(self, estimates, secrets, thres, sid):
        """
            Plots estimate and a threshold line.
            @input: estimates   estimates
                    thres       exponent of threshold
        """
        self.pltSetup(estimates, sid)

        # SECRETS
        for i in range(len(secrets)):
            if secrets[i] == 1:
                secrets[i] = estimates[i]

        plt.semilogy(secrets, alpha=1, label=r'peaks ' + str(sid), color='black', linewidth=1.2, linestyle=':')

        # vertical threshold line
        plt.axhline(y=10**thres, color='r', linestyle='-')
        self.pltFitting()
        plt.show()

    def pltEstimatesPeaksThres(self, estimates, peaks, thres, sid, show=True, stage='peakThres'):
        """
            Plots estimate and a threshold line.
            @input: estimates   estimates
                    thres       exponent of threshold
        """
        self.pltSetup(estimates, sid)

        # SECRETS
        plt_peaks = [0] * self.modulus_bitlength
        for p in peaks:
            plt_peaks[p] = estimates[p]

        plt.semilogy(plt_peaks, alpha=1, label=r'peaks ' + str(sid), color='gray', linewidth=1, linestyle='-')

        # vertical threshold line
        plt.axhline(y=10**thres, color='r', linestyle='-', label=r'threshold T')

        self.pltFitting()
        filename = self.plot_dir + str(stage) + '-' + str(sid) + '.png'
        plt.savefig(filename, bbox_inches='tight', dpi=500)

        print(f" (plotGraphs) pltEstimatesPeaksThres wrote file {filename}")
        # vertical threshold line
        if show:
            plt.show()


    def pltClassification(self, estimates, bit_range, sid, fname = "", show=True):
        """
            Plots estimates and bit range with classification

        """
        self.pltSetup(estimates, sid)

        """
            Prepare bit_range for plotting
        """
        plt_bit_range = self.pltPrepareBitRange(estimates, bit_range)

        plt.semilogy(plt_bit_range, alpha=1, label=r'bit range', color='green', linewidth=1.5)

        for p in range(len(bit_range)):

            start = bit_range[p][0]
            end = bit_range[p][1]
            empty = bit_range[p][2]

            if start < self.xstart or end > self.xend:
                continue

            h_text = max(estimates[start], estimates[end]) * 10**5

             # ADD TEXT WITH BEGINNING POS OF PART
            plt.text(start, h_text, "b" + "$_{" + str(p) + ",s}$")
            plt.text(end, h_text, "b" + "$_{" + str(p) + ",e}$")

            if end - start <= empty:
                plt.text(start + (end-start)/2, h_text, "$correct$")
            else:
                plt.text(start + (end-start)/2, h_text, "$sample$", color='red')

        self.pltFitting()
        filename = self.plot_dir + str(fname) + 'classification-' + str(sid) + '.png'
        plt.savefig(filename , bbox_inches='tight', dpi=500)

        print(f" (plotGraphs) pltClassification wrote file {filename}")
        
        if show:
            plt.show()

    def pltPartitioning(self, estimates, secrets, good, bad, sid, fname = "", show=True, showtxt=False):

        self.pltSetup(estimates, sid)

        # SECRETS
        for i in range(len(secrets)):
            if secrets[i] == 1:
                secrets[i] = estimates[i]

        plt.semilogy(secrets, alpha=1, label=r'secrets ' + str(sid), color='black', linewidth=1.2, linestyle=':')
        

        # GOOD
        plt_good = self.pltPrepareGood(estimates, good)
        plt.semilogy(plt_good, alpha=1, label=r'start of correct part', color='green', linewidth=1.5)
            

        # ADD CORRECT PART RANGE: [p_s, max_offset]
        dict_good_lowerHalf = self.getPartLowerHalfOfGood(good, bad)

        label_added = False
        for k,v in dict_good_lowerHalf.items():

            if not label_added:
                label_added = True
                plt.plot([k, k+v], [estimates[k], estimates[k]], color = 'red', linewidth=1.5, linestyle='--', label=r'lower half')
            else:
                plt.plot([k, k+v], [estimates[k], estimates[k]], color = 'red', linewidth=1.5, linestyle='--')

        # ADD SAMPLE RANGES
        # [START, END, EMPTY DISTANCE, PRECEDING EMPTY]
        plt_bad = self.pltPrepareBitRange(estimates, bad)

        plt.semilogy(plt_bad, alpha=1, label=r'sampling range', color='magenta', linewidth=.5)

        for r in range(len(bad)):
            start = bad[r][0]
            end = bad[r][1]
            empty = bad[r][2]

            triple=[start, end, empty]

            if start < self.xstart or end > self.xend:
                continue

            r_sample = end - start
            height = max(estimates[start], estimates[end])

            if showtxt:
                plt.text(start, height, "p" + "$_{" + str(r) + ",s}$")
                plt.text(end, height, "p" + "$_{" + str(r) + ",e}$")

        if showtxt:
            for p in good:
                plt.text(p[0], estimates[p[0]], "s: " + str(p[0]) + ", e: " + str(p[1]) + ", | |: " + str(p[1]-p[0]) + ", ety: " + str(p[2]))

            for k,v in dict_good_lowerHalf.items():
                plt.text(k+v, estimates[k], "offs: " + str(v))

            for r in range(len(bad)):
                plt.text(start, height, "s: " + str(start))
                plt.text(end, height, "e: " + str(end))
                plt.text(start + (end-start)/2, height, "||: " + str(r_sample) + ", epty: " + str(empty))

        # ADD RANGE WIDTH, EMPTY WIDTH
        self.pltFitting(bbox_anc=(0.9,1))
        filename = self.plot_dir + str(fname) + 'partitioning-' + str(sid) + '.png'
        plt.savefig(filename,  bbox_inches='tight', dpi=500)

        print(f" (plotGraphs) pltPartitioning wrote file {filename}")
        if show:
            plt.show()