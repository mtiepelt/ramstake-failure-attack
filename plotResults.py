# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt 
from cycler import cycler
import os

fontsize = 15
font = {
    'size': 10
}

class demonstratorPlots:

    def __init__(self):

        self.plot_dir = 'results/'
        os.system('mkdir -p ' + self.plot_dir)
        os.system('touch ' + self.plot_dir + '/PLACEHOLDER')

        # legend
        self.fontsize = 15

        # axis
        self.font = {
            'size': 14,
        }

    def pltSetup(self):
        plt.clf()

        plt.xscale('log', base=2)
        plt.yscale('log', base=2)

        plt.xlabel(u'decryption querries', fontdict=self.font)

        ylabel = "expected quantum steps"
        plt.ylabel(ylabel, fontdict=self.font)

        plt.rc('axes', prop_cycle=(cycler('color', ['r', 'g', 'b', 'y'])))


    def pltFitting(self):
        plt.legend(loc='upper right', prop={'size': self.fontsize})
        plt.tight_layout()
        plt.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.1)

    def pltResults(self, seeds):

        self.pltSetup()

        decryption_failures = [2**(9), 2**(10), 2**(11), 2**(12)]

        linestyles = [':', '-.', '--', '-']
        linestyle_cnt = 0

        if len(seeds) > len(linestyles):
            print ("More plots than different linestyles!")

        for k,v in seeds.items():
            plt.plot(decryption_failures, v, alpha=1, label=r'seed ' + str(k), linewidth=1.2, linestyle=linestyles[linestyle_cnt % len(linestyles)])
        
            linestyle_cnt += 1

        self.pltFitting()
        filename = self.plot_dir + 'results.png'
        plt.savefig(filename, bbox_inches='tight', dpi=500)
        print(f" (plotResults) pltResults wrote file {filename}")
        plt.show()




seeds = {#
    "c14553": [2**(66), 2**(46), 2**(41), 2**(40)],
    "195bc2": [2**(57), 2**(42), 2**(45), 2**(44)],
    "0d86a4": [2**(60), 2**(51), 2**(51), 2**(51)],
    "170784": [2**(68), 2**(53), 2**(42), 2**(40)],
}

seeds_pe = {#
    "c14553": [2**(74), 2**(51), 2**(46), 2**(45)], 
    "195bc2": [2**(62), 2**(46), 2**(50), 2**(48)],
    "0d86a4": [2**(65), 2**(53), 2**(52), 2**(52)], 
    "170784": [2**(70), 2**(55), 2**(42), 2**(39)], 
}

demoPlot = demonstratorPlots() 

#demoPlot.pltResults(seeds)
demoPlot.pltResults(seeds_pe)




