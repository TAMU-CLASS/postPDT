#! /usr/bin/env python

'''
Andrew Till
Spring 2015

Compare PDT to MCNP for the pin-cell tally
'''

#STDLIB
import os

#TPL
import numpy as np

#MINE
import plotutil as putil
import readmctal as mctal
import readpdtedit as rpe
import PDTXS as pdtxs

def define_defaults():
    '''Specify default parameters'''
    # Main parameters (plotOutput is always by default False)
    verbosity = False

    # What to do
    workOpt = 'condense'
    #workOpt = 'all'

    # Base filename for the mctal file
    basenameMCNP = 'pin'

    # Base filename for the PDT output file
    basenamePDT = 'out_pin'

    # What to compare
    talliesToCompare = [24,44,64]

    # Problem type for PDT
    problemStrPDT = 'triga'

    return {'verbosity': verbosity, 'workopt': workOpt, 'mcnpbasename': basenameMCNP, 'pdtbasename': basenamePDT, 'pdtproblemstr': problemStrPDT, 'tallieslist': talliesToCompare}

def do_all(inputDict):
    # Get inputs
    verbosity = inputDict['verbosity']
    workOpt = inputDict['workopt']
    talliesToCompare = [int(s) for s in inputDict['tallieslist']]
    basenameMCNP = inputDict['mcnpbasename']
    basenamePDT = inputDict['pdtbasename']
    mcnpDirr = inputDict['mcnpdir']
    pdtDirr = inputDict['pdtdir']
    figDirr = inputDict['outputdir']
    printTable = inputDict['table']
    useLaTeX = inputDict['latex']
    problemStrPDT = inputDict['pdtproblemstr']

    filenameMCNP = '{0}.mt'.format(basenameMCNP)
    filePathMCNP = os.path.join(mcnpDirr, filenameMCNP)
    mcnpDict = mctal.read_mctal(filePathMCNP)
    problemStrMCNP = mcnpDict['problemStr']
    mctal.process_mcnp_tallies_for_pincell(mcnpDict)

    filenamePDT = '{0}.txt'.format(basenamePDT)
    pdtEdits, volDict, Eg, coarseBdrs = read_pdt_rxn_rates(pdtDirr, filenamePDT)
    pdtDict = process_pdt_edits_for_pincell(pdtEdits, volDict, Eg, problemStrPDT, problemStrMCNP, verbosity)

    # Compare tallies
    talliesToCompare = [t for t in talliesToCompare if (t in mcnpDict and t in pdtDict)]
    if verbosity:
        print 'Comparing tallies', talliesToCompare
    talliesColorDict = mctal.get_tally_colors()
    talliesNameDict = mctal.get_tally_names(problemStrMCNP)
    talliesLongNameDict = mctal.get_long_tally_names(problemStrMCNP)

    # Compare k-eff, if present
    numGroups = len(Eg) - 1
    numRRR = numGroups - 136 # Hard-coded to SHEM-361!
    if 'k' in pdtDict and 'keff' in mcnpDict:
        kPDT = pdtDict['k']
        kMCNP = mcnpDict['keff']['val']
        kStdMCNP = mcnpDict['keff']['absStd']
        kErr = (kPDT / kMCNP - 1) * 1E5
        kStdErrPlus = (kPDT / (kMCNP + kStdMCNP) - 1) * 1E5
        kStdErrMinus = (kPDT / (kMCNP - kStdMCNP) - 1) * 1E5
        kStdErrLow = min(kStdErrMinus, kStdErrPlus)
        kStdErrHigh = max(kStdErrMinus, kStdErrPlus)
        kStdErrAvg = (kPDT / kMCNP) * (kStdMCNP / kMCNP) * 1E5
        efficiency = 1E5 / (kErr * numRRR)
        if verbosity:
            print 'k (PDT) {0:.5f}'.format(kPDT)
            print 'k (MCNP) {0:.5f} +/- {1:.5f}'.format(kMCNP, kStdMCNP)
            print 'Error (pcm) {0:.0f} -- {1:.0f}'.format(kStdErrLow, kStdErrHigh)
        print 'Error (pcm) {0:.0f} +/- {1:.0f}'.format(kErr, kStdErrAvg)
        print 'Efficiency {0:.1f}'.format(efficiency)

    if workOpt == 'all':
        # Plot difference before condensing
        mctal.plot_diff(figDirr, basenameMCNP, basenamePDT, pdtDict, mcnpDict, talliesToCompare, talliesColorDict, talliesNameDict, verbosity)

    if verbosity > 1:
        N = 10
        for tally in talliesToCompare:
            print 'MCNP tally {0}'.format(tally)
            rateVals = mcnpDict[tally]['avgTallies'][:N]
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr
            #
            rateVals = mcnpDict[tally]['relStds'][:N]
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr
            #
            print 'PDT "tally" {0}'.format(tally)
            rateVals = pdtDict[tally]['avgTallies'][:N]
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr

    # Change the lower energy bound to be larger for plotting
    coarseBinWidths = np.diff(coarseBdrs[::-1])
    coarseBdrsForPlotting = coarseBdrs.copy()[::-1]
    coarseBdrsForPlotting[0] = 1E-3

    # Condense tallies
    mcnpCondensedDict = {}
    pdtCondensedDict = {}
    volDict = mctal.get_tally_volumes(problemStrMCNP)
    mctal.set_tally_normalizations(mcnpDict, verbosity)
    for tally in talliesToCompare:
        # Save volume and normalization for later
        vol = volDict[tally]
        norm = 0
        if 'norm' in mcnpDict[tally]:
            norm = mcnpDict[tally]['norm']
        # Condense MCNP energy-averaged tally
        avgTally = mcnpDict[tally]['avgTallies'][::-1]
        mcnpBdrs = mcnpDict[tally]['energyBinEdges'][::-1]
        avgCondensedTally = rebin_flux(avgTally, mcnpBdrs, coarseBdrs)[::-1]
        # Condense the uncertainty by condensing the absolute uncertainties
        absStd = mcnpDict[tally]['relStds'][::-1] * mcnpDict[tally]['avgTallies'][::-1]
        relCondensedTally = rebin_flux(absStd, mcnpBdrs, coarseBdrs)[::-1] / avgCondensedTally
        intCondensedTally = avgCondensedTally * coarseBinWidths
        mcnpCondensedDict[tally] = {'avgTallies': avgCondensedTally,
                'intTallies': intCondensedTally,
                'energyBinEdges': coarseBdrsForPlotting,
                'relStds': relCondensedTally,
                'vol': vol,
                'norm': norm}
        # Condense PDT tally
        avgTally = pdtDict[tally]['avgTallies'][::-1]
        avgCondensedTally = rebin_flux(avgTally, Eg, coarseBdrs)[::-1]
        intCondensedTally = avgCondensedTally * coarseBinWidths
        zeroTally = np.zeros(avgCondensedTally.shape)
        pdtCondensedDict[tally] = {'avgTallies': avgCondensedTally,
                'intTallies': intCondensedTally,
                'energyBinEdges': coarseBdrsForPlotting,
                'relStds': zeroTally,
                'vol': vol,
                'norm': norm}

    if verbosity > 1:
        for tally in talliesToCompare:
            print 'MCNP tally {0}'.format(tally)
            rateVals = mcnpCondensedDict[tally]['avgTallies']
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr
            #
            rateVals = mcnpCondensedDict[tally]['relStds']
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr
            #
            print 'PDT "tally" {0}'.format(tally)
            rateVals = pdtCondensedDict[tally]['avgTallies']
            rateStr = ['{0:.5e}'.format(val) for val in rateVals]
            rateStr = '\t'.join(rateStr)
            print rateStr

    # Plot differences of condensed tallies
    basenameMCNP += '_c'
    basenamePDT += '_c'
    if not printTable:
        mctal.plot_diff(figDirr, basenameMCNP, basenamePDT, pdtCondensedDict, mcnpCondensedDict, talliesToCompare, talliesColorDict, talliesNameDict, verbosity)

    # Print differences of condensed tallies
    formatStr = 'xls'
    if useLaTeX:
        formatStr = 'tex'
    if printTable:
        print_diff(pdtCondensedDict, mcnpCondensedDict, talliesToCompare, talliesLongNameDict, formatStr)

###########################################################################
def print_diff(testDict, refDict, talliesToCompare, talliesNameDict, tableType='xls'):
    '''Print a LaTeX-friendly table with comparisons to errors. Only print out uncertainty in the ref values. Go in reverse order (high to low)'''
    percent = '%'
    romanNumerals = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII']
    if tableType == 'xls':
        # For copy-pasting into Excel
        d = '\t' # delimiter
        E = '' # line end
    elif tableType == 'tex':
        d = ' & '
        E = '\\\\'
        percent = '\\%'
    if tableType == 'tex':
        print '\hline \hline'
    print 'Group' + d + 'PDT' + d + 'MCNP' + d + '' + d + '(P-M)/M' + d + 'Rxn.' + E
    print ''      + d + ''    + d + 'mean' + d + 'sigma' + d + '' + d + 'Rate' + E
    print ''      + d + '(rxn. frac.)' + d + '(rxn. frac.)'+ d + '({0})'.format(percent)   + d + '({0})'.format(percent) + d + '(pcm)' + E
    if tableType == 'tex':
        print '\hline'
    for tally in talliesToCompare:
        name = talliesNameDict[tally]
        v1 = testDict[tally]
        v2 = refDict[tally]
        # Get tally normalization to absolute values
        norm = v2['norm']
        # Get edit energy bins
        e = v1['energyBinEdges'][1:]
        # Get test value
        y1 = v1['intTallies'] * norm
        y1T = np.sum(y1)
        # Get ref values
        y2 = v2['intTallies'] * norm
        y2T = np.sum(y2)
        u2 = v2['relStds'] * 100
        u2T = np.sqrt(np.sum(np.square(u2*y2))) / y2T
        # Get relative errors and their uncertainties
        yR = (v1['avgTallies'] / v2['avgTallies'] - 1.0) * 100
        yRT = (y1T / y2T - 1.0) * 100
        #uR = v1['avgTallies'] / v2['avgTallies'] * np.sqrt(
        #        np.square(v1['relStds']) + np.square(v2['relStds']) )
        # Get pcm errors and their uncertainties
        yP = (v1['intTallies'] - v2['intTallies']) * norm * 1E5
        yPT = np.sum(yP)
        #uP = np.sqrt(np.square(v1['intTallies'] * v1['relStds']) +
        #             np.square(v2['intTallies'] * v2['relStds']) ) * norm * 1E5
        nameStr = name
        if tableType == 'tex':
            nameStr = '\\multicolumn{{6}}{{c}}{{{0}}}'.format(name)
        print nameStr + E
        numGroups = len(e)
        for g in range(numGroups-1,-1,-1):
            gs = '{0}'.format(numGroups-g-1)
            if numGroups <= len(romanNumerals):
                gs = romanNumerals[numGroups-g-1]
            y1s = '{0:.5e}'.format(y1[g])
            y2s = '{0:.5e}'.format(y2[g])
            u2s = '{0:.2f}'.format(u2[g])
            yRs = '{0:.2f}'.format(yR[g])
            yPs = '{0:.0f}'.format(yP[g])
            if tableType == 'tex':
                y1s = convert_to_latex(y1s)
                y2s = convert_to_latex(y2s)
            print gs + d + y1s + d + y2s + d + u2s + d + yRs + d + yPs + E
        gs = 'Total'
        y1s = '{0:.5e}'.format(y1T)
        y2s = '{0:.5e}'.format(y2T)
        u2s = '{0:.2f}'.format(u2T)
        yRs = '{0:.2f}'.format(yRT)
        yPs = '{0:.0f}'.format(yPT)
        if tableType == 'tex':
            y1s = convert_to_latex(y1s)
            y2s = convert_to_latex(y2s)
        print gs + d + y1s + d + y2s + d + u2s + d + yRs + d + yPs + E
        if tableType == 'tex':
            print '\hline'
        else:
            print ''

def convert_to_latex(numStr):
    '''Convert a number string to LaTeX'''
    if 'e' in numStr:
        base, exponent = numStr.split('e')
        return '${0} \\times 10^{{{1}}}$'.format(base, int(exponent))
        #return '${0}$ & $10^{{{1}}}$'.format(base, int(exponent))
    else:
        return numStr

###########################################################################
def read_pdt_rxn_rates(dirr, filename):
    '''Read PDT reaction rates and group structures'''
    # Read reaction rates
    # Since this is a SS calculation, times should be empty
    filePath = os.path.join(dirr, filename)
    editDict, volDict, numGroups, times = rpe.read_edits(filePath)

    # Read in energy boundaries (descendingly sorted)
    filename = filename.replace('out', 'mesh')
    #filename = 'xs_kFUEL_{0}.data'.format(numGroups)
    filePath = os.path.join(dirr, filename)
    xsDict = pdtxs.read_PDT_xs_generally(filePath)
    Eg = xsDict.Eg

    # Read in condensed-group boundaries (descendingly sorted)
    groupBdrDirr = '../dat/energy_groups'
    # WARNING: PROBLEM-DEPENDENT
    groupBdrName = 'edits-12.txt'
    groupBdrPath = os.path.join(groupBdrDirr, groupBdrName)
    coarseBdrs = np.loadtxt(groupBdrPath, skiprows=2, usecols=[1])

    # Make the smallest energies zero for integration consistency
    Eg[-1] = 0.0
    coarseBdrs[-1] = 0.0

    return editDict, volDict, Eg, coarseBdrs

def process_pdt_edits_for_pincell(pdtDict, volDict, Eg, problemStrPDT, problemStrMCNP, verbosity):
    '''Map PDT edits to MCNP tally structures for the pincell problem. All energies should be descendingly sorted.'''

    # None means tally is not in current MCNP decks.
    # WARNING: PROBLEM-DEPENDENT
    PDTtoMCNP_triga = {
        ('inner_edit', 'rxrate_fiss'): 14,
        ('inner_edit', 'rxrate_abs'): 24,
        ('outer_edit', 'rxrate_fiss'): 34,
        ('outer_edit', 'rxrate_abs'): 44,
        ('fuel_edit', 'rxrate_fiss'): 54,
        ('fuel_edit', 'rxrate_abs'): 64,
    }

    if problemStrMCNP == 'triga':
        PDTtoMCNP = PDTtoMCNP_triga

    # In MCNP, z height is 100 cm. In PDT, it is effectively 1
    # WARNING: PROBLEM-DEPENDENT
    pdtVolumeFactor = 1/ 100.

    mcnpVolumes = mctal.get_tally_volumes(problemStrMCNP)

    # Remap from PDT-style edits to MCNP-style tallies
    remappedPDTdict = {}
    # Make energies ascendingly sorted
    Eg = Eg[::-1]
    energyBinSizes = np.diff(Eg)
    # Renormalize PDT source strength to match MCNP integrated source strength of 1.0
    pdtSourceStr = 1.0
    # PDT in 1D should be ~ the same as MCNP in 2D for the last group for tally 24 for NN5s
    isKeig = 'k' in pdtDict
    keys = [key for key in pdtDict if key != 'k']
    for (editName, editType) in keys:
        if (editName, editType) not in PDTtoMCNP:
            if verbosity:
                print (editName, editType), 'not in map from PDT to MCNP'
            continue
        tally = PDTtoMCNP[(editName, editType)]
        if PDTtoMCNP[(editName, editType)] is None:
            if verbosity:
                print (editName, editType), 'has no corresponding MCNP entry'
            continue
        vol = volDict[(editName, editType)]
        volMCNP = mcnpVolumes[tally]
        # Change from descending to ascending order; remove timestep index
        pdtEdit = pdtDict[(editName, editType)][0,::-1]
        pdtIntTally = pdtEdit * pdtVolumeFactor
        pdtAvgTally = pdtIntTally / energyBinSizes
        zeroVec = np.zeros(len(pdtIntTally))
        remappedPDTdict[tally] = {
                'energyBinEdges': Eg,
                'avgTallies': pdtAvgTally,
                'intTallies': pdtIntTally,
                'relStds': zeroVec,
                'vol': volMCNP}
    if isKeig:
        remappedPDTdict['k'] = pdtDict['k']
    remappedPDTdict['problemStr'] = problemStrMCNP
    return remappedPDTdict

###########################################################################
def rebin_flux(oldFlux, oldEnergyBins, newEnergyBins):
    '''Rebin flux or reaction rates in energy. Assume piecewise constant in each bin. Return newFlux, averaged over newEnergyBins.'''
    # rpe's rebin_flux assumes the bins are ascendingly sorted
    return rpe.rebin_flux(oldFlux[::-1], oldEnergyBins[::-1], newEnergyBins[::-1])[::-1]

###########################################################################
def define_input_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Reader of MCNP mctal files.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    defaults = define_defaults()
    #
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=defaults['verbosity'], choices=[0,1,2,3,4], type=int)
    parser.add_argument('-w', '--workopt', help="Work option. 'condense' means condense before comparing. 'all' means compare, condense, and re-compare", choices=['condense', 'all'], default=defaults['workopt'])
    parser.add_argument('-T', '--table', help='Instead of plotting, print a table', action='store_true', default=False)
    parser.add_argument('-l', '--latex', help='Print out table in LaTeX format', action='store_true', default=False)
    parser.add_argument('-t', '--tallieslist', help='List of tallies to compare between PDT and MCNP', nargs='+', default=defaults['tallieslist'])
    parser.add_argument('-m', '--mcnpbasename', help='MCNP base file name', default=defaults['mcnpbasename'])
    parser.add_argument('-p', '--pdtbasename', help='PDT base file name', default=defaults['pdtbasename'])
    parser.add_argument('-o', '--outputdir', help='Output directory', default='../figures')
    parser.add_argument('-I', '--inputdirmcnp', dest='mcnpdir', help='MCNP directory', default='.')
    parser.add_argument('-i', '--inputdirpdt', dest='pdtdir', help='PDT directory', default='.')
    parser.add_argument('-P', '--pdtproblemstr', help='String describing the PDT problem', default=defaults['pdtproblemstr'])
    return parser

###########################################################################
if __name__ == '__main__':
    parser = define_input_parser()
    inputDict = vars(parser.parse_args())
    if inputDict['verbosity'] > 2:
        print 'Summary of inputs:'
        print inputDict
    do_all(inputDict)
