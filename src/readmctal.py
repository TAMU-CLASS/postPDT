#! /usr/bin/env python

'''
Andrew Till
Spring 2015

Read an MCTAL file.
Not terribly general at the moment

'''

#STDLIB
import os

#TPL
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

#MINE
import plotutil as putil

mpl.rcParams.update({'font.size': 20, 'lines.linewidth': 2})
mpl.rcParams.update({'legend.fontsize': 16})

def define_defaults():
    '''Specify default parameters'''
    # Main parameters (plotOutput is always by default False)
    verbosity = True

    # What to do
    workOpt = 'single'
    #workOpt = 'compare'

    # Base filename for the mctal file
    basename = 'pin'

    # Comparison mctal file base name
    comparename = 'pin'

    return {'verbosity': verbosity, 'workopt': workOpt, 'basename': basename, 'comparename': comparename}

###########################################################################
def do_all(inputDict):
    # Read inputs
    verbosity = inputDict['verbosity']
    workOpt = inputDict['workopt']
    inDirr = inputDict['inputdir']
    basename = inputDict['basename']
    comparename = inputDict['comparename']
    figDirr = inputDict['outputdir']

    if workOpt == 'single':
        # Read tally file
        filename = '{0}.mt'.format(basename)
        filePath = os.path.join(inDirr, filename)
        talliesDict = read_mctal(filePath)

        # Plot contents of tally file
        talliesNameDict = get_tally_names()
        plot_tallies(figDirr, basename, talliesDict, talliesNameDict, verbosity)

    elif workOpt == 'compare':
        # Read first tally file
        basename1 = basename
        filename1 = '{0}.mt'.format(basename1)
        filePath1 = os.path.join(inDirr, filename1)
        firstDict = read_mctal(filePath1)

        # Read second tally file
        basename2 = comparename
        filename2 = '{0}.mt'.format(basename2)
        filePath2 = os.path.join(inDirr, filename2)
        secondDict = read_mctal(filePath2)

        #When comparing between MCNP runs with different total volumes, need to renormalize
        if firstDict['problemStr'] and secondDict['problemStr']:
            volumeNorm1 = get_volumes(firstDict['problemStr'])['fuel']
            volumeNorm2 = get_volumes(secondDict['problemStr'])['fuel']
            renormFactor = volumeNorm2 / volumeNorm1
            for tally in secondDict:
                if tally in ['keff', 'problemStr']:
                    continue
                secondDict[tally]['avgTallies'] *= renormFactor
                secondDict[tally]['intTallies'] *= renormFactor

        # Determine which tallies to compare
        talliesToCompare = set(firstDict.keys()) & set(secondDict.keys())
        if 'keff' in talliesToCompare:
            talliesToCompare.remove('keff')
        if 'problemStr' in talliesToCompare:
            talliesToCompare.remove('problemStr')
        talliesNameDict = get_tally_names()
        talliesColorDict = get_tally_colors()

        plot_diff(figDirr, basename1, basename2, firstDict, secondDict, talliesToCompare, talliesColorDict, talliesNameDict, verbosity)

###########################################################################
def plot_diff(figDirr, basename1, basename2, firstDict, secondDict, talliesToCompare, talliesColorDict, talliesNameDict, verbosity):
    # Plot tallies by themselves
    if verbosity > 3:
        plot_tallies(figDirr, basename1, firstDict, talliesNameDict)
        plot_tallies(figDirr, basename2, secondDict, talliesNameDict)

    # Plot rxnrate rel diff
    for tally in talliesToCompare:
        plt.figure(1)
        plt.clf()
        v1 = firstDict[tally]
        v2 = secondDict[tally]
        e = v1['energyBinEdges']
        # Change first bin to 1E-3 for prettier plotting
        e[0] = 1E-3
        y = v1['avgTallies'] / v2['avgTallies'] - 1.0
        # Apply 95% confidence
        uncer = v1['avgTallies'] / v2['avgTallies'] * np.sqrt(
            np.square(v1['relStds']) + np.square(v2['relStds']) )
        uncer *= 1.959964
        # Change to percent error
        y *= 100
        uncer *= 100
        x,y = putil.get_stairs(e, y)
        x,uncer = putil.get_stairs(e, uncer)
        color = talliesColorDict.get(tally, 'b')
        label = '{0}'.format(tally)
        if tally in talliesNameDict:
            label = talliesNameDict[tally]
        plt.semilogx(x, y, color=color, label=label)
        plt.fill_between(x, y-uncer, y+uncer, color=color, alpha=0.3)
        plt.axhline(0, linestyle='--', color='r')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Reaction rate difference\n (percent with 95% confidence)')
        # Make y-axis formatted as a '%'
        yAxisFormat = '%.1f%%'
        yticks = mtick.FormatStrFormatter(yAxisFormat)
        plt.gca().yaxis.set_major_formatter(yticks)
        #plt.legend(loc='best')
        plt.tight_layout()
        #
        figName = 'p_rxnrate_{0}_diff_{1}_{2}.pdf'.format(tally, basename1, basename2)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        #
        plt.xlim(2,1000)
        plt.ylim(-2.5, 3.0)
        plt.tight_layout()
        #
        figName = 'p_rxnrate_zoom_{0}_diff_{1}_{2}.pdf'.format(tally, basename1, basename2)
        figPath = os.path.join(figDirr, figName)
        if verbosity > 3:
            plt.savefig(figPath)

    # Plot rxnrate pcm diff
    for tally in talliesToCompare:
        plt.figure(1)
        plt.clf()
        v1 = firstDict[tally]
        v2 = secondDict[tally]
        if 'norm' in v2:
            norm = v2['norm']
        elif 'norm' in v1:
            norm = v1['norm']
        else:
            continue
        e = v1['energyBinEdges']
        # Change first bin to 1E-3 for prettier plotting
        e[0] = 1E-3
        norm = v2['norm']
        y = (v1['intTallies'] - v2['intTallies']) * norm * 1E5
        uncer = np.sqrt(np.square(v1['intTallies'] * v1['relStds']) +
                        np.square(v2['intTallies'] * v2['relStds']) ) * norm * 1E5
        # Apply 95% confidence
        uncer *= 1.959964
        cumulY = np.cumsum(y[::-1])[::-1]
        cumulYEnd = cumulY[0]
        cumulUncer = np.sqrt(np.cumsum(np.square(uncer[::-1])))[::-1]
        cumulUncerEnd = cumulUncer[0]
        x,y = putil.get_stairs(e, y)
        x,uncer = putil.get_stairs(e, uncer)
        x, cumulY = putil.get_stairs(e, cumulY)
        x, cumulUncer = putil.get_stairs(e, cumulUncer)
        color = talliesColorDict.get(tally, 'b')
        label = '{0}'.format(tally)
        if tally in talliesNameDict:
            label = talliesNameDict[tally]
        # Plot error
        plt.semilogx(x, y, color=color, label=label)
        plt.fill_between(x, y-uncer, y+uncer, color=color, alpha=0.3)
        # Plot cumulative error
        #plt.semilogx(x, cumulY, color=color, linestyle='--', label=r'Cumul. sum ({0:.0f}$\pm${1:.0f})'.format(cumulYEnd, cumulUncerEnd))
        #plt.fill_between(x, cumulY-cumulUncer, cumulY+cumulUncer, color=color, alpha=0.1)
        plt.axhline(0, linestyle='--', color='r')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Reaction rate difference\n (pcm with 95% confidence)')
        #plt.legend(loc='best')
        plt.tight_layout()
        #
        figName = 'p_pcm_rxnrate_{0}_diff_{1}_{2}.pdf'.format(tally, basename1, basename2)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        #
        plt.xlim(2,1000)
        plt.ylim(-100, 100)
        plt.tight_layout()
        #
        figName = 'p_pcm_rxnrate_zoom_{0}_diff_{1}_{2}.pdf'.format(tally, basename1, basename2)
        figPath = os.path.join(figDirr, figName)
        if verbosity > 3:
            plt.savefig(figPath)

###########################################################################
def plot_tallies(figDirr, basename, talliesDict, talliesNameDict, verbosity):
    # Plot flux
    plt.figure(1)
    plt.clf()
    talliesToPlot = [14, 54, 34]
    talliesToPlot = [tally for tally in talliesToPlot if tally in talliesDict]
    for tally in talliesToPlot:
        tallyDict = talliesDict[tally]
        x = tallyDict['energyBinEdges']
        # Change first bin to 1E-3 for prettier plotting
        x[0] = 1E-3
        y = tallyDict['avgTallies']
        x,y = putil.get_stairs(x, y)
        label = '{0}'.format(tally)
        if tally in talliesNameDict:
            label = talliesNameDict[tally]
        plt.loglog(x, y, label=label)
    plt.xlabel('Energy (eV)')
    plt.ylabel('Flux (arb./eV)')
    plt.legend(loc='best')
    #
    figName = 'p_flux_{0}.pdf'.format(basename)
    figPath = os.path.join(figDirr, figName)
    plt.savefig(figPath)
    #
    #
    plt.xlim(2,1000)
    plt.ylim(1E-8,1E-4)
    #
    figName = 'p_flux_zoom_{0}.pdf'.format(basename)
    figPath = os.path.join(figDirr, figName)
    plt.savefig(figPath)

    if 51 in talliesDict and 54 in talliesDict:
        # Plot J+, J-, etc.
        Jplus = talliesDict[51]['avgTallies'][1,:]
        Jminus = talliesDict[51]['avgTallies'][0,:]
        volume = 1E2 * 1E2 * 0.39218 # Hard-coded
        phi = talliesDict[54]['avgTallies'] * volume
        e = talliesDict[54]['energyBinEdges']
        # Change first bin to 1E-3 for prettier plotting
        e[0] = 1E-3
        #
        #
        plt.figure(1)
        plt.clf()
        x,y = putil.get_stairs(e, Jplus)
        plt.loglog(x, y, label=r'$J^+$')
        x,y = putil.get_stairs(e, Jminus)
        plt.loglog(x, y, label=r'$J^-$')
        x,y = putil.get_stairs(e, phi)
        plt.loglog(x, y, label=r'$\phi$')
        plt.xlabel('Energy (eV)')
        plt.ylabel('Current or Flux (arb./eV)')
        plt.legend(loc='best')
        #
        figName = 'p_current_{0}.pdf'.format(basename)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        #
        plt.xlim(2,1000)
        plt.ylim(volume*1E-8, volume*1E-4)
        #
        figName = 'p_current_zoom_{0}.pdf'.format(basename)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        #
        plt.figure(1)
        plt.clf()
        SigmaE = Jplus / phi
        avgSigmaE = np.median(SigmaE)
        x,y = putil.get_stairs(e, SigmaE)
        plt.loglog(x, y, label='Energy-binned')
        plt.axhline(avgSigmaE, color='r', linestyle='--', label='{0:.7f} 1/cm'.format(avgSigmaE))
        plt.xlabel('Energy (eV)')
        plt.ylabel('Escape cross section (1/cm)')
        plt.legend(loc='best')
        #
        figName = 'p_esc_{0}.pdf'.format(basename)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        #
        plt.xlim(2,1000)
        plt.ylim(0.1, 3.0)
        #
        figName = 'p_esc_zoom_{0}.pdf'.format(basename)
        figPath = os.path.join(figDirr, figName)
        plt.savefig(figPath)
        #
        if verbosity:
            print 'Average escape cross section is {0} 1/cm'.format(avgSigmaE)
            print 'Commensurate radius is {0} cm'.format(1/(2.*avgSigmaE))

    # Plot rxnrate
    plt.figure(1)
    plt.clf()
    #talliesToPlot = [24, 44, 64, 74, 84, 104, 124]
    talliesToPlot = [24, 44, 64]
    talliesToPlot = [tally for tally in talliesToPlot if tally in talliesDict]
    for tally in talliesToPlot:
        tallyDict = talliesDict[tally]
        x = tallyDict['energyBinEdges']
        # Change first bin to 1E-3 for prettier plotting
        x[0] = 1E-3
        y = tallyDict['avgTallies']
        x,y = putil.get_stairs(x, y)
        label = '{0}'.format(tally)
        if tally in talliesNameDict:
            label = talliesNameDict[tally]
        plt.loglog(x, y, label=label)
    plt.xlabel('Energy (eV)')
    plt.ylabel('Reaction rate (arb./eV)')
    plt.legend(loc='best')
    #
    figName = 'p_rxnrate_{0}.pdf'.format(basename)
    figPath = os.path.join(figDirr, figName)
    plt.savefig(figPath)
    #
    #
    plt.xlim(2,1000)
    plt.ylim(1E-9,1E-3)
    #
    figName = 'p_rxnrate_zoom_{0}.pdf'.format(basename)
    figPath = os.path.join(figDirr, figName)
    plt.savefig(figPath)

###########################################################################
def read_mctal(filePath):
    '''Reads an MCNP mctal file and returns the result in a dict indexed by tally number (or 'keff' for eigenvalue). Energy bin boundaries start at 0 eV'''
    with open(filePath, 'r') as fid:
        line = fid.readline()
        mcnpVersion = int(line.split()[1])

        # Read in tally numbers
        while not line.startswith('ntal'):
            line = fid.readline()
        numTallies = int(line.split()[1])
        tallyCount = 0
        tallies = set()
        while tallyCount < numTallies:
            line = fid.readline()
            tallies |= set([int(tally) for tally in line.split()])
            tallyCount = len(tallies)

        # Read tallies
        talliesDict = {}
        while line and set(talliesDict.keys()) != tallies:
            while line.split()[0] != 'tally':
                line = fid.readline()
            tally = int(line.split()[1])
            tallyType = tally % 10
            if tallyType == 1:
                numSkip = 7
            elif tallyType == 4:
                numSkip = 8
            if mcnpVersion == 6:
                # Untested
                numSkip += 1
            for i in range(numSkip):
                line = fid.readline()
            if tallyType == 1:
                # Read number of cosine bins
                numCosineBins = int(line.split()[1]) # Assume no total
                cosineBinEdges = [-1.0]
                numLines = get_num_lines(numCosineBins, valsPerLine=6)
                for i in range(numLines):
                    cosineBinEdges += [float(val) for val in fid.readline().strip().split()]
                cosineBinEdges = np.array(cosineBinEdges)
                line = fid.readline()
            # Read number of energy bins
            numEnergyBins = int(line.split()[1]) - 1 # Do not include the 'total' bin
            energyBinEdges = [0.0]
            numLines = get_num_lines(numEnergyBins, valsPerLine=6)
            for i in range(numLines):
                energyBinEdges += [float(val) for val in fid.readline().strip().split()]
            # Make energy bins in eV
            energyBinEdges = np.array(energyBinEdges) * 1E6
            # Read tally
            for i in range(2):
                line = fid.readline()
            if tallyType == 1:
                numVals = (numEnergyBins + 1) * numCosineBins # Energy has a total bin
            if tallyType == 4:
                numVals = numEnergyBins + 1 # Include the total bin
            tallyVals = []
            numLines = get_num_lines(numVals, valsPerLine=4)
            for i in range(numLines):
                tallyVals += [float(val) for val in fid.readline().strip().split()]
            if tallyType == 1:
                # Reshape into 2D array indexed by [cosine, energy]
                # Energies move faster, cosines move slower for MCNP as well
                tempMeans = np.array(tallyVals[::2]).reshape((numCosineBins, numEnergyBins+1))
                tempSigmas = np.array(tallyVals[1::2]).reshape((numCosineBins, numEnergyBins+1))
                # Don't include the total
                tallyMeans = tempMeans[:,:-1].copy()
                tallySigmas = tempSigmas[:,:-1].copy()
                # Calculate the average in energy (currently is integrated)
                energyWidths = np.diff(energyBinEdges)
                tallyAvgMeans = tallyMeans / energyWidths[np.newaxis,:]
            elif tallyType == 4:
                # Don't include the total
                tallyMeans = np.array(tallyVals[:-2:2])
                tallySigmas = np.array(tallyVals[1:-2:2])
                # Calculate the average in energy (currently is integrated)
                energyWidths = np.diff(energyBinEdges)
                tallyAvgMeans = tallyMeans / energyWidths
            # Bundle everything into a dict
            tallyDict = {
                    'energyBinEdges':energyBinEdges,
                    'avgTallies':tallyAvgMeans,
                    'intTallies':tallyMeans,
                    'relStds':tallySigmas }
            talliesDict[tally] = tallyDict
            line = fid.readline()

        # Get keff and its uncertainty
        while line and line.split()[0] != 'kcode':
            line = fid.readline()
        if not line:
            keff = -1
            keffAbsStd = 1
        else:
            numCycles = int(line.split()[1])
            if numCycles > 6500:
                numCycles = 6500
                print 'Warning! MCNP may only print 6500 cycles worth of information'
            linesPerCycle = 4
            numLinesToSkip = (numCycles - 1) * linesPerCycle
            for i in range(numLinesToSkip):
                fid.readline()
            for i in range(3):
                line = fid.readline()
            keff = float(line.split()[1])
            keffAbsStd = float(line.split()[2])
        tallyDict = {'val':keff, 'absStd':keffAbsStd}
        talliesDict['keff'] = tallyDict

    # Add comment to talliesDict, needed for later volume lookups
    talliesDict['problemStr'] = read_comment_string(filePath)

    return talliesDict

def read_comment_string(MTFilePath):
    '''Get the problem type from the comment line.
    Assumes a .mi input file exists with the same name as the .mt file.'''
    # Milan: Currently, which volumes to use relies on reading in the MCNP input file
    # and looking at the first word. Should we instead make this a user-specified parameter?
    with open(MTFilePath, 'r') as fid:
        fid.readline()
        problemStr = fid.readline().strip()
    if problemStr.startswith('TRIGA 2D fuel pin'):
        return 'triga'
    else:
        return 0

def get_num_lines(numVals, valsPerLine):
    numLines = numVals // valsPerLine
    if numVals % valsPerLine:
        numLines += 1
    return numLines

def print_tallies(talliesDict):
    for key in sorted(talliesDict):
        print '============================='
        print key
        tallyDict = talliesDict[key]
        for subKey in sorted(tallyDict):
            print '--------------------'
            print subKey
            print tallyDict[subKey]

###########################################################################
def get_tally_names(problemStr=0):
    if problemStr == 'triga':
        talliesNameDict = {
            14:'Inner fis.', 34: 'Outer fis.', 54: 'Fuel fis.',
            24:'Inner abs.', 44: 'Outer abs.', 64: 'Fuel abs.',
            51: 'Par. cur.',
            } 
    return talliesNameDict

def get_long_tally_names(problemStr=0):
    # WARNING: PROBLEM-DEPENDENT
    if problemStr == 'triga':
        talliesNameDict = {
            14:'Inner fission', 34: 'Outer fission', 54: 'Fuel fission',
            24:'Inner absorption', 44: 'Outer absorption', 64: 'Fuel absorption',
            51: 'Partial currents',
            } 
    return talliesNameDict

def get_tally_colors():
    purple = np.asarray([128,0,105])/255.
    orange = np.asarray([180,80,0])/255.
    darkRed = np.asarray([200,0,0])/255.
    darkishMagenta = np.asarray([150,90,150])/255.
    darkishPurple = np.asarray([70,0,50])/255.
    darkishBlue = np.asarray([0,0,100])/255.
    darkishGreen = np.asarray([0,100,0])/255.
    darkishTeal = np.asarray([0,100,100])/255.
    darkishOrange = np.asarray([90,40,0])/255.
    darkishRed = np.asarray([100,0,0])/255.
    darkerishRed = np.asarray([50,0,0])/255.
    grey2 = [0.2, 0.2, 0.2]
    grey3 = [0.3, 0.3, 0.3]
    grey4 = [0.4, 0.4, 0.4]
    grey5 = [0.5, 0.5, 0.5]
    grey6 = [0.6, 0.6, 0.6]
    # WARNING: PROBLEM-DEPENDENT
    talliesColorDict = {
            14: 'b',
            24: 'm',
            34: 'g',
            44: 'r',
            54: purple,
            64: orange,
            }
    return talliesColorDict

def process_mcnp_tallies_for_pincell(mcnpDict):
    '''Process the MCNP tallies for the pincell problem. Imbues tallies with volumes and normalizations.'''
    # MCNP edits should be:
    # tally      description
    #  14 -- inner fuel  fis
    #  24 -- inner fuel  abs
    #  34 -- outer fuel  fis
    #  44 -- outer fuel  abs
    #  54 -- fuel        fis
    #  64 -- fuel        abs
    #
    # Add tally volumes to the tally dicts
    problemStr = mcnpDict['problemStr']
    mcnpVolumes = get_tally_volumes(problemStr)
    for tally in mcnpDict:
        if tally not in ['keff', 'problemStr']:
            mcnpDict[tally]['vol'] = mcnpVolumes.get(tally, -1.0)
    #
    # Create the derived tallies (e.g., total abs rate) from components
    # (no derived tallies)

def get_volumes_triga():
    '''Return volumes for the true (i.e., 2D) pincell problem for one pincell (in cm^3)'''
    height = 100
    rZr = 0.3175
    rInnerFuel = 1.6970
    rOuterFuel = 1.7411
    rClad = 1.7920
    pitchX = 4.05003
    pitchY = 3.85445
    pi = np.pi
    # (surface area of fuel)
    area = height * 2 * pi * rOuterFuel
    Vzr = height * pi * rZr**2
    VinnerFuel = height * pi * (rInnerFuel**2 - rZr**2)
    VouterFuel = height * pi * (rOuterFuel**2 - rInnerFuel**2)
    Vclad = height * pi * (rClad**2 - rOuterFuel**2)
    Vfuel = height * pi * (rOuterFuel**2 - rZr**2)
    Vmod = height * (pitchX*pitchY - pi*rClad**2)
    Vall = height * pitchX * pitchY
    volDict = {
        'zr': Vzr,
        'in': VinnerFuel,
        'out': VouterFuel,
        'fuel': Vfuel,
        'clad': Vclad,
        'mod': Vmod,
        'all': Vall,
        'area': area}
    return volDict

def get_volumes(problemStr):
    '''Return dictionary of volumes'''
    # WARNING: PROBLEM-DEPENDENT
    if problemStr == 'triga':
        return get_volumes_triga()

def get_tally_volumes(problemStr):
    '''Returns volumes, indexed by tally'''
    Vs = get_volumes(problemStr)
    # WARNING: PROBLEM-DEPENDENT
    onePinVolumes = {
        # Reaction rates in the inner fuel
        14: Vs['in'],
        24: Vs['in'],
        # Reaction rates in the outer fuel
        34: Vs['out'],
        44: Vs['out'],
        # Reaction rates in the fuel
        54: Vs['fuel'],
        64: Vs['fuel'],
    }
    # WARNING: PROBLEM-DEPENDENT
    if problemStr == 'triga':
        mcnpVolumes = {}
        mcnpVolumes.update(onePinVolumes)
    # Surface area of a pin
    mcnpVolumes[51] = Vs['area']
    return mcnpVolumes

def set_tally_normalizations(mcnpDict, verbosity):
    '''Computes and stores a normalization for each tally based on the total absorption or fission rates. These normalizations are integral in both energy and space. This allows errors in different energy bins and tallies to be compared directly'''

    problemStr = mcnpDict['problemStr']
    volTot = get_volumes(problemStr)['all']
    # WARNING: PROBLEM-DEPENDENT
    fisTot = np.sum(mcnpDict[54]['intTallies'])
    absTot = np.sum(mcnpDict[64]['intTallies'])
    norms = {'abs': absTot, 'fis': fisTot}
    vols = get_tally_volumes(problemStr)
    relVols = {}
    for tally in vols:
        relVols[tally] = vols[tally] / volTot
    # WARNING: PROBLEM-DEPENDENT
    normToUse = {
            14: 'fis', 24: 'abs',
            34: 'fis', 44: 'abs',
            54: 'fis', 64: 'abs',
             }
    for tally in sorted(mcnpDict):
        if tally in normToUse:
            norm = relVols[tally] / norms[normToUse[tally]]
            mcnpDict[tally]['norm'] = norm
            if verbosity:
                print tally, normToUse[tally],
                print np.sum(mcnpDict[tally]['intTallies']) * norm
            if verbosity > 1:
                print tally, normToUse[tally],
                print np.sum(mcnpDict[tally]['intTallies']), norms[normToUse[tally]], relVols[tally], mcnpDict[tally]['norm']

def create_summed_tally(mcnpDict, inputTallies, outputTally):
    '''Create new tally that is the sum of the inputTallies. Assumes volumes add'''
    tallySize = mcnpDict[inputTallies[0]]['intTallies'].shape
    energyBinEdges = mcnpDict[inputTallies[0]]['energyBinEdges'].copy()
    outVol = 0.0
    outIntTally = np.zeros(tallySize)
    outAbsVar = np.zeros(tallySize)
    # Sum the volumes, tallies, and uncertainties for each component (assume independent)
    for tally in inputTallies:
        tallyDict = mcnpDict[tally]
        v = tallyDict['vol']
        outVol += v
        x = v * tallyDict['intTallies']
        outIntTally += x
        s = np.square(x * tallyDict['relStds'])
        outAbsVar += s
    outRelStd = np.sqrt(outAbsVar) / outIntTally
    # intTally is integrated over energy but averaged over volume (how MCNP does things)
    outIntTally /= outVol
    # avgTally is averaged over energy and volume (first bin edge changed for plotting purposes)
    outAvgTally = outIntTally / np.diff(energyBinEdges)
    # Save the results
    mcnpDict[outputTally] = {
        'energyBinEdges': energyBinEdges,
        'avgTallies': outAvgTally,
        'intTallies': outIntTally,
        'relStds': outRelStd,
        'vol': outVol}

def renorm_tally(mcnpDict, tally, renormalization):
    '''Renormalize a tally. avg and int refer to energy, not space.
    It is assumed that volume has independently been changed to be correct.'''
    mcnpDict[tally]['avgTallies'] *= renormalization
    mcnpDict[tally]['intTallies'] *= renormalization

###########################################################################
def define_input_parser():
    import argparse
    parser = argparse.ArgumentParser(description='Reader of MCNP mctal files.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #
    defaults = define_defaults()
    #
    # If nothing is specified, verbosity is False. If -v or --verbose is specified, verbosity is 1. If -v N is specified, verbosity is N.
    parser.add_argument('-v', '--verbose', dest='verbosity', nargs='?', const=1, default=defaults['verbosity'], choices=[0,1,2,3,4], type=int)
    #parser.add_argument('-p', '--plot', action='store_true', default=False)
    parser.add_argument('-w', '--workopt', help="Work option. 'single' means look at single mctal file. 'compare' means compare two mctal files.", choices=['single', 'compare'], default=defaults['workopt'])
    parser.add_argument('-i', '--inputdir', help='Input directory', default='.')
    parser.add_argument('-o', '--outputdir', help='Output directory', default='../figures')
    parser.add_argument('-I', '--inname', dest='basename', help='Input base filename that contains the MCNP mctal file.', default=defaults['basename'])
    parser.add_argument('-C', '--comparename', help='Input base filename that contains the MCNP mctal file to be compared against.', default=defaults['comparename'])
    return parser

###########################################################################
if __name__ == '__main__':
    parser = define_input_parser()
    inputDict = vars(parser.parse_args())
    if inputDict['verbosity'] > 1:
        print 'Summary of inputs:'
        print inputDict
    do_all(inputDict)
