#! /usr/bin/env python

'''
Andrew Till
Winter 2015

Read and plot reaction rates from PDT output
'''

#STDLIB
import os

#TPL
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

#MINE
import plotutil as putil
import PDTXS as pdtxs

def rebin_flux(oldFlux, oldTimeBins, newTimeBins):
    '''Rebin flux. oldFlux has been averaged by oldTimeBins. Assume piecewise constant in each bin. Return newFlux, averaged over newTimeBins.'''
    # Pad with zeros so the flux for times before oldTimeBins[0] and after oldTimeBins[-1] is 0.0
    oldFluxPadded = np.concatenate(([0.0], oldFlux, [0.0]))
    # allTimeBins is the union of all time bins
    allTimeBins = np.unique(np.concatenate((oldTimeBins, newTimeBins)))
    # oldBins contains the old time bin for each union time:
    # if time range [allTimeBins[i], allTimeBins[i+1]] is in [oldTimeBins[n], oldTimeBins[n+1]], oldBins[i] = n.
    oldBins = np.searchsorted(oldTimeBins, allTimeBins, 'right')[:-1]
    newBins = np.searchsorted(newTimeBins, allTimeBins, 'right')[:-1]
    #
    dtAll = np.diff(allTimeBins)
    allFlux = oldFluxPadded[oldBins] * dtAll
    # bincount scatters with the + reduction. The first and last bin are outside the range of newTimeBins
    newFlux = np.bincount(newBins, allFlux, len(newTimeBins)+1)[1:-1]
    dtNew = np.diff(newTimeBins)
    newFlux /= dtNew
    return newFlux

def read_edits(filePath):
    '''Reads PDT output in filePath. Returns editDict with keys of (editName, editType) and values of arrays of edit[timestep, group]. Also returns timestep times and editVolDict, which has edit volumes.'''
    # Keys should be (editName, editType)
    editDict = {}
    editVolDict = {}
    times = []
    numEdits = 0
    newFormat = False
    # Read edits in filePath
    with open(filePath, 'r') as fid:
        line = fid.readline()
        while line:
            if line.strip().startswith('PDT/STAPL Version'):
                # Read the PDT version, which determines what is printed
                versionPDT = int(line.split()[-1].split('M')[0])
                if versionPDT >= 1056:
                    newFormat = True
            if line.strip().startswith('<ngroups>'):
                # Read the number of groups
                numGroups = int(line.split('>')[1].split('<')[0])
            elif line.strip().startswith('<edit_reg.id>'):
                # Assumes <edit_reg.id> precedes the <edit_reg.edit>(s)
                editName = line.split('>')[1].split('<')[0]
                line = fid.readline()
                # This will fail if the edits in the XML are not ALL printed
                while not line.strip().startswith('</regions-edit_region>'):
                    # May not cover all input variations for edits
                    if newFormat and line.strip() == '<edit_reg.edit>':
                        editMaterial = None
                        line = fid.readline().strip()
                        while not line.startswith('</edit_reg.edit>'):
                            if line.startswith('<edit_reg.edit.type>'):
                                editType = line.split('>')[1].split('<')[0]
                            elif line.startswith('<edit_reg.edit.material>'):
                                editMaterial = line.split('>')[1].split('<')[0]
                            line = fid.readline().strip()
                        editNameFull = editName
                        if editMaterial:
                            editNameFull = '{0}_{1}'.format(editName, editMaterial)
                        editDict[(editNameFull, editType)] = []
                        numEdits += 1
                    elif line.strip().startswith('<edit_reg.edit>'):
                        editType = line.split('>')[1].split('<')[0]
                        editDict[(editName, editType)] = []
                        numEdits += 1
                    line = fid.readline()
            elif line.strip().find('*** Timestep') != -1:
                t = line.split()
                timestep = int(t[4].strip(','))
                if timestep == 0:
                    startTime = float(t[8])
                    times.append(startTime)
                endTime = float(t[18])
                times.append(endTime)
                # An edit region occurs once per timestep
            elif line.lstrip('[0]').strip().startswith('Edit Region'):
                # Skip line of '--------'
                line = fid.readline()
                if newFormat:
                    # New versions of PDT added another blank line
                    line = fid.readline()
                for iEdit in range(numEdits):
                    oneTimestepEdit = np.zeros(numGroups)
                    offset = 0
                    if newFormat:
                        offset += 1
                    for group in range(numGroups):
                        t = fid.readline().split()
                        editName = t[0+offset]
                        editType = t[1+offset]
                        editValue = float(t[3+offset])
                        oneTimestepEdit[group] = editValue
                    editDict[(editName, editType)].append(oneTimestepEdit)
                    # Skip line of group sum (recompute below)
                    line = fid.readline()
                    # Read edit volume
                    if newFormat:
                        t = fid.readline().split()
                        editVolume = float(t[3+offset])
                        editVolDict[(editName, editType)] = editVolume
            elif line.strip().find('keff =') != -1:
                t = line.split()
                editDict['k'] = float(t[3][:-1])
            line = fid.readline()
    # Convert editDict from a list of lists to a 2D array [time,group].
    numTimesteps = len(editDict.values()[0])
    newEditDict = {}
    for edit in editDict:
        if edit == 'k':
            newEditDict[edit] = editDict[edit]
            continue
        newEditDict[edit] = np.zeros((numTimesteps, numGroups))
        out = newEditDict[edit]
        for timestep in range(numTimesteps):
            out[timestep, :] = editDict[edit][timestep]
    # Make times an array
    times = np.array(times)
    # Return 2D-array version and times
    return newEditDict, editVolDict, numGroups, times
