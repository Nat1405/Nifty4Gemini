# MIT License

# Copyright (c) 2015, 2017 Marie Lemoine-Busserolle

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# STDLIB
# Fixes bug in graphical pyraf tasks
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
   
import getopt
import os, glob, shutil, logging, pkg_resources, json
import pexpect as p
import time, re
from pyraf import iraf
from pyraf import iraffunctions
import astropy.io.fits
import numpy

# LOCAL

# Import config parsing.
from ..configobj.configobj import ConfigObj
from ..nifsUtils import HeaderInfo, WavelengthError

# Define constants
# Paths to Nifty data.
RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
RUNTIME_DATA_PATH = pkg_resources.resource_filename('nifty', 'runtimeData/')


def run():
    """
    Merge final cubes.
    """
    # Store current working directory for later use.
    path = os.getcwd()

    # Set up iraf
    iraf.gemini()
    iraf.nifs()
    iraf.gnirs()
    iraf.gemtools()

    # Unlearn the used tasks.
    iraf.unlearn(iraf.gemini,iraf.gemtools,iraf.gnirs,iraf.nifs)

    # Prepare the package for NIFS
    iraf.nsheaders("nifs",logfile="Nifty.log")
    iraf.set(stdimage='imt2048')
    user_clobber=iraf.envget("clobber")
    iraf.reset(clobber='yes')

    # Set up the logging file.
    log = os.getcwd()+'/Nifty.log'

    logging.info('\n#################################################')
    logging.info('#                                               #')
    logging.info('#       Start the NIFS Final Cube Merging       #')
    logging.info('#                                               #')
    logging.info('#################################################\n')

    # Load reduction parameters from ./config.cfg.
    with open('./config.cfg') as config_file:
        config = ConfigObj(config_file, unrepr=True)
        # Read general pipeline config.
        manualMode = config['manualMode']
        over = config['over']
        scienceDirectoryList = config['scienceDirectoryList']
        # Read baselineCalibrationReduction specfic config.
        mergeConfig = config['mergeConfig']
        start = mergeConfig['mergeStart']
        stop = mergeConfig['mergeStop']
        mergeType = mergeConfig['mergeType']
        use_pq_offsets = mergeConfig['use_pq_offsets']
        im3dtran = mergeConfig['im3dtran']

    valindex = start
    while valindex <= stop:
        # There are three types of merging to choose from. You can:

        if valindex == 1:
            # Merge uncorrected cubes. These have the "ctfbrsn" prefix.
            mergeCubes(scienceDirectoryList, "uncorrected", mergeType, use_pq_offsets, im3dtran, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 1 - Merge Uncorrected Individual Observations - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        if valindex == 2:
            # Merge merged cubes from each observation.
            finalMergeCubes(mergeType, "uncorrected", scienceDirectoryList, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 2 - Merge Uncorrected Merged Observation Cubes - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        if valindex == 3:
            # Merge telluric corrected cubes. These have the "actfbrsn" prefix.
            mergeCubes(scienceDirectoryList, "telluricCorrected", mergeType, use_pq_offsets, im3dtran, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 3 - Merge Telluric Corrected Individual Observations - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        if valindex == 4:
            # Merge merged cubes from each observation.
            finalMergeCubes(mergeType, "telluricCorrected", scienceDirectoryList, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 4 - Merge Telluric Corrected Merged Observation Cubes - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        if valindex == 5:
            # Merge telluric corrected AND flux calibrated cubes. These have the "factfbrsn" prefix.
            mergeCubes(scienceDirectoryList, "telCorAndFluxCalibrated", mergeType, use_pq_offsets, im3dtran, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 5 - Merge Telluric Corrected and Flux Calibrated Cubes - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        if valindex == 6:
            # Merge merged cubes from each observation.
            finalMergeCubes(mergeType, "telCorAndFluxCalibrated", scienceDirectoryList, over)
            logging.info("\n##############################################################################")
            logging.info("")
            logging.info("  STEP 6 - Merge Telluric Corrected AND Flux Calibrated Cubes - COMPLETED ")
            logging.info("")
            logging.info("##############################################################################\n")

        valindex += 1

def mergeCubes(obsDirList, cubeType, mergeType, use_pq_offsets, im3dtran, over=""):
    """MERGE

    This module contains all the functions needed to merge
    the final data cubes.

    NOTE: If you wish to shift the cubes manually in QFits View
    you can combine them in this script by making sure that you
    attach the prefix "shif" to each shifted image and save them
    in the observation directory (ie. obs108). This is necessary
    for very faint objects.

    TODO(nat): I think what we want as a final product is one cube per grating,
    per object. I need to finish generalizing this so it really works with multiple
    gratings, observations and objects.

    We should loop like this: For each object, for each grating, for each observation.
    Merge all cubes in an observation. Then Merge all cubes of the same grating. Then
    do that for every object in the program.

    INPUT:
        - Reference data cubes
        - A list of paths where final data cubes are located
        - Transformed integral field spectra

    OUTPUT:
        - Merged cubes for each observation (ie. DATE_obs##(#).fits)
        - One final merged cube from entire observation program
    """

    # Store the current working directory so we can find our way back later on.
    path = os.getcwd()
    iraffunctions.chdir(path)

    # Set the default logfile for iraf tasks.
    # TODO: Set the logfile for all iraf tasks! Right now it is not logging their output because of im3dtran...
    # It seems im3dtran doesn't have a "log" parameter.
    log = "Nifty.log"

    # Set appropriate suffix for cube type.
    if cubeType == "uncorrected":
        suffix = "_uncorrected"
        unmergedDirectory = 'products_uncorrected'
    elif cubeType == "telluricCorrected":
        suffix = "_telluricCorrected"
        unmergedDirectory = 'products_telluric_corrected'
    elif cubeType == "telCorAndFluxCalibrated":
        suffix = "_telCorAndFluxCalibrated"
        unmergedDirectory = 'products_fluxcal_AND_telluric_corrected'
    else:
        logging.info("No suffix found!")
        return

    # Create some lists here.
    listsOfCubes = []        # List of lists of cubes (one list for each science observation directory).
    mergedCubes = []         # List of Merged cubes (one merged cube for each science observation directory).
    obsidlist = []           # List of science observation id s.

    # Store the merged directory
    targetDirectory = None

    # Pixel scale in arcseconds/pixel.
    pixScale = 0.05

    # TODO(nat): implement a way to read and save cubelists to textfiles. It would be nice for users to
    # be able to edit the list of cubes to merge by hand.
    # If no Merged directory exists that contains a textfile list of cubes:
    # Go to each science directory and copy cubes from there to a new directory called Merged.

    if len(obsDirList) == 0:
        logging.warning("NifsMerge called with an empty list of directories to look for cubes. Exiting.")
        os.chdir(path)
        return

    # TODO(nat): This code seems to work really well, but it could use some polishing. Feel free to refactor nicely!
    for obsDir in obsDirList:
        # Get date, obsid and obsPath by splitting each science directory name.
        # Eg: directory name is ""/Users/ncomeau/research/newer-nifty/hd165459/20160705/H/obs13", then:
        # temp1 == ('/Users/ncomeau/research/newer-nifty/hd165459/20160705/H', 'obs13')
        # temp2 == ('/Users/ncomeau/research/newer-nifty/hd165459/20160705', 'H')
        # temp3 == ('/Users/ncomeau/research/newer-nifty/hd165459', '20160705')
        # temp4 == ('/Users/ncomeau/research/newer-nifty', 'hd165459')

        if not obsDir:
            raise ValueError("nifsMerge: There was a problem with the science directory list.")

        # TODO: make this clearer.
        temp1 = os.path.split(obsDir)
        temp2 = os.path.split(temp1[0])
        temp3 = os.path.split(temp2[0])
        temp4 = os.path.split(temp3[0])
        objname = temp3[1]
        date = temp3[1]
        obsid = temp1[1]
        obsPath = temp3[0]
        targetDirectory = temp4[0]
        try:
            os.chdir(obsDir + '/'+unmergedDirectory)
        except OSError:
            raise OSError("nifsMerge: a science directory didn't exist.")

        obsidlist.append(obsPath+'/Merged'+suffix+'/'+date+'_'+obsid)

        # Create a directory called Merged and copy all the data cubes to this directory.
        if not os.path.exists(obsPath+'/Merged'+suffix+'/'):
            os.mkdir(obsPath+'/Merged'+suffix+'/')
            logging.info('\nI am creating a directory called Merged'+str(suffix))

        Merged = obsPath+'/Merged'+suffix

        if not os.path.exists(Merged+'/'+date+'_'+obsid):
            os.mkdir(Merged+'/'+date+'_'+obsid)
            logging.info('\nI am creating a directory with date and abs ID inside Merged.')

        # If a list called shiftedcubes already exists then just merge those shifted cubes and continue.
        if glob.glob("./shift*.fits"):
            if over:
                if os.path.exists('./'+obsid+'_merged.fits'):
                    os.remove('./'+obsid+'_merged.fits')
                    iraf.gemcube(input="shif*.fits[SCI]", output=obsid+'_merged', logfile = log)
            elif not os.path.exists('./'+obsid+'_merged.fits'):
                iraf.gemcube(input="shif*.fits[SCI]", output=obsid+'_merged', logfile = log)
            else:
                logging.info("Output exists and -over- not set - shifted cubes are not being merged")
            shutil.copy('./'+obsid+'_merged.fits', Merged)
            if obsDir==obsDirList[-1]:
                os.chdir(path)
                return
            else:
                continue

        # Create a list called cubes, which stores all the cubes from a particular night.
        # Store all the cubes lists in a list of lists called listsOfCubes.
        # TODO: syntax is fairly ugly; there may be a better way to do this.

        if cubeType == "uncorrected":
            cubes = glob.glob('ctfbrsnN*.fits')          # Cubes order at this point is arbitrary so we need to sort.
            if cubes:
                cubes.sort(key=lambda x: x[-8:-5])    # Sort cubes in increasing order by last three digits.
                listsOfCubes.append(cubes)
            else:
                logging.info("Warning: no uncorrected cubes found in " + str(os.getcwd()))
                logging.info("Skipping cube merging in this directory.")
                os.chdir(path)
                return
        elif cubeType == "telluricCorrected":
            cubes = glob.glob('actfbrsnN*.fits')
            if cubes:
                cubes.sort(key=lambda x: x[-8:-5])    # Sort cubes in increasing order by last three digits.
                listsOfCubes.append(cubes)
            else:
                logging.info("Warning: no telluric corrected cubes found in " + str(os.getcwd()))
                logging.info("Skipping cube merging in this directory.")
                os.chdir(path)
                return
        elif cubeType == "telCorAndFluxCalibrated":
            cubes = glob.glob('factfbrsnN*.fits')
            if cubes:
                cubes.sort(key=lambda x: x[-8:-5])    # Sort cubes in increasing order by last three digits.
                listsOfCubes.append(cubes)
            else:
                logging.info("Warning: no flux calibrated and telluric corrected cubes found in " + str(os.getcwd()))
                logging.info("Skipping cube merging in this directory.")
                os.chdir(path)
                return
        else:
            logging.info("Invalid cube type; skipping cube merge for " + str(os.getcwd()))
            os.chdir(path)
            return
        # Copy cubes to their respective data_obsid directory within Merged.
        for cube in cubes:
            shutil.copy(cube, Merged+'/'+date+'_'+obsid)

        os.chdir(Merged)

    n=0
    for cubes in listsOfCubes:

        shiftlist = []
        os.chdir(obsidlist[n])
        iraffunctions.chdir(obsidlist[n])

        if use_pq_offsets:
            # Set the zero point p and q offsets to the p and q offsets of the first cube in each list of cubes.
            header = astropy.io.fits.open(cubes[0])
            p0 = header[0].header['POFFSET']
            q0 = header[0].header['QOFFSET']
            foff = open('offsets.txt', 'w')
            foff.write('%d %d %d\n' % (0, 0, 0))
            foff.close()

        suffix = cubes[0][-8:-5]
        if im3dtran:
            if os.path.exists('transcube'+suffix+'.fits'):
                if not over:
                    logging.info('Output already exists and -over- not set - skipping im3dtran')
                if over:
                    os.remove('transcube'+suffix+'.fits')
                    iraf.im3dtran(input = cubes[0]+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 'transcube'+suffix)
            else:
                iraf.im3dtran(input = cubes[0]+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 'transcube'+suffix)
        else:
            iraf.imcopy(cubes[0]+'[SCI][*,*,*]', 'NONtranscube'+suffix+'.fits')
        shiftlist.append('cube'+suffix+'.fits')
        iraffunctions.chdir(os.getcwd())

        for i in range(len(cubes)):
            # Skip the first cube!
            if i == 0:
                if im3dtran:
                    prefix = 'transcube'
                else:
                    prefix = 'NONtranscube'
                continue
            header2 = astropy.io.fits.open(cubes[i])
            # Check to see if we are using ALTAIR. If we are, later we will invert the x offset
            # because of the different light path.
            ALTAIR = header2[0].header['AOFOLD'].strip() == 'IN'
            suffix = cubes[i][-8:-5]

            # If user wants to merge using p and q offsets, grab those from .fits headers.
            if use_pq_offsets:
                # find the p and q offsets of the other cubes in the sequence.
                xoff = header2[0].header['POFFSET']
                yoff = header2[0].header['QOFFSET']
                # calculate the difference between the zero point offsets and the offsets of the other cubes and convert that to pixels
                if ALTAIR:
                    xShift = round(-1*(xoff - p0)/pixScale)
                else:
                    xShift = round((xoff - p0)/pixScale)
                yShift = round((yoff - q0)/pixScale)
                # write all offsets to a text file (keep in mind that the x and y offsets use different pixel scales)
                foff = open('offsets.txt', 'a')
                if im3dtran:
                    # If we swap the y and lambda axis we must also write the offsets in x, lambda, y.
                    foff.write('%d %d %d\n' % (int(xShift), 0, int(yShift)))
                else:
                    # Write offsets in regular x, y, lambda.
                    foff.write('%d\t%d\t%d\n' % (xShift, yShift, 0.))
                foff.close()

            if im3dtran:
                prefix = 'transcube'
                if os.path.exists('transcube'+suffix+'.fits'):
                    if not over:
                        logging.info('Output already exists and -over- not set - skipping im3dtran')
                    if over:
                        os.remove('transcube'+suffix+'.fits')
                        iraf.im3dtran(input = cubes[i]+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 'transcube'+suffix)
                else:
                    iraf.im3dtran(input = cubes[i]+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 'transcube'+suffix)
            else:
                prefix = 'NONtranscube'
                iraf.imcopy(cubes[i]+'[SCI][*,*,*]', prefix+suffix+'.fits')
            shiftlist.append('cube'+suffix+'.fits')

        if not use_pq_offsets:
            # Before we combine make sure a suitable offsets.txt file exists.
            a = raw_input("\nPaused. Please provide a suitable offsets.txt file in ", obsidlist[n])
            while not os.path.exists('offsets.txt'):
                a = raw_input("No offsets.txt file found. Please try again.")
            logging.info('offsets.txt found successfully for', obsidlist[n])

        if os.path.exists('cube_merged.fits'):
            if over:
                os.remove('cube_merged.fits')
                iraf.imcombine(prefix+'*', output = 'cube_merged.fits',  combine = mergeType, offsets = 'offsets.txt')
            else:
                logging.info('Output already exists and -over- not set - skipping imcombine')
        else:
            iraf.imcombine(prefix+'*', output = 'cube_merged.fits',  combine = mergeType, offsets = 'offsets.txt')
        # TODO(nat): barf. This is pretty nasty... We should fix the overwrite statements here and
        # find a way to use less nesting
        if im3dtran:
            # Transpose the cube back to x, y, lambda.
            if os.path.exists('out.fits'):
                if over:
                    os.remove('out.fits')
                    iraf.im3dtran(input='cube_merged[*,-*,*]', new_x=1, new_y=3, new_z=2, output = 'out.fits')
                else:
                    logging.info('Output already exists and -over- not set - skipping final im3dtran')
            else:
                iraf.im3dtran(input='cube_merged[*,-*,*]', new_x=1, new_y=3, new_z=2, output = 'out.fits')
            iraf.fxcopy(input=cubes[0]+'[0], out.fits', output = obsidlist[n]+'_merged.fits')
        else:
            iraf.fxcopy(input=cubes[0]+'[0], cube_merged.fits', output = obsidlist[n]+'_merged.fits')
        mergedCubes.append(obsidlist[n]+'_merged.fits')

        n+=1
        os.chdir(Merged)

    os.chdir(path)
    # Save the data with JSON so we can grab it in the next step.
    mergedData = {}
    mergedData['mergedCubes'] = mergedCubes
    mergedData['Merged'] = Merged
    with open(targetDirectory+"/mergedInfo.txt", "w") as outfile:
        json.dump(mergedData, outfile)


def old_finalMergeCubes(mergeType, over):
    """
    Merge final merged cubes from all observations.
    """
    # Load data from the previous step.
    path = os.getcwd()
    try:
        with open('mergedInfo.txt') as data_file:
            mergedData = json.load(data_file)
    except IOError:
        logging.warning("No mergedInfo.txt file found in {}. Skipping merging of merged observation cubes to a single combined cube.".format(os.getcwd()))
        return
    mergedCubes = mergedData['mergedCubes']
    Merged = mergedData['Merged']

    if len(mergedCubes)>1:
        os.chdir(Merged)
        iraffunctions.chdir(Merged)
        gratlist = []
        for i in range(len(mergedCubes)):
            cubeheader = astropy.io.fits.open(mergedCubes[i])
            grat = cubeheader[0].header['GRATING']
            gratlist.append(grat)
        print("gratlist is: {}".format(gratlist))
        # TODO(nat): right now we do more final merges here than we have to. Eg, if there are three H
        # grating directories in gratlist, we will do a final merge three times. We should only be doing this once!
        # Right now it is only a problem when overwrite is turned on but it should still be fixed.
        for n in range(len(gratlist)): # For each unique grating
            # Grab the indices of the cubes associated with that grating.
            indices = [k for k, x in enumerate(gratlist) if x==gratlist[n]]
            newcubelist = []
            for ind in indices:
                newcubelist.append(mergedCubes[ind])
            print(newcubelist)
            # Do some housekeeping before the final cube merging.
            resizeAndCenterCubes(newcubelist, over)
            makeWavelengthOffsets(newcubelist, grat)
            for i in range(len(newcubelist)):
                # Build an input string containing all the cubes to combine.
                if i==0:
                    inputstring = newcubelist[i]+'[1]'
                else:
                    inputstring += ','+newcubelist[i]+'[1]'
            if os.path.exists('temp_merged'+gratlist[n][0]+'.fits'):
                if over:
                    iraf.delete('temp_merged'+gratlist[n][0]+'.fits')
                    iraf.imcombine(inputstring, output = 'temp_merged'+gratlist[n][0]+'.fits', combine = mergeType, offsets = 'waveoffsets'+grat[0]+'.txt')
                    iraf.fxcopy(input=newcubelist[0]+'[0], temp_merged'+gratlist[n][0]+'.fits', output = 'TOTAL_merged'+gratlist[0][0]+'.fits')
                else:
                    logging.info('Output exists and -over- not set - skipping final cube merge')
            else:
                iraf.imcombine(inputstring, output = 'temp_merged'+gratlist[n][0]+'.fits', combine = mergeType, offsets = 'waveoffsets'+grat[0]+'.txt')
                iraf.fxcopy(input=newcubelist[0]+'[0], temp_merged'+gratlist[n][0]+'.fits', output = 'TOTAL_merged'+gratlist[n][0]+'.fits')
    os.chdir(path)


def finalMergeCubes(mergeType, cubeType, scienceDirectoryList, over=""):
    """
    """
    mergedDirs = findMergeDirs(cubeType, scienceDirectoryList)

    for mergeDir in mergedDirs:
        for grating in ['K', 'H', 'J', 'Z']:
            try:
                cubes = findCubes(mergeDir, grating, cubeType)
            except WavelengthError:
                continue
            if len(cubes) >= 1:
                tmp_merge_path = os.path.join(mergeDir, "total_merge_products_"+grating)
                if os.path.exists(tmp_merge_path):
                    if over:
                        os.remove(tmp_merge_path)
                        os.mkdir(tmp_merge_path)
                    else:
                        logging.info('Output exists and -over- not set - skipping final cube merge in {}.'.format(mergeDir))
                        continue
                else:
                    os.mkdir(tmp_merge_path)
                for cube in cubes:
                    shutil.copy(cube, tmp_merge_path)
                mergeCubesInDir(tmp_merge_path, mergeType, over)

def findMergeDirs(cubeType, scienceDirectoryList):
    """
    Parse scienceDirectoryList to find merged directories.
    """
    mergeDirs = []
    for scienceDir in scienceDirectoryList:
        basePath = os.path.sep.join(scienceDir.split(os.path.sep)[:-3])
        if os.path.exists(os.path.join(basePath, "Merged_"+cubeType)) and os.path.join(basePath, "Merged_"+cubeType) not in mergeDirs:
            mergeDirs.append(os.path.join(basePath, "Merged_"+cubeType))
    return mergeDirs


def findCubes(mergeDir, grating, cubeType):
    if cubeType == "uncorrected":
        cubes = glob.glob(os.path.join(mergeDir, "*", "ctfbrsn*"))
    elif cubeType == "telluricCorrected":
        cubes = glob.glob(os.path.join(mergeDir, "*", "actfbrsn*"))
    elif cubeType == "telCorAndFluxCalibrated":
        cubes = glob.glob(os.path.join(mergeDir, "*", "factfbrsn*"))

    cubes = [cube for cube in cubes if HeaderInfo(cube).grat == grating]

    if len(cubes) == 0:
        return cubes

    try:
        crWav0 = HeaderInfo(cubes[0]).crWav
        assert all([abs(HeaderInfo(cube).crWav - crWav0) < 0.01 for cube in cubes])
    except AssertionError:
        logging.error("Cubes in {} didn't all have the same central wavelength! Skipping merge in this directory.".format(mergeDir))
        raise WavelengthError()

    return cubes

def mergeCubesInDir(tmp_merge_path, mergeType, over):
    path = os.getcwd()

    os.chdir(tmp_merge_path)
    iraffunctions.chdir(tmp_merge_path)

    cubes = glob.glob("*")

    cubes = makeOffsets(cubes)

    if len(cubes) == 0:
        return

    for cube in cubes:
        if os.path.exists('t'+cube):
            if over:
                os.remove('t'+cube)
                iraf.im3dtran(input = cube+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 't'+cube)
            else:
                logging.info('Output exists and -over- not set - skipping transpose of {}.'.format(cube))
        else:
            iraf.im3dtran(input = cube+'[SCI][*,*,-*]', new_x=1, new_y=3, new_z=2, output = 't'+cube)
    
    if os.path.exists('temp_merged.fits'):
        if over:
            os.remove('temp_merged.fits')
            iraf.imcombine(",".join(['t'+cube for cube in cubes]), output = 'temp_merged.fits', combine = mergeType, offsets = 'offsets.txt')
        else:
            logging.info('Output exists and -over- not set - skipping temp merge of {}.'.format(cubes))
    else:
        iraf.imcombine(",".join(['t'+cube for cube in cubes]), output = 'temp_merged.fits', combine = mergeType, offsets = 'offsets.txt')

    if os.path.exists('TOTAL_merged.fits'):
        if over:
            os.remove('TOTAL_merged.fits')
            iraf.im3dtran(input='temp_merged[*,-*,*]', new_x=1, new_y=3, new_z=2, output = 'TOTAL_merged.fits')
        else:
            logging.info('Output exists and -over- not set - skipping final merge of {}.'.format(cubes))
    else:
        iraf.im3dtran(input='temp_merged[*,-*,*]', new_x=1, new_y=3, new_z=2, output = 'TOTAL_merged.fits')
    os.chdir(path)


def makeOffsets(frames, pixScale=0.05, outfile='offsets.txt', im3dtran=True, over=False):
    """
    Inputs:
        frames: list of absolute paths to frames.

    Returns:
        frames: list of frames that were able to have their offsets created. May have length 0 or 1; if so, won't have written an offsets.txt file.
    """

    if os.path.exists(outfile) and not over:
        logging.info("{} exists and over not set; skipping creation of an offsets file.".format(outfile))
        return frames

    frames = [frame for frame in frames if checkPQHeader(frame)]
    
    if len(frames) == 0:
        logging.info("Number of usable cubes found was {}. Skipping creating of {}.".format(len(frames), outfile))
        return frames

    if len(frames) == 1:
        with open(outfile, 'w') as f:
            f.write("0 0 0\n")
        return frames

    # Make sure frames have same wavelength increment
    wdelt0 = astropy.io.fits.open(frames[0])[1].header['CD3_3']
    try:
        assert all([abs(astropy.io.fits.open(frame)[1].header['CD3_3'] - wdelt0) < 0.001 for frame in frames])
    except AssertionError:
        logging.warning("Cubes do not all have the same wavelength increment. Skipping creation of {}. Cubes:".format(outfile))
        for frame in frames:
            logging.warning("{}".format(frame))
        return []

    frame_headers = [HeaderInfo(frame) for frame in frames]

    try:
        assert len(frame_headers) > 1
        assert len(frame_headers) == len(frames)
    except AssertionError:
        logging.warning("There was a problem with cube headers. Skipping creating of {}.".format(outfile))
        return []

    p_offs = [float(HeaderInfo(frame).poff) for frame in frames]
    q_offs = [float(HeaderInfo(frame).qoff) for frame in frames]

    min_p = min(p_offs)
    max_p = max(p_offs)
    width_p = max_p - min_p

    min_q = min(q_offs)
    max_q = max(q_offs)
    width_q = max_q - min_q

    cubeheader = astropy.io.fits.open(frames[0])
    wstart0 = cubeheader[1].header['CRVAL3']
    wdelt0 = cubeheader[1].header['CD3_3']

    for i in range(len(frame_headers)):
        if frame_headers[i].ALTAIR:
            xShift = round((width_p - float(frame_headers[i].poff))/pixScale)
        else:
            xShift = round(-1*(width_p - float(frame_headers[i].poff))/pixScale)
        yShift = round((width_q + float(frame_headers[i].qoff))/pixScale)
        
        cubeheader = astropy.io.fits.open(frames[i])
        wstart = cubeheader[1].header['CRVAL3']
        wdelt = cubeheader[1].header['CD3_3']
        waveoff = round((wstart-wstart0)/wdelt)
        # write all offsets to a text file (keep in mind that the x and y offsets use different pixel scales)
        with open(outfile, 'a') as f:
            if im3dtran:
                # If we swap the y and lambda axis we must also write the offsets in x, lambda, y.
                f.write('%d %d %d\n' % (xShift, waveoff, yShift))
            else:
                # Write offsets in regular x, y, lambda.
                f.write('%d %d %d\n' % (int(xShift), int(yShift), waveoff))

    return frames


def checkPQHeader(frame):
    try:
        headers = HeaderInfo(frame)
        poff = float(headers.poff)
        qoff = float(headers.qoff)
        cubeheader = astropy.io.fits.open(frame)
        wstart0 = float(cubeheader[1].header['CRVAL3'])
        wdelt0 = float(cubeheader[1].header['CD3_3'])
        return True
    except Exception:
        logging.info("Problem with the headers of frame {}. Excluding that frame from future processing.".format(frame))
        return False

#####################################################################################
#                                        FUNCTIONS                                  #
#####################################################################################

def resizeAndCenterCubes(cubelist, over):
    """
    Resize and center a cube based on the size of the biggest cube.
    Based on https://stackoverflow.com/a/35751427/6571817
    """
    # TODO(nat): DANGER: This does not update the relevant x and y keywords.
    # That would be a good thing to implement.

    # TODO(nat): We don't save intermediate products yet. That would be good to
    # implement.

    # Find biggest x and y cube dimensions.
    maxXSize = 0
    maxYSize = 0
    cPix1Max = None
    cPix2Max = None
    for cube in cubelist:
        cubeData = astropy.io.fits.open(cube)[1].data
        cubeHeader = astropy.io.fits.open(cube)[1].header
        xSize = cubeData.shape[2]
        ySize = cubeData.shape[1]
        if xSize > maxXSize:
            maxXSize = xSize
            cPix1Max = cubeHeader['CRPIX1']
        if ySize > maxYSize:
            maxYSize = ySize
            cPix2Max = cubeHeader['CRPIX2']

    # Resize the cubes to the biggest, overwriting the old cube. Not optimal to overwrite the old one!
    for cube in cubelist:
        # Open the cubeData and header
        oldCube = astropy.io.fits.open(cube)
        cubeData = oldCube[1].data

        oldShape = cubeData.shape
        # Make a numpy zeros array as big as the largest cube.
        newShape = (oldShape[0], maxYSize, maxXSize)
        result = numpy.zeros(newShape)

        # Calculate X and Y offsets.
        #xOffset = int((maxXSize - cubeData.shape[2]) / 2.0)
        #yOffset = int((maxYSize - cubeData.shape[1]) / 2.0)

        #result[:, yOffset:oldShape[1]+yOffset, xOffset:oldShape[2]+xOffset] = cubeData
        result = cubeData

        oldCube[1].data = result
        #oldCube[1].header['CRPIX1'] = cPix1Max
        #oldCube[1].header['CRPIX2'] = cPix2Max

        if os.path.exists("temp.fits"):
            os.remove("temp.fits")
        oldCube.writeto("temp.fits")
        os.remove(cube)
        shutil.move("temp.fits", cube)


def makeWavelengthOffsets(cubelist, grat):
    """
    Get wavelength shift of cubes.
    """
    pixScale = 0.05
    cubeheader0 = astropy.io.fits.open(cubelist[0])
    wstart0 = cubeheader0[1].header['CRVAL3']
    p0 = cubeheader0[0].header['POFFSET']
    q0 = cubeheader0[0].header['QOFFSET']
    # If a user provided a waveoffsetsGRATING.txt, skip the creation of it.
    if os.path.exists('waveoffsets{0}.txt'.format(grat[0])):
        logging.info("\nwaveoffsets file exists; skipping creation of it.")
        return
    fwave = open('waveoffsets{0}.txt'.format(grat[0]), 'w')
    fwave.write('%d %d %d\n' % (0, 0, 0))
    for i in range(len(cubelist)):
        if i == 0:
            continue
        cubeheader = astropy.io.fits.open(cubelist[i])
        # Check to see if we are using ALTAIR. If we are, later we will invert the x offset
        # because of the different light path.
        ALTAIR = cubeheader[0].header['AOFOLD'].strip() == 'IN'
        # find the p and q offsets of the other cubes in the sequence.
        xoff = cubeheader[0].header['POFFSET']
        yoff = cubeheader[0].header['QOFFSET']
        # calculate the difference between the zero point offsets and the offsets of the other cubes and convert that to pixels
        if ALTAIR:
            xShift = round(-1*(xoff - p0)/pixScale)
        else:
            xShift = round((xoff - p0)/pixScale)
        yShift = round((yoff - q0)/pixScale)

        wstart = cubeheader[1].header['CRVAL3']
        wdelt = cubeheader[1].header['CD3_3']
        waveoff = round((wstart-wstart0)/wdelt)
        fwave.write('%d %d %d\n' % (xShift, yShift, waveoff))
    fwave.close()

#---------------------------------------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    pass
