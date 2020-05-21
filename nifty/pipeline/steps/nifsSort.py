#!/usr/bin/env python
# -*- coding: utf-8 -*-

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

################################################################################
#                Import some useful Python utilities/modules                   #
################################################################################

# STDLIB
# Fixes bug in graphical pyraf tasks
from sys import platform as sys_pf
if sys_pf == 'darwin':
    import matplotlib
    matplotlib.use("TkAgg")
   
from xml.dom.minidom import parseString
import urllib
from pyraf import iraf
import astropy.io.fits
import os, sys, shutil, glob, math, logging, pkg_resources, time, datetime, re
import numpy as np
# Import config parsing.
from ..configobj.configobj import ConfigObj

# LOCAL

# Import custom Nifty functions.
# TODO(nat): goodness, this is a lot of functions. It would be nice to split this up somehow.
from ..nifsUtils import getUrlFiles, getFitsHeader, FitsKeyEntry, stripString, stripNumber, \
datefmt, checkOverCopy, checkQAPIreq, checkDate, writeList, checkEntry, timeCalc, checkSameLengthFlatLists, \
rewriteSciImageList, datefmt, downloadQueryCadc

# Import NDMapper gemini data download, by James E.H. Turner.
from ..downloadFromGeminiPublicArchive import download_query_gemini

# Define constants
# Paths to Nifty data.
RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
RUNTIME_DATA_PATH = pkg_resources.resource_filename('nifty', 'runtimeData/')
    

def start():
    """
    nifsSort

    This module contains all the functions needed to copy and sort
    the NIFS raw data, where the data is located in a local directory.

    INPUT FILES:
    + rawPath files
      - Science frames
      - Calibration frames
      - Telluric frames
      - Acquisition frames (optional, but if data is copied from archive
        then acquisition frames will be copied and sorted)

    OUTPUT:
      - Sorted data
      - Lists of paths to the calibrations, science and telluric frames
      - Names of frames stored in text files for later use by the pipeline

    If -c True or a program or date is specified with -p or -d data will be copied from
    Gemini North internal network (used ONLY within Gemini).

    Args:
        rawPath (string):   Local path to raw files directory. Specified with -q at command line.
        over (boolean): If True old files will be overwritten during data reduction. Specified
                        with -o at command line. Default: False.
        program (string): OT observation id (used only within Gemini network). Specified with
                          -p at command line. "GN-2013B-Q-109".

    TODO(nat): add a way to open .tar.bz gemini public archive data, unzip it, and verify the md5 of
               each.

    """

    # Store current working directory for later use.
    path = os.getcwd()

    # Format logging options.
    FORMAT = '%(asctime)s %(message)s'
    DATEFMT = datefmt()

    """# Set up the logging file.
    logging.basicConfig(filename='Nifty.log',format=FORMAT,datefmt=DATEFMT,level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # This lets us logging.info(to stdout AND a logfile. Cool, huh?
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)"""

    # Set up the logging file.
    log = os.getcwd()+'/Nifty.log'

    logging.info('\n####################################')
    logging.info('#                                  #')
    logging.info('#  Start NIFS sorting and copying  #')
    logging.info('#                                  #')
    logging.info('####################################\n')

    # Load reduction parameters from runtimeData/config.cfg.
    with open('./config.cfg') as config_file:
        options = ConfigObj(config_file, unrepr=True)
        # Read general config.
        over = options['over']
        manualMode = options['manualMode']
        # Read sort specific config.
        sortConfig = options['sortConfig']
        rawPath = sortConfig['rawPath']
        program = sortConfig['program']
        dataSource = sortConfig['dataSource']
        proprietaryCookie = sortConfig['proprietaryCookie']
        skyThreshold = sortConfig['skyThreshold']
        sortTellurics = sortConfig['sortTellurics']
        telluricTimeThreshold = sortConfig['telluricTimeThreshold']

    # Check for invalid command line input. Cannot both copy from Gemini and sort local files.
    # Exit if -q <path to raw frame files> and -c True are specified at command line (cannot copy from
    # Gemini North internal network AND use local raw data).
    if rawPath and program:
        logging.info("\nError in sort: both a local path and a Gemini program ID (to download from Gemini Public Archive) were provided.\n")
        raise SystemExit

    # Download data from gemini public archive to ./rawData/.
    if program:
        if not os.path.exists('./rawData'):
            os.mkdir('./rawData')
        if dataSource == 'CADC':
            logging.info('\nDownloading data from the CADC archive to ./rawData. This will take a few minutes.')
            downloadQueryCadc(program, os.getcwd()+'/rawData')
        elif dataSource == 'GSA':
            if proprietaryCookie:
                download_query_gemini(program, './rawData', proprietaryCookie)
            else:
                download_query_gemini(program, './rawData')
        else:
            raise ValueError("Invalid dataSource in config file.")
        
        rawPath = os.getcwd()+'/rawData'

    ############################################################################
    ############################################################################
    #                                                                          #
    #              CASE 1: NIFS RAW DATA in local directory                    #
    #                                                                          #
    #    These conditions are used when a local path to raw files              #
    #    is specified with -q at command line.                                 #
    #                                                                          #
    #                                                                          #
    ############################################################################
    ############################################################################


    # IF a local raw directory path is provided, sort data.
    if rawPath:
        if manualMode:
            a = raw_input("About to enter makePythonLists().")
        allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, skyFrameList, telskyFrameList, obsidDateList, sciImageList = makePythonLists(rawPath, skyThreshold)
        if manualMode:
            a = raw_input("About to enter sortScienceAndTelluric().")
        objDirList, scienceDirectoryList, telluricDirectoryList = sortScienceAndTelluric(allfilelist, skyFrameList, telskyFrameList, sciImageList, rawPath)
        if manualMode:
            a = raw_input("About to enter sortCalibrations().")
        calibrationDirectoryList = sortCalibrations(arcdarklist, arclist, flatlist, flatdarklist, ronchilist, objectDateGratingList, objDirList, obsidDateList, sciImageList, rawPath, manualMode)
        # If a telluric reduction will be performed sort the science and telluric images based on time between observations.
        # This will NOT be executed if -t False is specified at command line.
        if sortTellurics:
            if manualMode:
                a = raw_input("About to enter matchTellurics().")
            matchTellurics(telluricDirectoryList, scienceDirectoryList, telluricTimeThreshold)

    # Exit if no or incorrectly formatted input is given.
    else:
        logging.info("\nERROR in sort: Enter a program ID, observation date, or directory where the raw files are located.\n")
        raise SystemExit

    ############################################################################
    ############################################################################
    #                                                                          #
    #               CASE 2: NIFS RAW DATA in GEMINI PRIVATE NETWORT            #
    #                                                                          #
    #     These conditions are used if a program, date or copy is specified    #
    #     with -p, -d or -c True at command line. Files can be copied from     #
    #     Gemini North internal network and sorted.                            #
    #                                                                          #
    #                                                                          #
    ############################################################################
    ############################################################################

    """# TODO(nat): implement private gemini archive downloads.
    elif copy or program or date:
        try:
            import gemini_sort
        except ImportError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: gemini_sort.py module is needed to find NIFS data")
            logging.info("                      within the Gemini Internal Network. This option")
            logging.info("                      is only available when the pipeline is run at")
            logging.info("                      the observatory.")
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            raise SystemExit
        else:
            gemini_sort.start(over, copy, program, date)"""


    os.chdir(path)

    # Update runtimeData/config.cfg with the paths to each:
    # 1) Science observation directory
    # 2) Calibration observation directory
    # 3) Telluric observation directory
    logging.info("\nnifsSort: writing scienceDirectoryList, calibrationDirectoryList and telluricDirectoryList in ./config.cfg.")
    with open('./config.cfg') as config_file:
        options = ConfigObj(config_file, unrepr=True)
    options['scienceDirectoryList'] = scienceDirectoryList
    options['telluricDirectoryList'] = telluricDirectoryList
    options['calibrationDirectoryList'] = calibrationDirectoryList
    with open('./config.cfg', 'w') as config_file:
        options.write(config_file)

##################################################################################################################
#                                                                                                                #
#                                                   FUNCTIONS                                                    #
#                                                                                                                #
##################################################################################################################

def makePythonLists(rawPath, skyThreshold):

    """Creates python lists of file names by type. No directories are created and no
    files are copied in this step."""

    allfilelist = [] # List of tellurics, science frames, aquisitions... But not calibrations!
    flatlist = [] # List of lamps on flat frames.
    flatdarklist = [] # List of lamps off flat frames.
    ronchilist = [] # List of ronchi flat frames.
    arclist = [] # List of arc frames.
    arcdarklist = [] # List of arc dark frames.

    objectDateGratingList = [] # 2D list of object (science or telluric) name, date pairs.

    skyFrameList = [] # List of sky frames.
    telskyFrameList = [] # List of telluric sky frames.

    obsidDateList = [] # 2D list of date, observation id pairs.
    sciDateList = [] # List of unique dates by science (including sky) frames.

    sciImageList = [] # List of science observation directories.

    # Store current working directory for later use.
    path = os.getcwd()

    # If files were copied from Gemini Internal network raw files directory will
    # be path+"/rawPath".
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path+'/rawPath'

    # Change to raw files directory, copy and sort FILE NAMES into lists of each type (Eg: science frames, ronchi flat frames).
    # Sort by opening the .fits headers, reading header data into variables and sorting based on those variables.
    # DOES NOT COPY RAW .fits DATA IN THIS STEP
    os.chdir(rawPath)
    logging.info("\nPath to raw file directory is: " + str(rawPath))

    logging.info("\nI am making lists of each type of file.")

    # Make a list of all the files in the rawPath directory.
    rawfiles = glob.glob('N*.fits')

    # Sort and copy each filename in the rawfiles directory into lists.
    for entry in rawfiles:

        # Open the .fits header.
        header = astropy.io.fits.open(entry)

        # Store information in variables.
        instrument = header[0].header['INSTRUME']
        if instrument != 'NIFS':
            # Only grab frames belonging to NIFS raw data!
            continue
        obstype = header[0].header['OBSTYPE'].strip()
        ID = header[0].header['OBSID']
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        aper = header[0].header['APERTURE']
        # If object name isn't alphanumeric, make it alphanumeric.
        objname = header[0].header['OBJECT']
        objname = re.sub('[^a-zA-Z0-9\n\.]', '', objname)
        poff = header[0].header['POFFSET']
        qoff = header[0].header['QOFFSET']

        # Make a list of science, telluric and acquisition frames.
        # Use the copied variable (1 not copied, 0 copied) to check later that
        # the file was copied correctly. allfilelist is a 2D list of
        # [[filename1, copied, obsclass1], [filename2, copied, obsclass2]] pairs.
        if obstype == 'OBJECT' and (obsclass == 'science' or obsclass == 'acq' or obsclass == 'acqCal' or obsclass == 'partnerCal'):

            # Append a list of [filename, copied, obsclass] to the list. copied is
            # 1 if not copied, 0 if copied. obsclass is used later for checks.
            templist = [entry, 1, obsclass]
            allfilelist.append(templist)

            # Create a list of science sky frames.
            # Differentiating between on target and sky frames.
            rad = math.sqrt(((poff)**2) + ((qoff)**2))

            # If the offsets are outside a circle of 5.0 units in radius, append to skyFrameList.
            if obsclass == 'science':
                sciImageList.append(entry)
                if rad  >= skyThreshold:
                    skyFrameList.append(entry)

            # Create a list of telluric sky frames.
            if obsclass == 'partnerCal':
                if rad >= skyThreshold:
                    telskyFrameList.append(entry)

            # Create sciDateList: list of unique dates of science observations.
            if obsclass == 'science':
                # Append if list is empty or not a duplicate of last entry.
                if not sciDateList or not sciDateList[-1]==date:
                    sciDateList.append(date)

        # Add arc frame names to arclist.
        if obstype == 'ARC':
            arclist.append(entry)

        # Add arc dark frame names to arcdarklist.
        if obstype == 'DARK':
            arcdarklist.append(entry)

        # Add lamps on flat frames to flatlist,
        # add lamps off flat frames to flatdarklist,
        # add ronchi flat frames to ronchilist.
        # Lamps on and lamps off flats, and lamps on and lamps off ronchis are
        # seperated by mean number of counts per pixel. Arbitrary threshold is
        # if mean_counts < 500, it is a lamps off flat or ronchi.
        if obstype == 'FLAT':

            if aper == 'Ronchi_Screen_G5615':
                # Only use lamps on ronchi flat frames.
                # Open the image and store pixel values in an array and
                # take the mean of all pixel values.
                array = astropy.io.fits.getdata(entry)
                mean_counts = np.mean(array)

                # Once the mean is stored in mean_counts we can check whether the
                # frame is a lamps off ronchi or a lamps on ronchi based on the counts.
                # 500.0 is an arbitrary threshold that appears to work well.
                if mean_counts < 500.0:
                    pass
                else:
                    ronchilist.append(entry)

            else:
                # Open the image and store pixel values in an array and
                # take the mean of all pixel values.
                array = astropy.io.fits.getdata(entry)
                mean_counts = np.mean(array)

                # Once the mean is stored in mean_counts we can check whether the
                # frame is a sky or an object based on the counts. 500.0 is an
                # arbitrary threshold that appears to work well.
                if mean_counts < 500.0:
                    flatdarklist.append(entry)
                else:
                    flatlist.append(entry)

    # Based on science (including sky) frames, make a list of unique [object, date] list pairs to be used later.
    for i in range(len(rawfiles)):

        header = astropy.io.fits.open(rawfiles[i])
        # Store information in variables.
        instrument = header[0].header['INSTRUME']
        if instrument != 'NIFS':
            # Only grab frames belonging to NIFS raw data!
            continue
        date = header[0].header['DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj= re.sub('[^a-zA-Z0-9\n\.]', '', obj)
        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID']
        grat = header[0].header['GRATING'][0:1]

        if obsclass == 'science':
            list1 = [obj, date, grat]
            # Append if list is empty or not a duplicate of last entry.
            if not objectDateGratingList or not list1 in objectDateGratingList:
                objectDateGratingList.append(list1)

    # Make list of unique [date, obsid] pairs from FLATS. If flat was taken on the same day as a science
    # frame, append that flat date. If not, append an arbitrary unique date from sciDateList.
    # This is so we can sort calibrations later by date and observation id.
    n = 0
    for flat in flatlist:
        header = astropy.io.fits.open(flat)
        # Store information in variables.
        instrument = header[0].header['INSTRUME']
        if instrument != 'NIFS':
            # Only grab frames belonging to NIFS raw data!
            continue
        obsid = header[0].header['OBSID']
        date = header[0].header['DATE'].replace('-','')
        # Make sure no duplicate dates are being entered.
        if flatlist.index(flat)==0 or not oldobsid==obsid:
            #if date in sciDateList:
            list1 = [date, obsid]
            obsidDateList.append(list1)
            #else:
                # Ugly fix, we have to check there aren't more flats than science dates.
            #    if n < len(sciDateList):
            #        list1 = [sciDateList[n], obsid]
             #       obsidDateList.append(list1)
            #n+=1
        oldobsid = obsid

    os.chdir(path)

    # Print information for user.
    logging.info("\nTotal number of files found by type.\n")
    logging.info("Length allfilelist (science and telluric frames): " +  str(len(allfilelist)))
    logging.info("Length sciImageList (science and science sky frames): " +  str(len(sciImageList)))
    logging.info("Length arclist (arc frames): " + str(len(arclist)))
    logging.info("Length arcdarklist (arc dark frames): " + str(len(arcdarklist)))
    logging.info("Length flatlist (lamps on flat frames): " + str(len(flatlist)))
    logging.info("Length flatdarklist (lamps off flat frames): "+str(len(flatdarklist)))
    logging.info("Length ronchilist (ronchi flat frames): "+str(len(ronchilist)))
    logging.info("Length skyFrameList (science sky frames): "+str(len(skyFrameList)))
    logging.info("Length telskyFrameList (telluric sky frames): "+str(len(telskyFrameList)))

    # Store number of telluric, telluric, sky, telluric sky and acquisition frames in number_files_to_be_copied.
    number_files_to_be_copied = len(allfilelist)

    # Store number of arcs, arc darks, lamps on flats, lamps off flats and ronchi flats in number_calibration_files_to_be_copied.
    number_calibration_files_to_be_copied = len(arclist) + len(arcdarklist) +\
                         len(flatlist) + len(flatdarklist) + len(ronchilist)

    logging.info("\nTotal number of frames to be copied: " + str(number_files_to_be_copied + number_calibration_files_to_be_copied))

    return allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, skyFrameList, telskyFrameList, obsidDateList, sciImageList

#----------------------------------------------------------------------------------------#

def sortScienceAndTelluric(allfilelist, skyFrameList, telskyFrameList, sciImageList, rawPath):

    """Sorts the science frames, tellurics and acquisitions into the appropriate directories based on date, grating, obsid, obsclass.
    """

    # Store number of science, telluric, sky, telluric sky and acquisition frames in number_files_to_be_copied
    number_files_to_be_copied = len(allfilelist)

    # Initialize a counter to track how many files were copied. If this is different than
    # number_files_to_be_copied logging.info("a warning for the user at the end of sortScienceAndTelluric.")
    number_files_that_were_copied = 0

    logging.info("\n\nMaking new directories and copying files. In this step I will process " + str(number_files_to_be_copied) + " files.")

    # List of paths sorted by object and date. ['path/to/object1/date1', 'path/to/object1/date2'].
    objDirList = []
    # Create a 2D list scienceDirList. Second part is a path to a science directory.
    # First part is a list of calculated times of each frame in that science directory.
    # Eg: [[[5400,6500,7200], '/path/to/first/science'], [[3400,4300,5200], '/path/to/second/science']...]
    scienceDirList = []
    # List of paths to telluric directories.
    telDirList = []

    path = os.getcwd()
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path+'/rawPath'

    # Make new sorted directories to copy files in to.
    logging.info("\nMaking new directories.\n")

    # All data is sorted by science frame data.
    # For each science frame, create a "science_object_name/date/" directory in
    # the current working directory.
    for entry in allfilelist:
        header = astropy.io.fits.open(rawPath+'/'+entry[0])

        objname = header[0].header['OBJECT']
        objname = re.sub('[^a-zA-Z0-9\n\.]', '', objname)
        obsclass = header[0].header['OBSCLASS']
        date = header[0].header[ 'DATE'].replace('-','')

        if obsclass=='science':
            if not os.path.exists(path+'/'+objname):
                os.mkdir(path+'/'+objname)
            if not os.path.exists(path+'/'+objname+'/'+date):
                os.mkdir(path+'/'+objname+'/'+date)
                objDir = path+'/'+objname+'/'+date
                if not objDirList or not objDirList[-1]==objDir:
                    objDirList.append(objDir)
            else:
                objDir = path+'/'+objname+'/'+date
                if not objDirList or not objDirList[-1]==objDir:
                    objDirList.append(objDir)

    # For each science frame, create a "science_object_name/date/grating/observationid/"
    # directory in the current working directory.
    for entry in allfilelist:
        header = astropy.io.fits.open(rawPath+'/'+entry[0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'][-3:].replace('-','')
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)

        if obsclass=='science':
            # Important- calculate the time of day in seconds that the science (and sky) frames
            # were taken. Used because one way to match telluric frames with science frames is
            # to pair the frames closest in time.
            time = timeCalc(rawPath+'/'+entry[0])

            objDir = path+'/'+obj
            # Create a directory for each observation date (YYYYMMDD) in objDir/.
            if not os.path.exists(objDir+'/'+date):
                os.mkdir(objDir+'/'+date)
            # Create a directory for each grating used in objDir/YYYYMMDD/.
            if not os.path.exists(objDir+'/'+date+'/'+grat):
                os.mkdir(objDir+'/'+date+'/'+grat)
            # Create a directory for each obsid (eg. obs25) in objDir/YYYYMMDD/grating/.
            if not os.path.exists(objDir+'/'+date+'/'+grat+'/obs'+obsid):
                os.mkdir(objDir+'/'+date+'/'+grat+'/obs'+obsid)
                # If a new directory append time of science (or sky) frame and directory name to scienceDirList.
                scienceDirList.append([[time], objDir+'/'+date+'/'+grat+'/obs'+obsid])
            # Else if a new list or not a duplicate of the previous entry append time and directory name to scienceDirList.
            elif not scienceDirList or not scienceDirList[-1][1]==objDir+'/'+date+'/'+grat+'/obs'+obsid:
                scienceDirList.append([[time], objDir+'/'+date+'/'+grat+'/obs'+obsid])
            # IF A DUPLICATE:
            # Append the time to an existing time list.
            #
            # Eg, before appending: [[[5400,6500,7200], '/path/to/first/science'], [[5200], '/path/to/second/science']]
            # If we are processing the second frame in /path/to/second/science/, scienceDirList[-1][1] will equal /path/to/second/science.
            # Then append the new time to the second entry.
            #
            # [[[5400,6500,7200], '/path/to/first/science'], [[5200, NEWTIMEHERE], '/path/to/second/science']]

            elif scienceDirList[-1][1] == objDir+'/'+date+'/'+grat+'/obs'+obsid:
                scienceDirList[-1][0].append(time)


    # Copy science and acquisition frames to the appropriate directory.
    logging.info("\nCopying Science and Acquisitions.\nCopying science frames and science acquisitions.\nNow copying: ")

    for i in range(len(allfilelist)):
        header = astropy.io.fits.open(rawPath+'/'+allfilelist[i][0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'][-3:].replace('-','')
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)

        # Only grab the most recent aquisition frame.
        if i!=len(allfilelist)-1:
            header2 = astropy.io.fits.open(rawPath+'/'+allfilelist[i+1][0])
            obsclass2 = header2[0].header['OBSCLASS']
            obj2 = header2[0].header['OBJECT']
            obj2 = re.sub('[^a-zA-Z0-9\n\.]', '', obj2)

        # Copy sky and science frames to appropriate directories. Write two text files in
        # those directories that store the names of the science frames and sky frames for later
        # use by the pipeline.
        if obsclass=='science':
            logging.info(allfilelist[i][0])
            objDir = path+'/'+obj
            shutil.copy(rawPath+'/'+allfilelist[i][0], objDir+'/'+date+'/'+grat+'/obs'+obsid+'/')
            number_files_that_were_copied += 1
            # Update status flag to show entry was copied.
            allfilelist[i][1] = 0
            # Create an scienceFrameList in the relevant directory.
            if allfilelist[i][0] not in skyFrameList:
                writeList(allfilelist[i][0], 'scienceFrameList', objDir+'/'+date+'/'+grat+'/obs'+obsid+'/')
            # Create a skyFrameList in the relevant directory.
            if allfilelist[i][0] in skyFrameList:
                writeList(allfilelist[i][0], 'skyFrameList', objDir+'/'+date+'/'+grat+'/obs'+obsid+'/')

        # Copy the most recent acquisition in each set to a new directory to be optionally
        # used later by the user for checks (not used by the pipeline).
        if obsclass=='acq' and obsclass2=='science':
            logging.info(allfilelist[i][0])
            # create an Acquisitions directory in objDir/YYYYMMDD/grating
            if not os.path.exists(path+'/'+obj2+'/'+date+'/'+grat+'/Acquisitions/'):
                os.makedirs(path+'/'+obj2+'/'+date+'/'+grat+'/Acquisitions/')
            shutil.copy(rawPath+'/'+allfilelist[i][0], path+'/'+obj2+'/'+date+'/'+grat+'/Acquisitions/')
            number_files_that_were_copied += 1
            allfilelist[i][1] = 0

    # Copy telluric frames to the appropriate folder.
    # Note: Because the 'OBJECT' of a telluric file header is different then the
    # science target, we need to sort by date, grating AND most recent time.
    logging.info("\nCopying telluric frames.\nNow copying: ")
    for i in range(len(allfilelist)):
        header = astropy.io.fits.open(rawPath+'/'+allfilelist[i][0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'][-3:].replace('-','')
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)
        telluric_time = timeCalc(rawPath+'/'+allfilelist[i][0])


        if obsclass=='partnerCal':
            logging.info(allfilelist[i][0])
            timeList = []
            for k in range(len(scienceDirList)):
                # Make sure date and gratings match.
                tempDir = scienceDirList[k][1].split(os.sep)
                if date in tempDir and grat in tempDir:
                    # Open the times of all science frames in science_directory.
                    times = scienceDirList[k][0]
                    # Find difference in each time from the telluric frame we're trying to sort.
                    diffList = []
                    for b in range(len(times)):
                        difference = abs(telluric_time-scienceDirList[k][0][b])
                        templist = []
                        templist.append(difference)
                        templist.append(scienceDirList[k][1])
                        diffList.append(templist)
                    # Find the science frame with the smallest difference.
                    minDiff = min(diffList)
                    # Pass that time and path out of the for loop.
                    timeList.append(minDiff)
            # Out of the for loop, compare min times from different directories.
            if timeList:
                closest_time = min(timeList)
                # Copy the telluric frame to the path of that science frame.
                path_to_science_dir = closest_time[1]
                path_to_tellurics = os.path.split(path_to_science_dir)[0]

                # Create a Tellurics directory in science_object_name/YYYYMMDD/grating.

                if not os.path.exists(path_to_tellurics + '/Tellurics'):
                    os.mkdir(path_to_tellurics + '/Tellurics')
                # Create an obsid (eg. obs25) directory in the Tellurics directory.
                if not os.path.exists(path_to_tellurics+'/Tellurics/obs'+obsid):
                    os.mkdir(path_to_tellurics+'/Tellurics/obs'+obsid)
                    telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
                elif not telDirList or not telDirList[-1]==path_to_tellurics+'/Tellurics/obs'+obsid:
                    telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
                shutil.copy(rawPath+'/'+allfilelist[i][0], path_to_tellurics+'/Tellurics/obs'+obsid+'/')
                number_files_that_were_copied += 1
                allfilelist[i][1] = 0
                # Create an scienceFrameList in the relevant directory.
                if allfilelist[i][0] not in telskyFrameList:
                    writeList(allfilelist[i][0], 'tellist', path_to_tellurics+'/Tellurics/obs'+obsid+'/')
                # Create a skyFrameList in the relevant directory.
                if allfilelist[i][0] in telskyFrameList:
                    writeList(allfilelist[i][0], 'skyFrameList', path_to_tellurics+'/Tellurics/obs'+obsid+'/')

    # Modify scienceDirList to a format telSort can use.
    tempList = []
    for i in range(len(scienceDirList)):
        tempList.append(scienceDirList[i][1])
    scienceDirList = tempList

    #------------------------------ TESTS -------------------------------------#

    # Check to see which files were not copied.
    logging.info("\nChecking for non-copied science, tellurics and acquisitions.\n")
    for i in range(len(allfilelist)):
        # Check the copied flag. If not 0, logging.info("the entry.")
        if allfilelist[i][1] != 0:
            logging.info(str(allfilelist[i][0]) + " " + str(allfilelist[i][2]) + " was not copied.")
    logging.info("\nEnd non-copied science, tellurics and acquisitions.\n")

    # Check that all science frames and sky frames were copied.
    count_from_raw_files = len(sciImageList) + len(skyFrameList)

    count = 0
    for science_directory in scienceDirList:
        for file in os.listdir(science_directory):
            if file.endswith('.fits'):
                count += 1

    if count_from_raw_files != count:
        logging.info("\nWARNING: " + str(count_from_raw_files - count) + " science frames (or sky frames) \
        were not copied.\n")
    else:
        logging.info("\nExpected number of science and sky frames copied.\n")

    logging.info("\nDone sorting and copying science and tellurics. Moving on to Calibrations.\n")

    # Test for telluric directories with mis-identified sky frames. For now, we identify sky frames
    # based on absolute P and Q offsets, not relative to a zero point. This can cause problems.
    # TODO(nat): look into if it is worth it to use relative P and Q offsets.
    for telluric_directory in telDirList:
        if os.path.exists(telluric_directory + '/skyFrameList') and not os.path.exists(telluric_directory + '/tellist'):
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: a telluric directory exists that Nifty thinks ")
            logging.info("                      contains only sky frames. Nifty uses absolute")
            logging.info("                      P and Q offsets to identify sky frames; a target")
            logging.info("                      not being at 0, 0 P and Q can cause this. If this is not")
            logging.info("                      the case you can try adjusting the skyThreshold")
            logging.info("                      parameter in Nifty's configuration.")
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            os.chdir(telluric_directory)
            rewriteSciImageList(2.0, "Telluric")
            try:
                sciImageList = open('tellist', "r").readlines()
                logging.info("\nSucceeded; a telluric frame list exists in " + str(os.getcwd()))
            except IOError:
                logging.info("\nWARNING: no telluric frames found in " + str(os.getcwd()) + ". You may have to adjust the skyThreshold parameter.")
                a = raw_input("Please make a tellist (list of telluric frames) in " + str(telluric_directory))
            os.chdir(path)

    os.chdir(path)

    return objDirList, scienceDirList, telDirList

#----------------------------------------------------------------------------------------#

def sortCalibrations(arcdarklist, arclist, flatlist, flatdarklist, ronchilist, objectDateGratingList, objDirList, obsidDateList, sciImageList, rawPath, manualMode):

    """Sort calibrations into appropriate directories based on date.
    """
    calDirList = []
    filelist = ['arclist', 'arcdarklist', 'flatlist', 'ronchilist', 'flatdarklist']

    # Save path for later use. The rawPath part is for Gemini North network sorting.
    path1 = os.getcwd()
    if rawPath:
        rawPath = rawPath
    else:
        rawPath = path1+'/rawPath'

    # Set up some tests and checks.
    count = 0
    expected_count = len(arcdarklist) + len(arclist) + len(flatlist)\
          + len(flatdarklist) + len(ronchilist)

    logging.info("\nI am attempting to sort " + str(expected_count) + " files.\n")

    # To make sure data was copied later in the pipeline:
    # Add a small copied flag to each frame in calibration file lists.
    new_flatlist = []
    for i in range(len(flatlist)):
        # Transform 1D list into 2D list of [[filename, 'copied']]
        # "copied" is 1 for not copied and 0 for copied.
        templist = []
        templist.append(flatlist[i])
        templist.append(1)
        new_flatlist.append(templist)
    flatlist = new_flatlist

    new_flatdarklist = []
    for i in range(len(flatdarklist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(flatdarklist[i])
        templist.append(1)
        new_flatdarklist.append(templist)
    flatdarklist = new_flatdarklist

    new_arclist = []
    for i in range(len(arclist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(arclist[i])
        templist.append(1)
        new_arclist.append(templist)
    arclist = new_arclist

    new_arcdarklist = []
    for i in range(len(arcdarklist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(arcdarklist[i])
        templist.append(1)
        new_arcdarklist.append(templist)
    arcdarklist = new_arcdarklist

    new_ronchilist = []
    for i in range(len(ronchilist)):
        # Transform 1D list into 2D list.
        templist = []
        templist.append(ronchilist[i])
        templist.append(1)
        new_ronchilist.append(templist)
    ronchilist = new_ronchilist

    os.chdir(rawPath)

    # Create Calibrations directories in each of the observation date directories based on existence of
    # lamps on flats. Eg: YYYYMMDD/Calibrations
    # Sort lamps on flats.
    logging.info("\nSorting flats:")
    # Create a flag so we only warn about non-standard gratings once.
    grating_warning_flag = False
    for i in range(len(flatlist)):
        header = astropy.io.fits.open(flatlist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        if grating not in ["K", "J", "H", "Z"]:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: non-standard (non K, J, H, K) grating encountered. ")
            logging.info("                      NIFTY has not been tested with non-standard")
            logging.info("                      gratings!")
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
        # TODO(nat): this is horrendous. Do this in a better way.
        # "Flat is better than nested."
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        for entry in objectDateGratingList:
                            if entry[1] == date and entry[2] == grating:
                                if not os.path.exists(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating):
                                    os.mkdir(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                    calDirList.append(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                else:
                                    if path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating not in calDirList:
                                        calDirList.append(path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                # Copy lamps on flats to appropriate directory.
                                shutil.copy('./'+flatlist[i][0], path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating)
                                flatlist[i][1] = 0
                                logging.info(flatlist[i][0])
                                count += 1
                                path = path1+'/'+entry[0]+'/'+entry[1]+'/Calibrations_'+grating+'/'
                                # Create a flatlist in the relevent directory.
                                # Create a text file called flatlist to store the names of the
                                # lamps on flats for later use by the pipeline.
                                writeList(flatlist[i][0], 'flatlist', path)

    # Sort lamps off flats.
    logging.info("\nSorting lamps off flats:")
    for i in range(len(flatdarklist)):
        os.chdir(rawPath)
        header = astropy.io.fits.open(flatdarklist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+flatdarklist[i][0], objDir+'/Calibrations_'+grating+'/')
                        flatdarklist[i][1] = 0
                        logging.info(flatdarklist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # Create a flatdarklist in the relevant directory.
                        writeList(flatdarklist[i][0], 'flatdarklist', path)

    # Sort ronchi flats.
    logging.info("\nSorting ronchi flats:")
    for i in range(len(ronchilist)):
        os.chdir(rawPath)
        header = astropy.io.fits.open(ronchilist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+ronchilist[i][0], objDir+'/Calibrations_'+grating+'/')
                        ronchilist[i][1] = 0
                        logging.info(ronchilist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # create a ronchilist in the relevant directory
                        writeList(ronchilist[i][0], 'ronchilist', path)

    # Sort arcs.
    logging.info("\nSorting arcs:")
    for i in range(len(arclist)):
        header = astropy.io.fits.open(arclist[i][0])
        date = header[0].header['DATE'].replace('-','')
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            if date in objDir:
                if not os.path.exists(objDir+'/Calibrations_'+grating):
                    os.mkdir(objDir+'/Calibrations_'+grating)
                shutil.copy('./'+arclist[i][0], objDir+'/Calibrations_'+grating+'/')
                arclist[i][1] = 0
                logging.info(arclist[i][0])
                count += 1
                path = objDir+'/Calibrations_'+grating+'/'
                # Create an arclist in the relevant directory.
                writeList(arclist[i][0], 'arclist', path)

    # Sort arc darks.
    logging.info("\nSorting arc darks:")
    for i in range(len(arcdarklist)):
        header = astropy.io.fits.open(arcdarklist[i][0])
        obsid = header[0].header['OBSID']
        grating = header[0].header['GRATING'][0:1]
        for objDir in objDirList:
            for item in obsidDateList:
                if obsid in item:
                    date = item[0]
                    if date in objDir:
                        if not os.path.exists(objDir+'/Calibrations_'+grating):
                            os.mkdir(objDir+'/Calibrations_'+grating)
                        shutil.copy('./'+arcdarklist[i][0], objDir+'/Calibrations_'+grating+'/')
                        arcdarklist[i][1] = 0
                        logging.info(arcdarklist[i][0])
                        count += 1
                        path = objDir+'/Calibrations_'+grating+'/'
                        # Create an arcdarklist in the relevant directory.
                        writeList(arcdarklist[i][0], 'arcdarklist', path)

    # Check that each file in flatlist was copied.
    for i in range(len(flatlist)):
        if flatlist[i][1] == 1:
            logging.info(str(flatlist[i][0])+ " was not copied.")


    # ---------------------------- Tests ------------------------------------- #

    # Check to see how many calibrations were copied.
    if expected_count - count == 0:
        logging.info("\nI sorted the " + str(expected_count) + " expected calibrations.\n")
    else:
        logging.info("\nI did not copy " + str(expected_count - count) + " calibration file(s).\n")

    # Check each calibration file list to see which ones were not copied.
    # Check that each file in flatlist was copied.
    for i in range(len(flatlist)):
        if flatlist[i][1] == 1:
            logging.info(str(flatlist[i][0])+ " from flatlist was not copied.")

    # Check that each file in flatdarklist was copied.
    for i in range(len(flatdarklist)):
        if flatdarklist[i][1] == 1:
            logging.info(str(flatdarklist[i][0])+ " from flatdarklist was not copied.")

    # Check that each file in ronchilist was copied.
    for i in range(len(ronchilist)):
        if ronchilist[i][1] == 1:
            logging.info(str(ronchilist[i][0])+ " from ronchilist was not copied.")

    # Check that each file in arclist was copied.
    for i in range(len(arclist)):
        if arclist[i][1] == 1:
            logging.info(str(arclist[i][0])+ " from arclist was not copied.")

    # Check that each file in arcdarklist was copied.
    for i in range(len(arcdarklist)):
        if arcdarklist[i][1] == 1:
            logging.info(str(arcdarklist[i][0])+ " from arcdarklist was not copied.")


    # Change back to original working directory.
    os.chdir(path1)

    # Check that each science directory exists and has associated calibration data.
    # Pseudocode (repeated below with actual code):
    # For each science directory, make sure that:
    # a calibrations directory is present.
    # flatlist exists and has more than one file.
    # flatdarklist exists and has more than one file.
    # arclist exists and has more than one file.
    # arcdarklist exists and has more than one file.
    # ronchilist exists and has more than one file.


    # TODO(nat): This is horrifying. Good grief, wrap these repetitive calls in a function!
    logging.info("\nChecking that each science image has required calibration data. ")
    # For each science image, read its header data and try to change to the appropriate directory.
    # Check that:
    for i in range(len(sciImageList)):
        header = astropy.io.fits.open(rawPath+'/'+sciImageList[i])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'][-3:].replace('-','')
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)

        # a science and Calibrations directory are present.
        try:
            os.chdir(path1+'/'+obj+'/'+date+'/'+grat+'/obs'+obsid+'/')
            os.chdir('../../Calibrations_'+grat+'/')
        except OSError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no Calibrations directory found for ")
            logging.info("                      science frame "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            continue

        # flatlist exists and has more than one file.
        try:
            flatListFile = open('flatlist', "r").readlines()
            if len(flatListFile) <= 1:
                logging.info("\n#####################################################################")
                logging.info("#####################################################################")
                logging.info("")
                logging.info("     WARNING in sort: only 1 lamps on flat frame found for science")
                logging.info("                      frame "+str(sciImageList[i]))
                logging.info("")
                logging.info("#####################################################################")
                logging.info("#####################################################################\n")
        except IOError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no flatlist found for science frame")
            logging.info("                      "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")

            if not manualMode:
                # Sometimes flats can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some flats.
                foundflatFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through flatlist and see if there is an flat taken on this date
                for i in range(len(flatlist)):
                    header = astropy.io.fits.open(rawPath+'/'+flatlist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an flatlist.
                        shutil.copy(rawPath + '/' + flatlist[i][0], './')
                        writeList(flatlist[i][0], 'flatlist', path)
                        logging.info("\n#####################################################################")
                        logging.info("#####################################################################")
                        logging.info("")
                        logging.info("     WARNING in sort: found a flat taken one day after a science frame.")
                        logging.info("                      "+str(sciImageList[i]))
                        logging.info("                       using that.")
                        logging.info("")
                        logging.info("#####################################################################")
                        logging.info("#####################################################################\n")
                        foundflatFlag = True
                        flatlist[i][1] = 0
                if not foundflatFlag:
                    # If that quick check fails, give user a chance to try and provide an flat file.
                    a = raw_input("\n Please provide a textfile called flatlist in " + str(os.getcwd()))


        # flatdarklist exists and has more than one file.
        try:
            flatDarkListFile = open('flatdarklist', "r").readlines()
            if len(flatDarkListFile) <= 1:
                logging.info("\n#####################################################################")
                logging.info("#####################################################################")
                logging.info("")
                logging.info("     WARNING in sort: only 1 lamps off flat frame found for science")
                logging.info("                      frame "+str(sciImageList[i]))
                logging.info("")
                logging.info("#####################################################################")
                logging.info("#####################################################################\n")
        except IOError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no flatdarklist found for science frame")
            logging.info("                      "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            if not manualMode:
                # Sometimes flatdarks can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some flatdarks.
                foundflatdarkFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through flatdarklist and see if there is an flatdark taken on this date
                for i in range(len(flatdarklist)):
                    header = astropy.io.fits.open(rawPath+'/'+flatdarklist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an flatdarklist.
                        shutil.copy(rawPath + '/' + flatdarklist[i][0], './')
                        writeList(flatdarklist[i][0], 'flatdarklist', path)
                        logging.info("\n#####################################################################")
                        logging.info("#####################################################################")
                        logging.info("")
                        logging.info("     WARNING in sort: found a flatdark taken one day after a science frame.")
                        logging.info("                      "+str(sciImageList[i]))
                        logging.info("                       using that.")
                        logging.info("")
                        logging.info("#####################################################################")
                        logging.info("#####################################################################\n")
                        foundflatdarkFlag = True
                        flatdarklist[i][1] = 0
                if not foundflatdarkFlag:
                    # If that quick check fails, give user a chance to try and provide an flatdark file.
                    a = raw_input("\n Please provide a textfile called flatdarklist in " + str(os.getcwd()) + \
                    " or be sure not to attempt a wavelength calibration for this directory.")
        # Make sure flatlist and flatdarklist are the same length. nsflat() complains otherwise.
        checkSameLengthFlatLists()

        # arclist exists.
        try:
            arcListFile = open('arclist', "r").readlines()
        except IOError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no arclist found for science frame")
            logging.info("                      "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")

            if not manualMode:
                # Sometimes arcs can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some arcs.
                foundArcFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through arclist and see if there is an arc taken on this date
                for i in range(len(arclist)):
                    header = astropy.io.fits.open(rawPath+'/'+arclist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an arclist.
                        shutil.copy(rawPath + '/' + arclist[i][0], './')
                        writeList(arclist[i][0], 'arclist', path)
                        logging.info("\n#####################################################################")
                        logging.info("#####################################################################")
                        logging.info("")
                        logging.info("     WARNING in sort: found an arc taken one day after a science frame.")
                        logging.info("                      "+str(sciImageList[i]))
                        logging.info("                       using that.")
                        logging.info("")
                        logging.info("#####################################################################")
                        logging.info("#####################################################################\n")
                        foundArcFlag = True
                        arclist[i][1] = 0
                if not foundArcFlag:
                    # If that quick check fails, give user a chance to try and provide an arc file.
                    a = raw_input("\n Please provide a textfile called arclist in " + str(os.getcwd()) + \
                    " or be sure not to attempt a wavelength calibration for this directory.")

        # arcdarklist exists.
        try:
            arcDarkListFile = open('arcdarklist', "r").readlines()
        except IOError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no arcdarklist found for science frame")
            logging.info("                      "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            if not manualMode:
                # Sometimes arcdarks can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some arcdarks.
                foundarcdarkFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through arcdarklist and see if there is an arcdark taken on this date
                for i in range(len(arcdarklist)):
                    header = astropy.io.fits.open(rawPath+'/'+arcdarklist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an arcdarklist.
                        shutil.copy(rawPath + '/' + arcdarklist[i][0], './')
                        writeList(arcdarklist[i][0], 'arcdarklist', path)
                        logging.info("\n#####################################################################")
                        logging.info("#####################################################################")
                        logging.info("")
                        logging.info("     WARNING in sort: found an arcdark taken one day after a science frame.")
                        logging.info("                      "+str(sciImageList[i]))
                        logging.info("                       using that.")
                        logging.info("")
                        logging.info("#####################################################################")
                        logging.info("#####################################################################\n")
                        foundarcdarkFlag = True
                        arcdarklist[i][1] = 0
                if not foundarcdarkFlag:
                    # If that quick check fails, give user a chance to try and provide an arcdark file.
                    a = raw_input("\n Please provide a textfile called arcdarklist in " + str(os.getcwd()) + \
                    " or be sure not to attempt a wavelength calibration for this directory.")
        # ronchilist exists and has more than one file.
        try:
            ronchiListFile = open('ronchilist', "r").readlines()
            if len(ronchiListFile) <= 1:
                logging.info("\n#####################################################################")
                logging.info("#####################################################################")
                logging.info("")
                logging.info("     WARNING in sort: only 1 ronchi flat frame found for science frame")
                logging.info("                      "+str(sciImageList[i]))
                logging.info("")
                logging.info("#####################################################################")
                logging.info("#####################################################################\n")
        except IOError:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: no ronchilist found for science frame")
            logging.info("                      "+str(sciImageList[i]))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
            if not manualMode:
                # Sometimes ronchis can be taken a day after the observing night. First
                # look for these, and if they are not found, ask the user to provide some ronchis.
                foundronchiFlag = False
                # Get date after the science observation
                t=time.strptime(date,'%Y%m%d')
                newdate=datetime.date(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(1)
                # Loop through ronchilist and see if there is an ronchi taken on this date
                for i in range(len(ronchilist)):
                    header = astropy.io.fits.open(rawPath+'/'+ronchilist[i][0])
                    date = header[0].header[ 'DATE'].replace('-','')
                    if str(date) == newdate.strftime('%Y%m%d'):
                        # If so, copy it to the appropriate calibrations directory and write an ronchilist.
                        shutil.copy(rawPath + '/' + ronchilist[i][0], './')
                        writeList(ronchilist[i][0], 'ronchilist', path)
                        logging.info("\n#####################################################################")
                        logging.info("#####################################################################")
                        logging.info("")
                        logging.info("     WARNING in sort: found a ronchi taken one day after a science frame.")
                        logging.info("                      "+str(sciImageList[i]))
                        logging.info("                       using that.")
                        logging.info("")
                        logging.info("#####################################################################")
                        logging.info("#####################################################################\n")
                        foundronchiFlag = True
                        ronchilist[i][1] = 0
                if not foundronchiFlag:
                    # If that quick check fails, give user a chance to try and provide an ronchi file.
                    a = raw_input("\n Please provide a textfile called ronchilist in " + str(os.getcwd()) + \
                    " or be sure not to attempt a wavelength calibration for this directory.")
        os.chdir(path1)

    # Change back to original working directory.
    os.chdir(path1)

    # ---------------------------- End Tests --------------------------------- #

    return calDirList

#----------------------------------------------------------------------------------------#

def matchTellurics(telDirList, obsDirList, telluricTimeThreshold):

    """Matches science images with the telluric frames that are closest in time. 
    Assumes the science frames and telluric frames were observed on the same UT day.
    Creates a file in each telluric observation directory called scienceMatchedTellsList.
    scienceMatchedTellsList lists the obsid of the science images (ie. obs123) and then the
    science images with this obsid that match the telluric observation.

    EXAMPLE:    obs28
                N20130527S0264
                N20130527S0266
                obs30
                N201305727S0299

    """

    logging.info("\nI am matching science images with tellurics closest in time.\n")
    # For each science frame, find the telluric observation that is closest in time to it.

    # First, do it for each date/grating pair.

    # Get list of unique target, date, grating triples
    target_date_grates = getTargetDateGrates(obsDirList)

    basePath = os.path.sep.join(os.path.normpath(obsDirList[0]).split(os.path.sep)[:-4])
    for target, date, grating in target_date_grates:
        targetDateGratePath = os.path.join(basePath, target, date, grating)
        # Get a list of paths to the science frames by opening the scienceFrameLists
        science_frame_paths = getScienceFrames(targetDateGratePath, full_path=True)
        
        # and get a list of paths to the telluric frames by opening up the tellists
        telluric_frame_paths = getTelluricFrames(targetDateGratePath, full_path=True, listToOpen='tellist')

        # Remove old scienceMatchedTellsLists if they exist
        removeOldScienceMatchedTellsList(telluric_frame_paths)

        # For each science frame, find the telluric frame it is closest in time to.
        # Append that frame name and observation to that directories scienceMatchedTellsList.

        for science_frame in science_frame_paths:

            telluric_frame = matchScienceTelluric(science_frame, telluric_frame_paths)

            try:
                assert abs(timeCalc(telluric_frame) - timeCalc(science_frame)) < telluricTimeThreshold
            except AssertionError as e:
                logging.error("The closest telluric frame in time ({}) for science frame {} has a calculated time delta of {} seconds. This is over the allowed telluric time threshold of {} seconds specified in the config file. Terminating.".format(science_frame, telluric_frame, abs(timeCalc(telluric_frame) - timeCalc(science_frame)), telluricTimeThreshold))
                raise e

            appendToscienceMatchedTellsList(science_frame, telluric_frame)


    # ---------------------------- Tests ------------------------------------- #

    # Don't use tests if user doesn't want them
    tests = True
    if tests:
        # Check that each science observation has valid telluric data.

        # Make list of science frames for a given grating and data
        
        # Get list of unique target, date, grating triples
        target_date_grates = getTargetDateGrates(obsDirList)

        # Now loop over each target, date, grating triple.
        # Each target, date, grating triple needs all science frames to be matched with a telluric in some scienceMatchedTellsList.
        basePath = os.path.sep.join(os.path.normpath(obsDirList[0]).split(os.path.sep)[:-4])
        for target, date, grating in target_date_grates:
            targetDateGratePath = os.path.join(basePath, target, date, grating)
            # Now get all science frames by opening up the scienceFrameLists
            science_frames = getScienceFrames(targetDateGratePath)
            
            # and get all telluric frames by opening up the scienceMatchedTellsLists
            tellurics_frames = getTelluricFrames(targetDateGratePath)

            # Now ensure the two lists are equal
            science_frames.sort()
            tellurics_frames.sort()
            if len(science_frames) != len(tellurics_frames):
                logging.error("A telluric correction was requested but not enough tellurics were found in {}. Terminating.".format(targetDateGratePath))
                raise ValueError
            for science, telluric in zip(science_frames, tellurics_frames):
                if science != telluric:
                    logging.error("Science frame and scienceMatchedTellsLists differ for {}. Terminating.".format(targetDateGratePath))


    # ---------------------------- End Tests --------------------------------- #

    logging.info("\nI am finished matching science images with telluric frames.")
    return

def getTargetDateGrates(obsDirs):
    target_date_grates = []
    for science_directory in obsDirs:
        target = science_directory.split(os.sep)[-4]
        date = science_directory.split(os.sep)[-3]
        grating = science_directory.split(os.sep)[-2]
        if (target, date, grating) not in target_date_grates:
            target_date_grates.append((target, date, grating))
    return target_date_grates


def getScienceFrames(path, full_path=False):
    science_frames = []
    for science_obs_dir in glob.glob(os.path.join(path, "obs*")):
        try:
            with open(os.path.join(path, science_obs_dir, 'scienceFrameList'), 'r') as f:
                for frame in f.readlines():
                    if frame not in science_frames:
                        if full_path:
                            science_frames.append(os.path.join(path, science_obs_dir, frame.rstrip('\n')+'.fits'))
                        else:
                            science_frames.append(frame)
        except IOError as e:
            logging.error("Science directory {} didn't have a scienceFrameList! Terminating sort.".format(os.path.join(path, science_obs_dir)))
            raise e
    try:
        assert len(science_frames) != 0
    except AssertionError as e:
        logging.error("No science frames were found for a given grating and date ({}). This is most likely an error.".format(path))
        raise e

    return science_frames

def getTelluricFrames(path, full_path=False, listToOpen='scienceMatchedTellsList'):
    tellurics_frames = []
    for telluric_obs_dir in glob.glob(os.path.join(path, "Tellurics", "obs*")):
        try:
            with open(os.path.join(path, telluric_obs_dir, listToOpen), 'r') as f:
                for frame in list(filter(lambda x: 'obs' not in x, f.readlines())): # Strips out "obs_n" entries
                    if frame not in tellurics_frames:
                        if full_path:
                            tellurics_frames.append(os.path.join(path, telluric_obs_dir, frame.rstrip('\n')+'.fits'))
                        else:
                            tellurics_frames.append(frame)
        except IOError as e:
            # Not usually a problem if a scienceMatchedTellsList isn't there.
            pass
    return tellurics_frames

def matchScienceTelluric(science_frame_path, telluric_frame_paths):
    # Calculate time of the science frame
    science_time = timeCalc(science_frame_path)

    telluric_time_deltas = [abs(timeCalc(x) - science_time) for x in telluric_frame_paths]

    try:
        assert len(telluric_time_deltas) == len(telluric_frame_paths)
    except AssertionError as e:
        logging.error("timeCalc seems to have failed on a frame. One example is {}. Terminating.".format(telluric_frame_paths[0]))

    # Find closest telluric frame in time.

    return sorted(zip(telluric_time_deltas, telluric_frame_paths))[0][1]

def appendToscienceMatchedTellsList(science_frame_path, telluric_frame_path):
    # Get telluric observation directory where we should find the scienceMatchedTellsList
    tell_obs_directory = os.path.sep.join(os.path.normpath(telluric_frame_path).split(os.path.sep)[:-1])
    science_frame = os.path.normpath(science_frame_path).split(os.path.sep)[-1].rstrip('.fits')
    science_obs = os.path.normpath(science_frame_path).split(os.path.sep)[-2]
    # Append the science frame name and observation to the scienceMatchedTellsList
    try:
        with open(os.path.join(tell_obs_directory, 'scienceMatchedTellsList'), 'a+') as f:
            f.seek(0)
            lines = [x.rstrip('\n') for x in f.readlines()]
    except IOError as e:
        logging.error("Failed to open a scienceMatchedTellsList in {}. Terminating.".format(tell_obs_directory))
        raise e
        
    if science_frame in lines:
        logging.error("A duplicate entry tried to be inserted into a scienceMatchedTellsList in {}. Terminating.".format(tell_obs_directory))
        raise AssertionError

    if science_obs in lines: # Insert somewhere after the occurence of obs in the list
        lines.insert(lines.index(science_obs)+1, science_frame)
    else:
        lines.append(science_obs)
        lines.append(science_frame)

    try:
        with open(os.path.join(tell_obs_directory, 'scienceMatchedTellsList'), 'w') as f:
            for line in lines:
                f.write(line+'\n')
    except IOError as e:
        logging.error("Failed to append to a scienceMatchedTellsList in {}. Terminating.".format(tell_obs_directory))
        raise e

def removeOldScienceMatchedTellsList(telluric_frame_paths):
    # Get a list of telluric observatory directories from the paths
    tell_obs_dirs = [os.path.sep.join(os.path.normpath(x).split(os.path.sep)[:-1]) for x in telluric_frame_paths]

    for tell_obs_dir in tell_obs_dirs:
        if (os.path.exists(os.path.join(tell_obs_dir, 'scienceMatchedTellsList'))):
            os.remove(os.path.join(tell_obs_dir, 'scienceMatchedTellsList'))



#--------------------------- End of Functions ---------------------------------#

if __name__ == '__main__':
    # Set up logging if called as a standalone script.
    # Format logging options.
    FORMAT = '%(asctime)s %(message)s'
    DATEFMT = datefmt()

    # Set up the logging file.
    logging.basicConfig(filename='nifsSort.log',format=FORMAT,datefmt=DATEFMT,level=logging.DEBUG)
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # This lets us logging.info(to stdout AND a logfile. Cool, huh?
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(message)s')
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Start nifsSort from the beginning!
    start()
