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
rewriteSciImageList, datefmt, downloadQueryCadc, CalibrationsNotFoundError, CalibrationsError, TelluricsNotFoundError, \
ScienceObservationError, SkyFrameError, ObservationDirError, WavelengthError, checkForMDFiles

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
            query = "SELECT observationID, publisherID, productID FROM caom2.Observation AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID WHERE instrument_name='NIFS' AND proposal_id={}".format("'" + program + "'")
            downloadQueryCadc(query, os.getcwd()+'/rawData')
            json_query = "https://archive.gemini.edu/jsonsummary/not_site_monitoring/NotFail/{}/notengineering/canonical/present".format(program)
            checkForMDFiles(os.path.join(os.getcwd(), 'rawData'), json_query)
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
        allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, obsidDateList, sciImageList = makePythonLists(rawPath, skyThreshold)
        try:
            objDirList, scienceDirectoryList, telluricDirectoryList = sortScienceAndTelluric(allfilelist, sciImageList, rawPath, skyThreshold)
        except TelluricsNotFoundError:
            logging.warning("Insufficient tellurics found. Turning off telluric correction.", exc_info=True)
            turnOffTelluricCorrectionFluxCalibration()

        scienceDirectoryList, telluricDirectoryList, calibrationDirectoryList = sortCalibrations(arcdarklist, arclist, flatlist, flatdarklist, ronchilist, objectDateGratingList, objDirList, obsidDateList, sciImageList, rawPath, manualMode, dataSource, scienceDirectoryList, telluricDirectoryList)

        # If a telluric reduction will be performed sort the science and telluric images based on time between observations.
        # This will NOT be executed if -t False is specified at command line.
        if sortTellurics:
            try:
                matchTellurics(telluricDirectoryList, scienceDirectoryList, telluricTimeThreshold)
            except TelluricsNotFoundError:
                logging.warning("Insufficient tellurics found. Turning off telluric correction.", exc_info=True)
                turnOffTelluricCorrectionFluxCalibration()

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

    # Fix for arc dark identification, could be done better...
    arc_exp_times = []

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
        ID = header[0].header['OBSID'].split('-')[-1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        aper = header[0].header['APERTURE']
        # If object name isn't alphanumeric, make it alphanumeric.
        objname = header[0].header['OBJECT']
        objname = re.sub('[^a-zA-Z0-9\n\.]', '', objname)
        poff = header[0].header['POFFSET']
        qoff = header[0].header['QOFFSET']
        try:
            exptime = float(header[0].header['EXPTIME'])
        except KeyError:
            logging.warning("EXPTIME keyword missing for a frame, frame {} will not be used.".format(entry))
            continue

        # Make a list of science, telluric and acquisition frames.
        # Use the copied variable (1 not copied, 0 copied) to check later that
        # the file was copied correctly. allfilelist is a 2D list of
        # [[filename1, copied, obsclass1], [filename2, copied, obsclass2]] pairs.
        if obstype == 'OBJECT' and (obsclass == 'science' or obsclass == 'acq' or obsclass == 'acqCal' or obsclass == 'partnerCal'):

            # Append a list of [filename, copied, obsclass] to the list. copied is
            # 1 if not copied, 0 if copied. obsclass is used later for checks.
            templist = [entry, 1, obsclass]
            allfilelist.append(templist)

            # Create sciDateList: list of unique dates of science observations.
            if obsclass == 'science':
                sciImageList.append(entry)
                # Append if list is empty or not a duplicate of last entry.
                if not sciDateList or not sciDateList[-1]==date:
                    sciDateList.append(date)

        # Add arc frame names to arclist and save exp times to identify arc darks.
        elif isArc(obstype):
            arclist.append(entry)
            arc_exp_times.append(exptime)

        # Add lamps on flat frames to flatlist,
        # add lamps off flat frames to flatdarklist,
        # add ronchi flat frames to ronchilist.
        # Lamps on and lamps off flats, and lamps on and lamps off ronchis are
        # seperated by mean number of counts per pixel. Arbitrary threshold is
        # if mean_counts < 500, it is a lamps off flat or ronchi.
        elif isRonchiFlat(obstype, aper, entry):
            ronchilist.append(entry)

        elif isFlat(obstype, aper, entry):
            flatlist.append(entry)

        elif isFlatDark(obstype, aper, entry):
            flatdarklist.append(entry)

    # Do arc dark frames outside of loop b/c need to match by exptime
    # Add arc dark frame names to arcdarklist.
    for entry in rawfiles:
        headers = HeaderInfo(entry)
        if isArcDark(headers.obstype, headers.exptime, arc_exp_times):
            arcdarklist.append(entry)

    # Based on science (including sky) frames, make a list of unique [object, date] list pairs to be used later.
    for i in range(len(rawfiles)):
        headers = HeaderInfo(rawfiles[i])
        
        if headers.obsclass == 'science':
            list1 = [headers.objname, headers.date, headers.grat]
            # Append if list is empty or not a duplicate of last entry.
            if not objectDateGratingList or not list1 in objectDateGratingList:
                objectDateGratingList.append(list1)

    # Make list of unique [date, obsid] pairs from FLATS. If flat was taken on the same day as a science
    # frame, append that flat date. If not, append an arbitrary unique date from sciDateList.
    # This is so we can sort calibrations later by date and observation id.
    n = 0
    for flat in flatlist:
        headers = HeaderInfo(flat)
        # Make sure no duplicate dates are being entered.
        if flatlist.index(flat)==0 or not oldobsid==headers.ID:
            #if date in sciDateList:
            list1 = [headers.date, headers.ID]
            obsidDateList.append(list1)
            #else:
                # Ugly fix, we have to check there aren't more flats than science dates.
            #    if n < len(sciDateList):
            #        list1 = [sciDateList[n], obsid]
             #       obsidDateList.append(list1)
            #n+=1
        oldobsid = headers.ID

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

    # Store number of telluric, telluric, sky, telluric sky and acquisition frames in number_files_to_be_copied.
    number_files_to_be_copied = len(allfilelist)

    # Store number of arcs, arc darks, lamps on flats, lamps off flats and ronchi flats in number_calibration_files_to_be_copied.
    number_calibration_files_to_be_copied = len(arclist) + len(arcdarklist) +\
                         len(flatlist) + len(flatdarklist) + len(ronchilist)

    logging.info("\nTotal number of frames to be copied: " + str(number_files_to_be_copied + number_calibration_files_to_be_copied))

    return allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, obsidDateList, sciImageList

#----------------------------------------------------------------------------------------#

def sortScienceAndTelluric(allfilelist, sciImageList, rawPath, skyThreshold):

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
                if not objDirList or not objDir in objDirList :
                    objDirList.append(objDir)
            else:
                objDir = path+'/'+objname+'/'+date
                if not objDirList or not objDir in objDirList:
                    objDirList.append(objDir)

    # For each science frame, create a "science_object_name/date/grating/observationid/"
    # directory in the current working directory.
    for entry in allfilelist:
        header = astropy.io.fits.open(rawPath+'/'+entry[0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'].split('-')[-1]
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
            # Else if a new list or not a duplicate append time and directory name to scienceDirList.
            if not scienceDirList or not filter(lambda x: x[1] == objDir+'/'+date+'/'+grat+'/obs'+obsid, scienceDirList):
                scienceDirList.append([[time], objDir+'/'+date+'/'+grat+'/obs'+obsid])
            # IF A DUPLICATE:
            # Append the time to an existing time list.
            #
            # Eg, before appending: [[[5400,6500,7200], '/path/to/first/science'], [[5200], '/path/to/second/science']]
            # If we are processing the second frame in /path/to/second/science/, scienceDirList[-1][1] will equal /path/to/second/science.
            # Then append the new time to the second entry.
            #
            # [[[5400,6500,7200], '/path/to/first/science'], [[5200, NEWTIMEHERE], '/path/to/second/science']]

            else:
                index = scienceDirList.index(filter(lambda x: x[1] == objDir+'/'+date+'/'+grat+'/obs'+obsid, scienceDirList)[0])
                scienceDirList[index][0].append(time)


    # Copy science and acquisition frames to the appropriate directory.
    logging.info("\nCopying Science and Acquisitions.\nCopying science frames and science acquisitions.\nNow copying: ")

    for i in range(len(allfilelist)):
        header = astropy.io.fits.open(rawPath+'/'+allfilelist[i][0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'].split('-')[-1]
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

        # Copy the most recent acquisition in each set to a new directory to be optionally
        # used later by the user for checks (not used by the pipeline).
        if obsclass=='acq':
            logging.info(allfilelist[i][0])
            # Find appropriate directory to copy to by doing a time calculation
            try:
                basePath = getBasePathWithTimes(scienceDirList, os.path.join(rawPath, allfilelist[i][0]))
                # create an Acquisitions directory in objDir/YYYYMMDD/grating
                if not os.path.exists(os.path.join(basePath, 'Acquisitions')):
                    os.makedirs(os.path.join(basePath, 'Acquisitions'))
                #if copyMostRecentAcquisition(rawPath +'/'+ allfilelist[i][0], path+'/'+obj2+'/'+date+'/'+grat+'/Acquisitions/' + allfilelist[i][0]):
                shutil.copy(rawPath+'/'+allfilelist[i][0], os.path.join(basePath, 'Acquisitions'))
                number_files_that_were_copied += 1
                allfilelist[i][1] = 0
            except IndexError:
                pass

    # Copy telluric frames to the appropriate folder.
    # Note: Because the 'OBJECT' of a telluric file header is different then the
    # science target, we need to sort by date, grating AND most recent time.
    logging.info("\nCopying telluric frames.\nNow copying: ")
    for i in range(len(allfilelist)):
        header = astropy.io.fits.open(rawPath+'/'+allfilelist[i][0])

        obstype = header[0].header['OBSTYPE'].strip()
        obsid = header[0].header['OBSID'].split('-')[-1]
        grat = header[0].header['GRATING'][0:1]
        date = header[0].header[ 'DATE'].replace('-','')
        obsclass = header[0].header['OBSCLASS']
        obj = header[0].header['OBJECT']
        obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)
        telluric_time = timeCalc(rawPath+'/'+allfilelist[i][0])


        if isTelluric(obstype, obsclass):
            logging.info(allfilelist[i][0])
            headers = HeaderInfo(os.path.join(rawPath, allfilelist[i][0]))
            try:
                path_to_tellurics = getBasePathWithTimes(scienceDirList, os.path.join(rawPath, allfilelist[i][0]))
            except IndexError:
                continue

            # Create a Tellurics directory in science_object_name/YYYYMMDD/grating.

            if not os.path.exists(path_to_tellurics + '/Tellurics'):
                os.mkdir(path_to_tellurics + '/Tellurics')
            # Create an obsid (eg. obs25) directory in the Tellurics directory.
            if not os.path.exists(path_to_tellurics+'/Tellurics/obs'+obsid):
                os.mkdir(path_to_tellurics+'/Tellurics/obs'+obsid)
                telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
            if not telDirList or not path_to_tellurics+'/Tellurics/obs'+obsid in telDirList:
                telDirList.append(path_to_tellurics+'/Tellurics/obs'+obsid)
            shutil.copy(rawPath+'/'+allfilelist[i][0], path_to_tellurics+'/Tellurics/obs'+obsid+'/')
            number_files_that_were_copied += 1
            allfilelist[i][1] = 0

    """
    # Run another copy loop to make sure telluric sky frames get copied over.
    for frame_obj in allfilelist:
        frame = frame_obj[0]
        headers = HeaderInfo(os.path.join(rawPath, frame))
        if isTelluricSky(headers.obstype, headers.obsclass, headers.poff, headers.qoff, 2.0):
            # Find the telluric dir it originally got copied to.
            for telluric_directory in telDirList:
                if os.path.exists(os.path.join(telluric_directory, frame)):
                    tells_prefix = os.path.sep.join(os.path.normpath(telluric_directory).split(os.path.sep)[:-1])
                    break
            else:
                logging.warning("Telluric sky frame {} not found to have been copied.".format(frame))
            
            # Copy this telluric frame to all other telluric observation directories with the same target/date/grating triple.
            telluric_obs = glob.glob(os.path.join(tells_prefix, "obs*"))
            for tell_obs_dir in telluric_obs:
                shutil.copy(os.path.join(rawPath, frame), os.path.join(tells_prefix, tell_obs_dir))
    """
    # Modify scienceDirList to a format telSort can use.
    tempList = []
    for i in range(len(scienceDirList)):
        tempList.append(scienceDirList[i][1])
    scienceDirList = tempList

    # make skyFrameList/scienceFrameList, and that directory contains a non-empty skyFrameList and scienceFrameList.
    for science_directory in list(scienceDirList):
        try:
            makeSkyLists(science_directory, skyThreshold, science=True)
            checkSkyFrameDivision(science_directory, science=True)
        except ObservationDirError as e:
            logging.error("Science directory {} has a problem with the scienceFrameList or skyFrameList. Removing that directory from the list of directories to reduce.".format(science_directory), exc_info=True)
            scienceDirList = [x for x in scienceDirList if x != science_directory]
        try:
            checkStandardWavelength(science_directory)
        except WavelengthError as e:
            logging.error("Science directory {} has frames with a non-standard wavelength. Terminating.".format(science_directory), exc_info=True)
            raise e

    # Check that telluric directory contains a non-empty skyFrameList and scienceFrameList.
    for telluric_directory in list(telDirList):
        try:
            makeSkyLists(telluric_directory, skyThreshold, science=False)
            checkSkyFrameDivision(telluric_directory, science=False)
        except SkyFrameError as e:
            logging.warning("Possibly no telluric sky frames found in {}. Turning off telluric sky subtraction for all telluric observations.".format(telluric_directory))
            turnOffTelluricSkySub()
        except ObservationDirError as e:
            logging.error("Telluric observation directory {} has a problem with the tellist or skyFrameList. Removing that directory from the list of directories to reduce.".format(telluric_directory), exc_info=True)
            telDirList = [x for x in telDirList if x != telluric_directory]
        try:
            checkStandardWavelength(telluric_directory)
        except WavelengthError as e:
            logging.error("Telluric directory {} has frames with a non-standard wavelength. Terminating.".format(science_directory), exc_info=True)
            raise e

    #------------------------------ TESTS -------------------------------------#

    # Check that each science frame was copied to a directory.
    for science_frame in sciImageList:
        headers = HeaderInfo(os.path.join(rawPath, science_frame))
        science_dir = os.path.join(os.getcwd(), headers.objname, headers.date, headers.grat, 'obs'+headers.obsid)
        if not os.path.isdir(science_dir):
            logging.error("Science frame {} wasn't copied to a directory; {} was supposed to exist but wasn't found.".format(science_frame, science_dir))

    # Check to see which files were not copied.
    logging.info("\nChecking for non-copied science, tellurics and acquisitions.\n")
    for i in range(len(allfilelist)):
        # Check the copied flag. If not 0, logging.info("the entry.")
        if allfilelist[i][1] != 0:
            logging.info(str(allfilelist[i][0]) + " " + str(allfilelist[i][2]) + " was not copied.")
    logging.info("\nEnd non-copied science, tellurics and acquisitions.\n")

    # Check that all science frames were copied.
    count_from_raw_files = len(sciImageList)

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


    # Test for telluric directories with mis-identified sky frames. For now, we identify sky frames
    # based on absolute P and Q offsets, not relative to a zero point. This can cause problems.
    # TODO(nat): look into if it is worth it to use relative P and Q offsets.
    # If there's a problem with the tellurics directory, try to fix it; if we can't, log a warning and turn off the telluric reduction.
    """
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
                tellImageList = open('tellist', "r").readlines()
                if len(tellImageList) > 0:
                    logging.info("\nSucceeded; a telluric frame list exists in " + str(os.getcwd()))
                else:
                    raise IOError()
            except IOError:
                logging.error("\nWARNING: no telluric frames found in " + str(os.getcwd()) + ". You may have to adjust the skyThreshold parameter.")
                raise TelluricsNotFoundError()
            os.chdir(path)
    """
    logging.info("\nDone sorting and copying science and tellurics. Moving on to Calibrations.\n")


    return objDirList, scienceDirList, telDirList

#----------------------------------------------------------------------------------------#

def sortCalibrations(arcdarklist, arclist, flatlist, flatdarklist, ronchilist, objectDateGratingList, objDirList, obsidDateList, sciImageList, rawPath, manualMode, dataSource, scienceDirectoryList, telluricDirectoryList):

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
        obsid = header[0].header['OBSID'].split('-')[-1]
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
        obsid = header[0].header['OBSID'].split('-')[-1]
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
        obsid = header[0].header['OBSID'].split('-')[-1]
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
        obsid = header[0].header['OBSID'].split('-')[-1]
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


    logging.info("\nChecking that each science image has required calibration data. ")

    for science_frame in sciImageList:
        try:
            checkCalibrationsPresent(rawPath, science_frame, flatlist, flatdarklist, arclist, arcdarklist, ronchilist, dataSource)
        except CalibrationsError:
            logging.error("Calibrations not found for science frame {}. Unmarking all affected science and telluric telluric directories for reduction.".format(science_frame), exc_info=True)
            scienceDirectoryList, telluricDirectoryList, calDirList = removeAffectedDirectories(os.path.join(rawPath, science_frame), scienceDirectoryList, telluricDirectoryList, calDirList)
        os.chdir(path1)

    # Change back to original working directory.
    os.chdir(path1)

    # ---------------------------- End Tests --------------------------------- #

    return scienceDirectoryList, telluricDirectoryList, calDirList

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

    try:
        basePath = os.path.sep.join(os.path.normpath(obsDirList[0]).split(os.path.sep)[:-4])
    except IndexError:
        logging.error("No tellurics directories were provided.")
        raise TelluricsNotFoundError()
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

            try:
                telluric_frame = matchScienceTelluric(science_frame, telluric_frame_paths)
            except TelluricsNotFoundError as e:
                logging.error("Telluric frame not found for science frame {}. A possible solution is to raise the telluricTimeThreshold in the config file (it is currently {} seconds).".format(science_frame, telluricTimeThreshold))
                raise e


            if abs(timeCalc(telluric_frame) - timeCalc(science_frame)) > telluricTimeThreshold:
                raise TelluricsNotFoundError("The closest telluric frame in time ({}) for science frame {} has a calculated time delta of {} seconds. This is over the allowed telluric time threshold of {} seconds specified in the config file.".format(science_frame, telluric_frame, abs(timeCalc(telluric_frame) - timeCalc(science_frame)), telluricTimeThreshold))

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
                raise TelluricsNotFoundError()
            for science, telluric in zip(science_frames, tellurics_frames):
                if science != telluric:
                    logging.error("Science frame and scienceMatchedTellsLists differ for {}. Terminating.".format(targetDateGratePath))
                    raise TelluricsNotFoundError()


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
    if len(zip(telluric_time_deltas, telluric_frame_paths)) == 0:
        raise TelluricsNotFoundError()

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


def checkListExists(science_frame, in_list, list_name, list_description, rawPath):
    # $list_name exists and has more than one file.
    try:
        listFile = open(list_name, "r").readlines()
        if len(listFile) == 1:
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING in sort: only 1 {} frame found for science".format(list_description))
            logging.info("                      frame "+str(science_frame))
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")
    except IOError:
        logging.info("\n#####################################################################")
        logging.info("#####################################################################")
        logging.info("")
        logging.info("     WARNING in sort: no {} found for science frame".format(list_name))
        logging.info("                      "+str(science_frame))
        logging.info("")
        logging.info("#####################################################################")
        logging.info("#####################################################################\n")

        # Sometimes calibration frames can be taken a day after the observing night. First
        # look for these, and if they are not found, ask the user to provide some calibration frames.
        # Get date before and after the science observation
        newdate_before, newdate_after = getPlusMinusDays(os.path.join(rawPath, science_frame), 1)
        # Loop through in_list and see if there is a frame taken on this date, optionally checking that exposure times match
        for i in range(len(in_list)):
            header = astropy.io.fits.open(rawPath+'/'+in_list[i][0])
            date = header[0].header['DATE'].replace('-','')

            if (str(date) == newdate_after.strftime('%Y%m%d')) or (str(date) == newdate_before.strftime('%Y%m%d')):
                # If so, copy it to the appropriate calibrations directory and write in the text list.
                shutil.copy(rawPath + '/' + in_list[i][0], './')
                writeList(in_list[i][0], list_name, os.getcwd())
                logging.info("\n#####################################################################")
                logging.info("#####################################################################")
                logging.info("")
                logging.info("     WARNING in sort: found a {} taken one day before/after a science frame.".format(list_description))
                logging.info("                      "+str(science_frame))
                logging.info("                       using that.")
                logging.info("")
                logging.info("#####################################################################")
                logging.info("#####################################################################\n")
                break
                in_list[i][1] = 0
        else:
            raise RuntimeError("Error: no {} found in {}; some {} frames are missing.".format(list_name, os.getcwd(), list_description))


def tryDownloadPlusMinusOneDay(rawPath, science_frame, list_name, data_type, dataSource):
    # Query data source for calibrations plus or minus one day of a given science frame.
    if dataSource == 'CADC':
        day_before, day_after = getPlusMinusDays(os.path.join(rawPath, science_frame), 1)
        newdate_before_mjd = astropy.time.Time(day_before).mjd
        # CADC is exclusive of upper bound so need to add a day.
        newdate_after_mjd = astropy.time.Time(day_after+datetime.timedelta(1)).mjd

        query = "SELECT observationID, publisherID, productID FROM caom2.Observation " + \
                         "AS o JOIN caom2.Plane AS p ON o.obsID=p.obsID " + \
                         "WHERE  ( type = '{}' AND instrument_name = 'NIFS' ".format(data_type) + \
                         "AND collection = 'GEMINI' " + \
                         "AND INTERSECTS( INTERVAL( {}, {} ), ".format(newdate_before_mjd, newdate_after_mjd) + \
                         "p.time_bounds ) = 1 )"
        try:
            downloadQueryCadc(query, directory=os.getcwd())
        except Exception as e:
            logging.warning("A problem occured while trying to download calibrations one day before/after science frame {} in {}.".format(science_frame, os.getcwd()))
            raise e
    else:
        logging.error("Querying for calibrations plus/minus one day is unsupported with the dataSource={} option. Current support is for the 'CADC' option.".format(dataSource))
        raise ValueError

    # Now see if the downloaded files can be used to augment the list that's missing the files.
    # rewriteCalibrationList() returns num calibrations found
    if not rewriteCalibrationList(list_name):
        logging.error("No new calibrations were found for science frame {} in {} by downloading calibrations one day ahead and behind.".format(science_frame, os.getcwd()))
        raise RuntimeError


def getPlusMinusDays(science_frame, num_days):
    """
    science_frame: absolute path to science frame.
    num_days: number of UT days in radius to search.
    """
    sci_header = astropy.io.fits.open(science_frame)
    sci_date = sci_header[0].header['DATE'].replace('-','')
    t = time.strptime(sci_date,'%Y%m%d')
    newdate_after=datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday)+datetime.timedelta(num_days)
    newdate_before=datetime.datetime(t.tm_year,t.tm_mon,t.tm_mday)-datetime.timedelta(num_days)
    return (newdate_before, newdate_after)


def rewriteCalibrationList(list_name):
    # rewrites a calibration list based on the files in the current working directory, counting number of files written.
    count = 0
    frames = glob.glob("N2*")
    if list_name == 'flatlist':
        with open('flatlist', 'w') as f:
            for frame in frames:
                headers = HeaderInfo(frame)
                if isFlat(headers.obstype, headers.aper, frame):
                    writeList(frame, 'flatlist', os.getcwd())
                    count += 1
    elif list_name == 'arclist':
        with open('arclist', 'w') as f:
            for frame in frames:
                headers = HeaderInfo(frame)
                if isArc(headers.obstype):
                    writeList(frame, 'arclist', os.getcwd())
                    count += 1
    elif list_name == 'arcdarklist':
        # First open the list of arcs to get exptimes of them
        try:
            with open('arclist') as f:
                lines = f.readlines()
                frames = [x.rstrip('\n')+'.fits' for x in frames]
            arc_exp_times = [HeaderInfo(frame).exptime for frame in frames]
        except Exception:
            arc_exp_times = None
        with open('arcdarklist', 'w') as f:
            for frame in frames:
                headers = HeaderInfo(frame)
                if isArcDark(headers.obstype, headers.exptime, arc_exp_times):
                    writeList(frame, 'arcdarklist', os.getcwd())
                    count += 1
    elif list_name == 'flatdarklist':
        with open('flatdarklist', 'w') as f:
            for frame in frames:
                headers = HeaderInfo(frame)
                if isFlatDark(headers.obstype, headers.aper, frame):
                    writeList(frame, 'flatdarklist', os.getcwd())
                    count += 1
    elif list_name == 'ronchilist':
        with open('ronchilist', 'w') as f:
            for frame in frames:
                headers = HeaderInfo(frame)
                if isRonchiFlat(headers.obstype, headers.aper, frame):
                    writeList(frame, 'ronchilist', os.getcwd())
                    count += 1
    else:
        logging.error("Invalid list name.")
        raise ValueError
    assert count <= len(frames)
    return count


def isFlat(obstype, aper, frame):
    if (obstype == 'FLAT') and (aper != 'Ronchi_Screen_G5615'):
        array = astropy.io.fits.getdata(frame)
        mean_counts = np.mean(array)
        return mean_counts > 500
    return False

def isArc(obstype):
    return obstype == 'ARC'

def isFlatDark(obstype, aper, frame):
    if (obstype == 'FLAT') and (aper != 'Ronchi_Screen_G5615'):
        array = astropy.io.fits.getdata(frame)
        mean_counts = np.mean(array)
        return mean_counts < 500
    return False

def isArcDark(obstype, exptime, exptimes=None):
    """
    Arc darks are hard because there's no header info that distinguishes them from biases.
    Thus, we need to compare exptimes with lamps on arcs (assuming no collisions?) to identify.
    """
    if exptimes:
        return (obstype == 'DARK') and (exptime in exptimes)
    # Don't do an in-depth check if no exptimes.
    return obstype == 'DARK'

def isRonchiFlat(obstype, aper, frame):
    # Once the mean is stored in mean_counts we can check whether the
    # frame is a lamps off ronchi or a lamps on ronchi based on the counts.
    # 500.0 is an arbitrary threshold that appears to work well.
    if (obstype == 'FLAT') and (aper == 'Ronchi_Screen_G5615'):
        array = astropy.io.fits.getdata(frame)
        mean_counts = np.mean(array)
        return mean_counts > 500
    return False

def isTelluric(obstype, obsclass):
    return (obstype=='OBJECT') and (obsclass=='partnerCal')

def isTelluricSky(obstype, obsclass, poff, qoff, telluricSkyThreshold):
    return (obstype == 'OBJECT') and (obsclass == 'partnerCal') and math.sqrt((poff**2)+(qoff**2)) < telluricSkyThreshold

class HeaderInfo(object):
    def __init__(self, frame):
        try:
            header = astropy.io.fits.open(frame)

            # Store information in variables.
            self.instrument = header[0].header['INSTRUME']
            if self.instrument != 'NIFS':
                # Only grab frames belonging to NIFS raw data!
                raise ValueError("Data isn't from NIFS!")
            self.obstype = header[0].header['OBSTYPE'].strip()
            self.ID = header[0].header['OBSID'].split('-')[-1]
            self.grat = header[0].header['GRATING'][0:1]
            self.date = header[0].header[ 'DATE'].replace('-','')
            self.obsclass = header[0].header['OBSCLASS']
            self.aper = header[0].header['APERTURE']
            # If object name isn't alphanumeric, make it alphanumeric.
            self.objname = re.sub('[^a-zA-Z0-9\n\.]', '', header[0].header['OBJECT'])
            self.obsid = header[0].header['OBSID'].split('-')[-1]
            self.poff = header[0].header['POFFSET']
            self.qoff = header[0].header['QOFFSET']
            self.exptime = float(header[0].header['EXPTIME'])
            self.crWav = float(header[0].header['GRATWAVE'])
        except Exception as e:
            logging.error("Error getting header info for frame {}.".format(frame))
            raise e

def turnOffTelluricCorrectionFluxCalibration():
    with open('./config.cfg') as config_file:
        options = ConfigObj(config_file, unrepr=True)
    options['nifsPipelineConfig']['telluricReduction'] = False
    options['nifsPipelineConfig']['telluricCorrection'] = False
    options['nifsPipelineConfig']['fluxCalibration'] = False
    with open('./config.cfg', 'w') as config_file:
        options.write(config_file)

def checkCalibrationsPresent(rawPath, science_frame, flatlist, flatdarklist, arclist, arcdarklist, ronchilist, dataSource):
    header = astropy.io.fits.open(rawPath+'/'+science_frame)

    obstype = header[0].header['OBSTYPE'].strip()
    obsid = header[0].header['OBSID'].split('-')[-1]
    grat = header[0].header['GRATING'][0:1]
    date = header[0].header[ 'DATE'].replace('-','')
    obsclass = header[0].header['OBSCLASS']
    obj = header[0].header['OBJECT']
    obj = re.sub('[^a-zA-Z0-9\n\.]', '', obj)

    # a science and Calibrations directory are present.
    try:
        os.chdir(os.getcwd()+'/'+obj+'/'+date+'/'+grat+'/obs'+obsid+'/')
        os.chdir('../../Calibrations_'+grat+'/')
    except OSError as e:
        logging.info("\n#####################################################################")
        logging.info("#####################################################################")
        logging.info("")
        logging.info("     WARNING in sort: no Calibrations directory found for ")
        logging.info("                      science frame "+str(science_frame))
        logging.info("")
        logging.info("#####################################################################")
        logging.info("#####################################################################\n")
        raise CalibrationsNotFoundError("Calibrations directory not found in {}.".format(os.path.join(os.getcwd(), obj, date, grat, 'obs'+obsid)))

    try:
        checkListExists(science_frame, flatlist, 'flatlist', 'lamps on flat', rawPath)
    except RuntimeError as e_list:
        try:
            tryDownloadPlusMinusOneDay(rawPath, science_frame, 'flatlist', 'FLAT', dataSource)
        except Exception as e_download:
            raise CalibrationsNotFoundError(
                [
                    e_list,
                    e_download
                ]
            )
    try:
        checkListExists(science_frame, flatdarklist, 'flatdarklist', 'lamps off flat', rawPath)
    except RuntimeError as e_list:
        try:
            tryDownloadPlusMinusOneDay(rawPath, science_frame, 'flatdarklist', 'DARK', dataSource)
        except Exception as e_download:
            raise CalibrationsNotFoundError(
                [
                    e_list,
                    e_download
                ]
            )
    # Make sure flatlist and flatdarklist are the same length. nsflat() complains otherwise.
    checkSameLengthFlatLists()

    try:
        checkListExists(science_frame, arclist, 'arclist', 'arc', rawPath)
    except RuntimeError as e_list:
        try:
            tryDownloadPlusMinusOneDay(rawPath, science_frame, 'arclist', 'ARC', dataSource)
        except Exception as e_download:
            raise CalibrationsNotFoundError(
                [
                    e_list,
                    e_download
                ]
            )
    try:
        checkListExists(science_frame, arcdarklist, 'arcdarklist', 'arc dark', rawPath)
    except RuntimeError as e_list:
        try:
            tryDownloadPlusMinusOneDay(rawPath, science_frame, 'arcdarklist', 'DARK', dataSource)
        except Exception as e_download:
            raise CalibrationsNotFoundError(
                [
                    e_list,
                    e_download
                ]
            )
    try:
        checkListExists(science_frame, ronchilist, 'ronchilist', 'ronchi flat', rawPath)
    except RuntimeError as e_list:
        try:
            tryDownloadPlusMinusOneDay(rawPath, science_frame, 'ronchilist', 'RONCHI', dataSource)
        except Exception as e_download:
            raise CalibrationsNotFoundError(
                [
                    e_list,
                    e_download
                ]
            )

def removeAffectedDirectories(science_frame_absolute, scienceDirectoryList, telluricDirectoryList, calDirList):
    """ Takes an absolute path to a science frame and removes all relevant science observation, telluric observation,
    and calibration directories from the list of things to reduce.
    """
    headers = HeaderInfo(science_frame_absolute)

    tells_and_science_path = os.path.join(os.path.split(os.getcwd())[0], headers.grat)
    calibrations_path = os.getcwd()

    # Remove science and telluric observation directories
    scienceDirectoryList = [x for x in scienceDirectoryList if tells_and_science_path not in x]
    telluricDirectoryList = [x for x in telluricDirectoryList if tells_and_science_path not in x]
    calDirList = [x for x in calDirList if calibrations_path not in x]

    return scienceDirectoryList, telluricDirectoryList, calDirList


def checkSkyFrameDivision(observation_dir, science=True):
    if science:
        framelist = 'scienceFrameList'
    else:
        framelist = 'tellist'

    try:
        with open(os.path.join(observation_dir, framelist), 'r') as f:
            object_frames = f.readlines()
        assert len(object_frames) > 0
    except (IOError,  AssertionError):
        raise ObservationDirError("There was a problem with the object frame list (scienceFrameList or tellist) in {}.".format(observation_dir))
    
    try:
        with open(os.path.join(observation_dir, 'skyFrameList'), 'r') as f:
            sky_frames = f.readlines()
    except IOError:
        raise SkyFrameError("There was a problem with the skyFrameList in {}.".format(observation_dir))
    if science:
        try:
            assert len(sky_frames) > 0
        except AssertionError:
            raise SkyFrameError("There was a problem with the skyFrameList in {}.".format(observation_dir))


def makeSkyLists(science_dir, skyThreshold, science=True):
    frames = glob.glob(os.path.join(science_dir, "N2*"))

    object_frame_count = 0
    sky_frame_count = 0

    for frame in frames:
        headers = HeaderInfo(frame)
        radii = math.sqrt((headers.poff**2)+(headers.qoff**2))
        if science:
            if radii < skyThreshold:
                object_frame_count += 1
                writeList(os.path.split(frame)[1].rstrip('.fits'), 'scienceFrameList', os.path.split(frame)[0])
            else:
                sky_frame_count += 1
                writeList(os.path.split(frame)[1].rstrip('.fits'), 'skyFrameList', os.path.split(frame)[0])
        else:
            if radii < skyThreshold:
                object_frame_count += 1
                writeList(os.path.split(frame)[1].rstrip('.fits'), 'tellist', os.path.split(frame)[0])
            else:
                sky_frame_count += 1
                writeList(os.path.split(frame)[1].rstrip('.fits'), 'skyFrameList', os.path.split(frame)[0])

    try:
        assert sky_frame_count < int(3.0*object_frame_count)
    except AssertionError:
        logging.error("Sky frame count was {} which is >= 3.0* Object Frame Count ({}) in {}.".format(sky_frame_count, object_frame_count, science_dir))
        raise SkyFrameError

def turnOffTelluricSkySub():
    with open('./config.cfg') as config_file:
        options = ConfigObj(config_file, unrepr=True)
    options['telluricReductionConfig']['telluricSkySubtraction'] = False
    with open('./config.cfg', 'w') as config_file:
        options.write(config_file)

def copyMostRecentAcquisition(rawPath, destPath):
    # Determine if this acq is the most recent in the given directory.
    aqcDir = os.path.split(destPath)[0]
    aqcs = glob.glob("N2*")
    times = [timeCalc(os.path.join(aqcDir, x)) for x in aqcs]

    # Should only be most recent acq in there
    assert len(times) <= 1
    new_time = timeCalc(rawPath)

    if new_time > times:
        shutil.copy(rawPath, destPath)

def getBasePathWithTimes(scienceDirList, framePath):
    # Look at all science dir lists and find the appropriate one based on time.
    
    # First narrow by date and grating.
    headers = HeaderInfo(framePath)
    filtered_scienceDirList = [x for x in scienceDirList if ((headers.date in x[1].split(os.path.sep)[-3]) and (headers.grat in x[1].split(os.path.sep)[-2]))]
    
    # Find time of frame we're looking at
    frame_time = timeCalc(framePath)

    # Find science dir with the closest time and return that base path.
    for entry in filtered_scienceDirList:
        entry[0] = [abs(frame_time - x) for x in entry[0]]
        entry[0].sort()

    try:
        scienceDirList = sorted(filtered_scienceDirList, key=lambda x: x[0][0])
        basePath = os.path.sep.join(scienceDirList[0][1].split(os.path.sep)[:-1])
    except IndexError as e:
        logging.warning("getBasePathWithTimes() failed to find a base path for frame {}.".format(framePath))
        raise e

    return basePath


def checkStandardWavelength(observationDirectory):
    frames = glob.glob(os.path.join(observationDirectory, "N20*"))
    for frame in frames:
        headers = HeaderInfo(frame)

        if  (('Z' in headers.grat) and (abs(headers.crWav-1.05) < 0.01)) or \
            (('J' in headers.grat) and (abs(headers.crWav-1.25) < 0.01)) or \
            (('H' in headers.grat) and (abs(headers.crWav-1.65) < 0.01)) or \
            (('K' in headers.grat) and (abs(headers.crWav-2.20) < 0.01)):
            return
    raise WavelengthError("A non standard wavelength configuration of {} grating and {} central wavelength was detected in frame {}.".format(headers.grat, headers.crWav, frame))


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
