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

import time, sys, calendar, astropy.io.fits, urllib, shutil, glob, os, fileinput, logging, smtplib, pkg_resources, math, re, collections, requests, hashlib, tempfile, json
import astropy.io.fits as fits
import numpy as np
from xml.dom.minidom import parseString
from pyraf import iraf
from astroquery.cadc import Cadc

# LOCAL

# Import config parsing.
from configobj.configobj import ConfigObj

# Define constants
# Paths to Nifty data.
RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
RUNTIME_DATA_PATH = pkg_resources.resource_filename('nifty', 'runtimeData/')

# Let us print colors to the terminal!
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

#--------------------------------------------------------------------#
#                                                                    #
#     DEFS                                                           #
#                                                                    #
#    Library of non reduction or sorting specific functions          #
#                                                                    #
#    The following functions were taken from the IPM scripts from    #
#    globaldefs.py:                                                  #
#    datefmt, getFitsHeader, FitsKeyEntry, stripString               #
#    stripNumber, getURLFiles                                        #
#                                                                    #
#--------------------------------------------------------------------#


#-----------------------------------------------------------------------------#

def interactiveNIFSInput():
    """
    Get NIFS configuration interactively. This is based on Sphinx's interactive input session.

    """

    logging.info("\nWelcome to Nifty! The current mode is NIFS data reduction.\n\nPress enter to accept default data reduction options. Type 'yes' or 'no' when prompted.")

    fullReduction = getParam(
                "Do a full data reduction with default parameters loaded from recipes/defaultConfig.cfg? [no]: ",
                False,
                "Type yes to start Nifty with data reduction input parameters loaded from recipes/defaultConfig.cfg file."
    )
    if fullReduction == False:
        # "Select in". User has to turn individual steps on.
        # TODO(nat): Implement these steps.

        manualMode = getParam(
        "Run in manualMode? [no]: ",
        False,
        "nifsPipeline can be run as a fully automatic pipeline or pausing after each step, giving you "
        "a chance to review intermediate products."
        )

        over = getParam(
        "Overwrite old files? [no]: ",
        False,
        "Would you like to overwrite any old results Nifty finds?"
        )

        extractionXC = getParam(
        "Extraction X Coordinate? [15.0]: ",
        15.0,
        "You can set the nfextract x coordinate here."
        )
        extractionYC = getParam(
        "Extraction Y Coordinate? [33.0]: ",
        33.0,
        "You can set the nfextract y coordinate here."
        )
        extractionRadius = getParam(
        "Extraction y coordinate? [2.5]: ",
        2.5,
        "You can specify a radius for nfextract here."
        )
        scienceOneDExtraction = getParam(
        "Extract one D spectra from science cubes? [yes]: ",
        'yes',
        "Nifty can extract summed, combined one-D spectra from science cubes as well as telluric cubes. The extraction "
        "radius, x coordinate and y coordinate will be the same as in the telluric one D extraction."
        )
        scienceDirectoryList = getParam(
        "Science Directory List? []: ",
        [],
        "You can specify a Python list of science observation directories here."
        )
        telluricDirectoryList = getParam(
        "Telluric Directory List? []: ",
        [],
        "You can specify a Python list of telluric observation directories here."
        )
        calibrationDirectoryList = getParam(
        "Calibration Directory List? []: ",
        [],
        "You can specify a Python list of calibration observation directories here."
        )

        sort = getParam(
        "Sort raw data? [yes]: ",
        'yes',
        "Nifty needs raw data to be in a specific directory structure with text files helping it. " + \
        "You can have it build a structure automatically."
        )
        # See if we want to reduce the baseline calibrations. And if so, which substeps
        # to perform.
        calibrationReduction = getParam(
        "Reduce baseline calibrations? [yes]: ",
        'yes',
        "NIFS data comes with a set of calibrations that should be reduced and used."
        )
        # Check and get params for a telluric data reduction.
        telluricReduction = getParam(
        "Reduce telluric data? [yes]: ",
        'yes',
        "You can specify to do a telluric data reduction."
        )
        # Check for science as well.
        scienceReduction = getParam(
        "Reduce science data? [yes]: ",
        'yes',
        "Nifty can reduce science data producing UNMERGED 3D data cubes."
        )
        telluricCorrection = getParam(
        "Do a telluric correction? [yes]: ",
        'yes',
        "Nifty can derive and apply a telluric correction using various methods. Currently the only method "
        "implemented is that of Mason et al. 2015. in the XDGNIRS pipeline."
        )
        fluxCalibration = getParam(
        "Do a flux calibration? [yes]: ",
        'yes',
        "Nifty can derive and apply a flux calibration using various methods. Currently the only method "
        "implemented is that of Mason et al. 2015. in the XDGNIRS pipeline."
        )
        merge = getParam(
        "Do a final data cube merging? [yes]: ",
        'yes',
        "Nifty can merge final data cubes by object, grating and correction type."
        )
        telluricCorrectionMethod = getParam(
        "Telluric correction method? [gnirs]: ",
        "gnirs",
        "Specify a telluric correction method or not to apply one. The options are \"gnirs\", \"none\" and \"iraf\"."
        )
        fluxCalibrationMethod = getParam(
        "Flux calibration method? [gnirs]: ",
        "gnirs",
        "Nifty can do a flux calibration. The available options are \"gnirs\" and \"none\"."
        )
        mergeMethod = getParam(
        "Merge final data cubes? []: ",
        '',
        "Not yet implemented; Nifty has support for multiple merge methods."
        )

        rawPath = getParam(
        "Path to raw files directory? []: ",
        "",
        "An example of a valid raw files path string: \"/Users/nat/data/\""
        )
        program = getParam(
        "Gemini Program ID? []: ",
        "",
        "Nifty can either download raw data from the Gemini Public archive (to ./rawData) or copy it from " + \
        "a local directory. An example of a valid program string: \"GN-2013A-Q-62\". "
        )
        proprietaryCookie = getParam(
        "Cookie for proprietary downloads? []: ",
        '',
        "You can provide a cookie from your Gemini public archive login session to automatically " + \
        "download proprietary data."
        )
        dataSource = getParam(
        "Select a raw data source; 'GSA' for Gemini Science Archive, 'CADC' for Canadian Astronomy Data Centre. [GSA]: ",
        'GSA',
        "Automatic downloads can happen from either the Gemini Science Archive or the Canadian Astronomy Data Centre."
        )
        skyThreshold = getParam(
        "Sky threshold? [2.0]: ",
        2.0,
        "Nifty differentiates between sky and science frames from the telescope P and Q offsets. " + \
        "If sqrt(Poffset^2 + Qoffset^2) is more than the given threshold, a frame is marked as a science frame. " + \
        "Nifty also tries to take the telescope P and Q offset zero point into account if the first past seems to only identify sky frames."
        )
        sortTellurics = getParam(
        "Sort telluric observations? [yes]: ",
        'yes',
        "Nifty can sort telluric observations into the directory structure."
        )
        telluricTimeThreshold = getParam(
        "Max time between matched science and standard star frames? [5400]",
        5400,
        "Nifty will try to match science frames with the telluric frame closest to it in UT start time, " + \
        "within a certain threshold (in seconds). The default is 5400 seconds (1.5 hours)."
        )

        # By default do all of them.
        baselineCalibrationStart = getParam(
        "Starting point of baseline calibration reductions? [1]: ",
        1,
        "Specify the start step for the reduction of calibrations here."
        )
        baselineCalibrationStop = getParam(
        "Stopping point of baseline calibration reductions? [4]: ",
        4,
        "Specify the stop step for the reduction of calibrations here."
        )

        telStart = getParam(
        "Start point? [1]: ",
        1,
        "Starting point of science and telluric reductions."
        )
        telStop = getParam(
        "Stop point? [6]: ",
        6,
        "Stopping point of science and telluric reductions"
        )
        telluricSkySubtraction = getParam(
        "Do a telluric sky subtraction? [yes]: ",
        'yes',
        "Specify to subtract sky frames from telluric frames."
        )

        sciStart = getParam(
        "Starting point of science and telluric reductions? [1]: ",
        1,
        "Starting point of science reduction."
        )
        sciStop = getParam(
        "Stopping point of science and telluric reductions? [6]: ",
        6,
        "Stopping point of science reduction."
        )
        scienceSkySubtraction = getParam(
        "Subtract sky frames from science frames? [yes]: ",
        'yes',
        "Nifty can subtract a sky frame of the same exposure duration from each science frame."
        )

        telluricCorrectionStart = getParam(
        "Start step of telluric correction? [1]: ",
        1,
        "Choose a start step of the telluric correction."
        )
        telluricCorrectionStop = getParam(
        "Stop step of telluric correction? [5]: ",
        5,
        "Choose a stop step of the telluric correction."
        )
        hLineMethod = getParam(
        "H-line removal method? [vega]: ",
        "none",
        "Nifty can attempt to remove H-lines from a telluric correction spectrum. The available options are \"vega\" and \"none\"."
        )
        # Some of these are disabled (for now!) because of bugs in interactive Pyraf tasks.
        # TODO(nat): when interactive is fixed re-enable this.
        # Temp fix:
        hLineInter = getParam(
        "Interative H-line removal? [no]: ",
        False,
        "WARNING: This is currently broken due to bugs in interactive PyRAF tasks. Use with caution."
        )
        continuumInter = getParam(
        "Interative telluric continuum fitting? [no]: ",
        False,
        "WARNING: This is currently broken due to bugs in interactive PyRAF tasks. Use with caution."
        )
        telluricInter = getParam(
        "Interative telluric correction with iraf.telluric()? [no]: ",
        False,
        "WARNING: This is currently broken due to bugs in interactive PyRAF tasks. Use with caution."
        )
        tempInter = getParam(
        "Plot intermediate results using matplotlib? [no]: ",
        False,
        "As a short term fix, you can choose to plot intermediate results of the telluric correction at runtime "
        "using matplotlib."
        )
        standardStarSpecTemperature = getParam(
        "Effective temperature in kelvin of telluric standard star? [""]: ",
        "",
        "You can specify the temperature of the telluric standard star; if not Nifty will attempt " + \
        "a SIMBAD query to find Teff."
        )
        standardStarMagnitude = getParam(
        "Magnitude of standard star? [""]: ",
        "",
        "You can specify the magnitude of the telluric standard star. If not Nifty will attempt "+ \
        "a SIMBAD query to look it up."
        )
        standardStarRA = getParam(
        "RA of standard star? [""]: ",
        "",
        "You can specify the RA of the telluric standard star. If not Nifty will attempt "+ \
        "a SIMBAD query to look it up."
        )
        standardStarDec = getParam(
        "Dec of standard star? [""]: ",
        "",
        "You can specify the Dec of the telluric standard star. If not Nifty will attempt "+ \
        "a SIMBAD query to look it up."
        )
        standardStarBand = getParam(
        "Band of standard star? [""]: ",
        "",
        "You can specify the spectral band of the telluric standard star. If not Nifty will attempt "+ \
        "to look it up from the directory structure."
        )

        fluxCalibrationStart = getParam(
        "Start step of flux calibration? [1]: ",
        1,
        "Choose a start point for the flux calibration."
        )
        fluxCalibrationStop = getParam(
        "Stop step of flux calibration? [6]",
        6,
        "Choose a stop point for the flux calibration."
        )

        mergeStart = getParam(
        "Start step of the final cube merging? [1]: ",
        1,
        "Choose a start point for the final cube merging."
        )
        mergeStop = getParam(
        "Stop step of final cube merging? [3]:",
        6,
        "Choose a stop point for the final cube merging."
        )
        mergeType = getParam(
        "Type of merging to do? [median]: ",
        "median",
        "What type of merging would you like imcombine to do in merging? The options are the same as those "
        "available in imcombine; median, average, sum, etc."
        )
        use_pq_offsets = getParam(
        "Use pq offsets to merge data cubes? [yes]: ",
        "yes",
        "Nifty can merge cubes blindly using telescope P and Q offsets. If not, Nifty will pause and "+ \
        "ask you to first shift cubes by hand (say, in QFitsView) before merging cubes."
        )
        im3dtran = getParam(
        "Transpose cubes for faster merging? [yes]: ",
        'yes',
        "Nifty can transpose cubes to work around a bug in iraf.imcombine(). If not using this, note Nifty " + \
        "will take over 25 minutes to merge each cube."
        )

        # Save the options as a .cfg file.
        config = ConfigObj(RECIPES_PATH+'defaultConfig.cfg', unrepr=True)

        # General config used by all scripts.
        config['manualMode'] = manualMode
        config['over'] = over
        config['extractionXC'] = extractionXC
        config['extractionYC'] = extractionYC
        config['extractionRadius'] = extractionRadius
        config['scienceOneDExtraction'] = scienceOneDExtraction
        config['scienceDirectoryList'] = scienceDirectoryList
        config['calibrationDirectoryList'] = calibrationDirectoryList
        config['telluricDirectoryList'] = telluricDirectoryList

        config['nifsPipelineConfig'] = {}
        config['nifsPipelineConfig']['sort'] = sort
        config['nifsPipelineConfig']['calibrationReduction'] = calibrationReduction
        config['nifsPipelineConfig']['telluricReduction'] = telluricReduction
        config['nifsPipelineConfig']['scienceReduction'] = scienceReduction
        config['nifsPipelineConfig']['telluricCorrection'] = telluricCorrection
        config['nifsPipelineConfig']['fluxCalibration'] = fluxCalibration
        config['nifsPipelineConfig']['merge'] = merge
        config['nifsPipelineConfig']['telluricCorrectionMethod'] = telluricCorrectionMethod
        config['nifsPipelineConfig']['fluxCalibrationMethod'] = fluxCalibrationMethod
        config['nifsPipelineConfig']['mergeMethod'] = mergeMethod

        config['sortConfig'] = {}
        config['sortConfig']['rawPath'] = rawPath
        config['sortConfig']['program'] = program
        config['sortConfig']['proprietaryCookie'] = proprietaryCookie
        config['sortConfig']['dataSource'] = dataSource
        config['sortConfig']['skyThreshold'] = skyThreshold
        config['sortConfig']['sortTellurics'] = sortTellurics
        config['sortConfig']['telluricTimeThreshold'] = telluricTimeThreshold

        config['calibrationReductionConfig'] = {}
        config['calibrationReductionConfig']['baselineCalibrationStart']= baselineCalibrationStart
        config['calibrationReductionConfig']['baselineCalibrationStop'] = baselineCalibrationStop

        config['telluricReductionConfig'] = {}
        config['telluricReductionConfig']['telStart'] = telStart
        config['telluricReductionConfig']['telStop'] = telStop
        config['telluricReductionConfig']['telluricSkySubtraction'] = telluricSkySubtraction

        config['scienceReductionConfig'] = {}
        config['scienceReductionConfig']['sciStart'] = sciStart
        config['scienceReductionConfig']['sciStop'] = sciStop
        config['scienceReductionConfig']['scienceSkySubtraction'] = scienceSkySubtraction

        config['telluricCorrectionConfig'] = {}
        config['telluricCorrectionConfig']['telluricCorrectionStart'] = telluricCorrectionStart
        config['telluricCorrectionConfig']['telluricCorrectionStop'] = telluricCorrectionStop
        config['telluricCorrectionConfig']['hLineMethod'] = hLineMethod
        config['telluricCorrectionConfig']['hLineInter'] = hLineInter
        config['telluricCorrectionConfig']['continuumInter'] = continuumInter
        config['telluricCorrectionConfig']['telluricInter'] = telluricInter
        config['telluricCorrectionConfig']['tempInter'] = tempInter
        config['telluricCorrectionConfig']['standardStarSpecTemperature'] = standardStarSpecTemperature
        config['telluricCorrectionConfig']['standardStarMagnitude'] = standardStarMagnitude
        config['telluricCorrectionConfig']['standardStarRA'] = standardStarRA
        config['telluricCorrectionConfig']['standardStarDec'] = standardStarDec
        config['telluricCorrectionConfig']['standardStarBand'] = standardStarBand

        config['fluxCalbrationConfig'] = {}
        config['fluxCalbrationConfig']['fluxCalibrationStart'] = fluxCalibrationStart
        config['fluxCalbrationConfig']['fluxCalibrationStop'] = fluxCalibrationStop

        config['mergeConfig'] = {}
        config['mergeConfig']['mergeStart'] = mergeStart
        config['mergeConfig']['mergeStop'] = mergeStop
        config['mergeConfig']['mergeType'] = mergeType
        config['mergeConfig']['use_pq_offsets'] = use_pq_offsets
        config['mergeConfig']['im3dtran'] = im3dtran

        # Convert yes/no responses to True/False
        def update(u):
            for k, v in u.iteritems():
                if isinstance(v, collections.Mapping):
                    u[k] = update(u.get(k))
                else:
                    if u[k] == 'yes':
                        u[k] = True
                    elif u[k] == 'no':
                        u[k] = False
            return u

        update(config)


        with open('./config.cfg', 'w') as outfile:
            config.write(outfile)

    return fullReduction

#-----------------------------------------------------------------------------#

def datefmt():
    datefmt = '%Y/%m/%d %H:%M:%S '
    return datefmt

#-----------------------------------------------------------------------------#

def copyCalibration(inputFile, outputFile, grating, over):
    """
    Copy calibrations over to science directories.

    Copies inputFile to outputFile in all ../grating/scienceObservation/obs*/calibrations/ directories.
    """
    # Also copy it over to the relevant science directories.
    for scienceDirectory in glob.glob('../'+grating+'/obs*'):
        # Make sure the cals directory exists.
        if not os.path.exists(scienceDirectory+'/calibrations'):
            os.mkdir(scienceDirectory+'/calibrations')
        # Check if the file exists; if not, copy it over.
        if os.path.exists(scienceDirectory+'/calibrations/'+outputFile):
            if over:
                os.remove(scienceDirectory+'/calibrations/'+outputFile)
                shutil.copy(inputFile, scienceDirectory+'/calibrations/'+outputFile)
            else:
                logging.info("\nOutput exists and -over not set - skipping copy of shift file to science directory")
        else:
            shutil.copy(inputFile, scienceDirectory+'/calibrations/'+outputFile)
    # Also copy to telluric directories.
    for telluricDirectory in glob.glob('../'+grating+'/Tellurics/obs*'):
        # Make sure the cals directory exists.
        if not os.path.exists(telluricDirectory+'/calibrations'):
            os.mkdir(telluricDirectory+'/calibrations')
        # Check if the shift file exists; if not, copy it over.
        if os.path.exists(telluricDirectory+'/calibrations/'+outputFile):
            if over:
                os.remove(telluricDirectory+'/calibrations/'+outputFile)
                shutil.copy(inputFile, telluricDirectory+'/calibrations/'+outputFile)
            else:
                logging.info("\nOutput exists and -over not set - skipping copy of shift file to science directory")
        else:
            shutil.copy(inputFile, telluricDirectory+'/calibrations/'+outputFile)

def copyCalibrationDatabase(inputPrefix, grating, over):
    """
    Copy calibrations over to science directories.

    Copies inputPrefix* to all ../grating/scienceObservation/obs*/calibrations/database directories,
    under a new name, inputPrefix+fileType+_SCI_+sliceNumber+_
    """
    # Also copy it over to the relevant science directories.
    for scienceDirectory in glob.glob('../'+grating+'/obs*'):
        # Make sure the cals directory exists.
        if not os.path.exists(scienceDirectory+'/calibrations'):
            os.mkdir(scienceDirectory+'/calibrations')
        # Make sure the database directory exists
        if os.path.isdir(scienceDirectory+"/calibrations/database"):
            if glob.glob(scienceDirectory+"/calibrations/database/"+inputPrefix+"*"):
                if over:
                    for item in glob.glob(scienceDirectory+"/calibrations/database/"+inputPrefix+"*"):
                        os.remove(item)
                    for item in glob.glob('database/'+inputPrefix+'*'):
                        shutil.copy(item, scienceDirectory+"/calibrations/database/")
                else:
                    print "\nOutput exists and -over not set - skipping copy of database directory"
            else:
                for item in glob.glob('database/'+inputPrefix+'*'):
                    shutil.copy(item, scienceDirectory+"/calibrations/database/")
        else:
            os.mkdir(scienceDirectory+"/calibrations/database")
            for item in glob.glob('database/'+inputPrefix+'*'):
                shutil.copy(item, scienceDirectory+"/calibrations/database/")
    # Also copy to telluric directories
    for telluricDirectory in glob.glob('../'+grating+'/Tellurics/obs*'):
        # Make sure the cals directory exists.
        if not os.path.exists(telluricDirectory+'/calibrations'):
            os.mkdir(telluricDirectory+'/calibrations')
        # Make sure the database directory exists
        if os.path.isdir(telluricDirectory+"/calibrations/database"):
            if glob.glob(telluricDirectory+"/calibrations/database/"+inputPrefix+"*"):
                if over:
                    for item in glob.glob(telluricDirectory+"/calibrations/database/"+inputPrefix+"*"):
                        os.remove(item)
                    for item in glob.glob('database/'+inputPrefix+'*'):
                        shutil.copy(item, telluricDirectory+"/calibrations/database/")
                else:
                    print "\nOutput exists and -over not set - skipping copy of database directory"
            else:
                for item in glob.glob('database/'+inputPrefix+'*'):
                    shutil.copy(item, telluricDirectory+"/calibrations/database/")
        else:
            os.mkdir(telluricDirectory+"/calibrations/database")
            for item in glob.glob('database/'+inputPrefix+'*'):
                shutil.copy(item, telluricDirectory+"/calibrations/database/")


#-----------------------------------------------------------------------------#

def replaceNameDatabaseFiles(inputFile, oldFileName, newFileName):
    """
    Replace old filenames in database text files with new file names.
    """
    with open(inputFile, "r") as f:
        lines = f.read()
    modifiedLines = re.sub(oldFileName, newFileName, lines)
    with open("temp.txt", "w") as f:
        f.write(modifiedLines)
    os.remove(inputFile)
    shutil.move("temp.txt", re.sub(oldFileName, newFileName, inputFile))


#-----------------------------------------------------------------------------#
def copyResultsToScience(inputFile, outputFile, over):
    """

    """
    # Also copy it over to the relevant science directories.
    for scienceDirectory in glob.glob('../../obs*'):
        # Make sure the cals directory exists.
        if not os.path.exists(scienceDirectory+'/products_telluric_corrected'):
            os.mkdir(scienceDirectory+'/products_telluric_corrected')
        # Check if the file exists; if not, copy it over.
        if os.path.exists(scienceDirectory+'/products_telluric_corrected/'+outputFile):
            if over:
                os.remove(scienceDirectory+'/products_telluric_corrected/'+outputFile)
                shutil.copy(inputFile, scienceDirectory+'/products_telluric_corrected/'+outputFile)
            else:
                logging.info("\nOutput exists and -over not set - skipping copy of file to science directory")
        else:
            shutil.copy(inputFile, scienceDirectory+'/products_telluric_corrected/'+outputFile)

#-----------------------------------------------------------------------------#

def rewriteSciImageList(threshold, kind):
    """
    Find zero point from first image in skyFrameList
    calculate difference; if larger than threshold, append to new skyFrameList.
    Else append to new scienceFrameList.
    Write out both lists at end.
    """
    skyFrameList = open('skyFrameList', "r").readlines()
    skyFrameList = [image.strip() for image in skyFrameList]
    firstImage = astropy.io.fits.open(skyFrameList[0]+'.fits')
    p0 = firstImage[0].header['POFFSET']
    q0 = firstImage[0].header['QOFFSET']
    r0 = math.sqrt(((p0)**2) + ((q0)**2))
    newSkyFrameList = []
    newScienceFrameList = []
    for item in skyFrameList:
        image = astropy.io.fits.open(item+'.fits')
        poff = image[0].header['POFFSET']
        qoff = image[0].header['QOFFSET']
        pdiff = abs(p0 - poff)
        qdiff = abs(q0 - qoff)
        rdiff = math.sqrt(((pdiff)**2) + ((qdiff)**2))
        if rdiff >= threshold:
            newSkyFrameList.append(item)
        else:
            newScienceFrameList.append(item)
    os.remove("skyFrameList")
    for item in newSkyFrameList:
        writeList(item, 'skyFrameList', './')
        # Create a skyFrameList in the relevant directory.
    if kind == "Science":
        for item in newScienceFrameList:
            writeList(item, 'scienceFrameList', './')
    elif kind == "Telluric":
        for item in newScienceFrameList:
            writeList(item, 'tellist', './')

#-----------------------------------------------------------------------------#

def printDirectoryLists():
    """Print paths to science, telluric and calibration observations.

    Useful for:
        - Making sure sorting worked properly.
        - Making sure pipeline is loading runtimeData/scienceDirectoryList.txt,
          runtimeData/telluricDirectoryList.txt and runtimeData/calibrationDirectoryList.txt
          correctly.
    """
    # Print the current directory of data being reduced.
    logging.info("\nThe following values were found in ./config.cfg: ")

    with open('./config.cfg') as config_file:
        options = ConfigObj(config_file, unrepr=True)
    logging.info("\nScience Directory List: ")
    for i in range(len(options['scienceDirectoryList'])):
        logging.info(options['scienceDirectoryList'][i])
    logging.info("\nTelluric Directory List: ")
    for i in range(len(options['telluricDirectoryList'])):
        logging.info(options['telluricDirectoryList'][i])
    logging.info("\nCalibration Directory List: ")
    for i in range(len(options['calibrationDirectoryList'])):
        logging.info(options['calibrationDirectoryList'][i])

#-----------------------------------------------------------------------------#

def getParam(prompt, default, helpMessage="No help implemented yet!"):
    """Get a parameter from the user interactively. Also sets a default and
    displays a help message if a user types "h" or "help".
    """

    print "\n" + bcolors.OKBLUE + helpMessage + bcolors.ENDC + "\n"
    param = raw_input(prompt)
    if param == "no" or "No" or "N" or "n":
        param == False
    if param == "yes" or "Yes" or "Y" or "n":
        parm = True
    param = param or default

    return param

#-----------------------------------------------------------------------------#

def getFitsHeader(fitsFile, fitsKeyWords):
    """ imported from /astro/sos/da/scisoft/das/daLog/MakeDaDataCheckLogDefs.py """
    selection2 ="fullheader/"+fitsFile
    url2 ="http://fits/" + selection2
    u2 = urllib.urlopen(url2)
    xml2 = u2.read()
    u2.close()
    fitsHeaderList = [fitsFile[:-5]]
    for entry in fitsKeyWords:
        myOut = FitsKeyEntry(entry, xml2)
        fitsHeaderList.append(myOut)
    #
    return fitsHeaderList

#-----------------------------------------------------------------------------#

def FitsKeyEntry(fitsKeyWd, fullheader):
    """ imported from /astro/sos/da/scisoft/das/daLog/MakeDaDataCheckLogDefs.py """
    selectEntry ="none found"
    fullList = fullheader.splitlines()
    checkKeyWd = fitsKeyWd.ljust(8,' ')
    for index in range(len(fullList)):
        if fullList[index][:8] == checkKeyWd:
            if fullList[index][10] == "'":
                selectEntry = stripString(fullList[index])
            else:
                selectEntry = stripNumber(fullList[index])
    return selectEntry

#-----------------------------------------------------------------------------#

def stripString(inputString):
    """ imported from /astro/sos/da/scisoft/das/daLog/MakeDaDataCheckLogDefs.py """

    delimString ="'"
    delimList = []
    for index in range(len(inputString)):
        if inputString[index] == delimString:
            delimList.append(index)
    outFull = inputString[delimList[0]+1:delimList[-1]]
    outPut = outFull.replace(" ","")
    #
    return outPut

#-----------------------------------------------------------------------------#

def stripNumber(inputString):
    """ imported from /astro/sos/da/scisoft/das/daLog/MakeDaDataCheckLogDefs.py
    """
    delim1 ="="
    delim2 ="/"
    delimList = []
    for index in range(len(inputString)):
        if inputString[index] == delim1:
            delimList.append(index)
        if inputString[index] == delim2:
            delimList.append(index)
    if len(delimList) == 1:
        delimList.append(index)
    outFull = inputString[delimList[0]+1:delimList[1]]
    outPut = float(outFull)
    #
    return outPut

#-----------------------------------------------------------------------------#

def getUrlFiles(url,tag):
    """ imported from IPM scripts globaldefs.py """
    u = urllib.urlopen(url)
    xml = u.read()
    u.close()
    dom = parseString(xml)

    # Get file list:
    fileList = []
    previousFilename =""
    for fe in dom.getElementsByTagName(tag):
        fitsFile = str(fe.getElementsByTagName('filename')[0].childNodes[0].data)
        # exclude consecutive duplicates:
        if fitsFile != previousFilename:
            fileList.append(fitsFile)
            previousFilename = fitsFile

    #Return file list:
    return fileList

#-----------------------------------------------------------------------------#

def checkOverCopy(filelist, path, over):
    """ checks if over is True or False and copies files from /net/mko-nfs/sci/dataflo
    based on this.
    """

    rawfiles = []
    missingRaw = []

    raw = '/net/mko-nfs/sci/dataflo'


    for entry in filelist:
        if glob.glob(path+'/'+entry):
            rawfiles.append(glob.glob(path+'/'+entry))
        else:
            missingRaw.append(entry)

    if rawfiles:
        if over:
            for entry in rawfiles:
                if os.path.exists(entry[0]):
                    os.remove(entry[0])
            # copy all science images from a given night into ./Raw/
            for entry in filelist:
                if os.path.exists(raw+'/'+entry):
                    shutil.copy(raw+'/'+entry, path)
                else:
                    logging.info('SKIPPED ', entry)
        else:
            for entry in missingRaw:
                if os.path.exists(raw+'/'+entry):
                    shutil.copy(raw+'/'+entry, path)
                else:
                    logging.info('SKIPPED ', entry)

    else:
        for entry in filelist:
            if os.path.exists(raw+'/'+entry):
                shutil.copy(raw+'/'+entry, path)
            else:
                logging.info('SKIPPED ', entry)

    return

#-----------------------------------------------------------------------------#

def checkQAPIreq(alist):
    """ checks to make sure that the arcs meet the PI and QA requirements """

    blist = []
    for entry in alist:
        blist.append(entry)
    for i in range(len(alist)):
        fitsKeyWords = ['RAWPIREQ', 'RAWGEMQA']
        headerList = getFitsHeader(alist[i], fitsKeyWords)
        rawPIreq = headerList[1]
        rawGemQA = headerList[2]
        if rawPIreq in ["YES","UNKNOWN"] and rawGemQA in ["USABLE","UNKNOWN"]:
            logging.info(alist[i]+' added for processing')
        else:
            logging.info(alist[i]+' excluded, set to USABLE/FAIL')
            blist.remove(alist[i])

    return blist

#-----------------------------------------------------------------------------#

def listit(list, prefix):
    """ Returns a string where each element of list is prepended with prefix """

    l = []
    for x in list:
        l.append(prefix+(x.strip()).rstrip('.fits'))
    return ",".join(l)

#-----------------------------------------------------------------------------#

def checkDate(list):
    """ check the dates on all the telluric and acquisition files to make sure that
        there are science images on the same night
    """

    removelist = []
    datelist = []

    for entry in list:
        fitsKeyWords = ['DATE', 'OBSCLASS']
        header = getFitsHeader(entry, fitsKeyWords)
        date = header[1]
        obsclass = header[2]
        if obsclass=='science':
            if not datelist or not datelist[-1]==date:
                datelist.append(date)

    for entry in list:
        fitsKeyWords = ['DATE', 'OBSCLASS']
        header = getFitsHeader(entry, fitsKeyWords)
        date = header[1]
        obsclass = header[2]
        if not obsclass=='science' and date not in datelist:
            removelist.append(entry)

    return removelist

#-----------------------------------------------------------------------------#

def writeList(image, file, path):
    """ write image name into a file """
    homepath = os.getcwd()
    os.chdir(path)
    image = image.rstrip('.fits')
    if os.path.exists(file):
        with open(file, 'r') as f:
            filelist = f.readlines()
        if image+'\n' in filelist:
            pass
        else:
            with open(file, 'a') as f:
                f.write(image+'\n')
    else:
        with open(file, 'a') as f:
            f.write(image+'\n')
    os.chdir(homepath)

#-----------------------------------------------------------------------------#

def checkEntry(entry, entryType, filelist):
    """ checks to see the that program ID given matches the OBSID in the science headers
        checks to see that the date given matches the date in the science headers
    """

    if entryType == 'program':
        header = getFitsHeader(filelist[0], ['OBSID'])
        if entry in header[1]:
            pass
        else:
            logging.info("\n Program number was entered incorrectly.\n")
            raise SystemExit

    if entryType == 'date':
        header = getFitsHeader(filelist[0], ['DATE'])
        if entry in header[1].replace('-',''):
            pass
        else:
            logging.info("\n Date was entered incorrectly or there is no NIFS data for the date given. Please make sure the date has been entered as such: YYYYDDMM.\n")
            raise SystemExit

#-----------------------------------------------------------------------------#

def checkLists(original_list, path, prefix, suffix, outlist=None):
    """Check that all files made it through an iraf step. """

    new_list = []
    for image in original_list:
        image = image.strip()
        if os.path.exists(path+'/'+prefix+image+suffix):
            new_list.append(image)
        else:
            logging.info('\n' +  str(image)+ '.fits not being processed due to error in image.\n')
            logging.info("\n#####################################################################")
            logging.info("#####################################################################")
            logging.info("")
            logging.info("     WARNING: " + str(image) + " .fits was removed from a list after a checkLists call.")
            logging.info("               An iraf task may have failed. ")
            logging.info("")
            logging.info("#####################################################################")
            logging.info("#####################################################################\n")


            pass

    if outlist:
        # Write outlist in path
        with open(os.path.join(path, outlist), 'w') as f:
            for frame in new_list:
                f.write(frame+"\n")

    return new_list

#-----------------------------------------------------------------------------#

def checkSameLengthFlatLists():
    """Reads two textfile lists of filenames. If not the same length,
    removes last entry from longer list until they are. Prints loud warnings to
    tell people it is doing this.
    This is done because it was found nsflat() complains when the gain of Combined
    lamps-on flats does not match the gain of combined lamps-off flats.
    This is caused when different numbers of flats combined into one combined frame.
    It seemed simplest to exclude one of the darks. Feel free to attempt a more complicated
    fix!
    """
    # Read the flatlist into a python list.
    with open('./flatlist', 'r') as f:
        flatlist = f.readlines()
    # Read the flatdarklist into a python list.
    with open('./flatdarklist', 'r') as f:
        flatdarklist = f.readlines()
    # Check that both lists are the same length.
    if len(flatlist) != len(flatdarklist):
        # Print a nice loud warning.
        logging.info("\n#####################################################################")
        logging.info("#####################################################################")
        logging.info("")
        logging.info("     WARNING in sort: flatlist and flatdarklist are not the same ")
        logging.info("                      length. Removing extra entries from the")
        logging.info("                      longer list. Original lists can be found in")
        logging.info("                      original_flatlist and original_flatdarklist")
        logging.info("     in " + str(os.getcwd()))
        logging.info("")
        logging.info("#####################################################################")
        logging.info("#####################################################################\n")
        # Copy the original flatlist and flatdarklist to backup files.
        shutil.copy2('./flatlist', './original_flatlist')
        shutil.copy2('./flatdarklist', './original_flatdarklist')
        # while they are not the same length:
        while len(flatlist) != len(flatdarklist):
            # remove the last entry from the longer list.
            if len(flatlist) > len(flatdarklist):
                del flatlist[-1]
            else:
                del flatdarklist[-1]
            if len(flatlist) == 0 or len(flatdarklist) == 0:
                raise CalibrationsError("Error trying to make a flatlist and a flatdarklist the same length in {}.".format(os.getcwd()))
        # Write the new flatlist to the flatlist textfile, overwriting anything already there.
        with open('./flatlist', 'w') as f:
            for item in flatlist:
                f.write(item)
        # Write the new flatdarklist to the flatdarklist textfile, overwriting anything already there.
        with open('./flatdarklist', 'w') as f:
            for item in flatdarklist:
                f.write(item)

#-----------------------------------------------------------------------------#

def makeSkyList(skyFrameList, sciencelist, obsDir):
    """ Makes a skyFrameList equal in length to object list with sky and object
    frames closest in time at equal indices.

    Writes the original skyFrameList to original_skyFrameList.

    Writes results to skyFrameList textfile.

    Returns:
        skyFrameList (list): list of sky frames organized so each science frame has subtracted
                        the sky frame closest in time.

    Eg:
        Observations had an ABA ABA pattern:
            obs1
            sky1
            obs2

            obs3
            sky2
            obs4

            obs5
            sky3
            obs6

        sciencelist was:    skyFrameList was:   Output skyFrameList will be:
            obs1                    sky1            sky1
            obs2                    sky2            sky1
            obs3                    sky3            sky2
            obs4                                    sky2
            obs5                                    sky3
            obs6                                    sky3

        Observations had an AB AB AB pattern:
            obs1
            sky1
            obs2
            sky2
            obs3
            sky3

        sciencelist was:    skyFrameList was:   Output skyFrameList will be:
            obs1                    sky1            sky1
            obs2                    sky2            sky1
            obs3                    sky3            sky2
    """
    logging.info("\n#############################################################")
    logging.info("#                                                           #")
    logging.info("#  Matching science frames with sky frames closest in time  #")
    logging.info("#                                                           #")
    logging.info("#############################################################\n")
    # Do some tests first.
    # Check that data is either:
    # ABA ABA ABA- one sky frame per two science frames.
    # AB AB AB- one sky frame per one two science frames.
    #
    # If it is neither warn user to verify that sky frames were matched with science frames correctly.
    if len(skyFrameList) != len(sciencelist)/2 and len(skyFrameList) != len(sciencelist):
        logging.info("\n#####################################################################")
        logging.info("#####################################################################")
        logging.info("")
        logging.info("     WARNING in reduce: it appears science frames and sky frames were not")
        logging.info("                        taken in an ABA ABA or AB AB pattern.")
        logging.info("")
        logging.info("#####################################################################")
        logging.info("#####################################################################\n")
    skytimes = []
    prepared_sky_list = []
    # Calculate time of each sky frame. Store the calculated time and the frame name in skytimes, a
    # 2D list of [skyframe_time, skyframe_name] pairs.
    # Eg: [[39049.3, 'N20130527S0241'], [39144.3, 'N20130527S0244'], [39328.8, 'N20130527S0247'], [39590.3, 'N20130527S0250']]
    for item in skyFrameList:
        # Strip off the trailing newline.
        item = str(item).strip()
        # Calculate the time of the sky frame.
        skytime = timeCalc(item+'.fits')
        # Store the sky frame time and corresponding sky frame name in skytimes.
        templist = [skytime, item]
        skytimes.append(templist)
    logging.info("scienceframelist:      skyFrameList:      time delta (between observation UT start times from .fits headers):")
    for item in sciencelist:
        # Calculate time of the science frame in seconds.
        item = str(item).strip()
        sciencetime = timeCalc(item+'.fits')
        # Sort the 2D list of [skyframe_time, skyframe_name] pairs by absolute science_frame_time - skyframe_time.
        # Eg: [[39049.3, 'N20130527S0241'], [39144.3, 'N20130527S0244'], [39328.8, 'N20130527S0247'], [39590.3, 'N20130527S0250']]
        sorted_by_closest_time = sorted(skytimes, key=lambda x: (abs(sciencetime - x[0])))
        # Append the name corresponding to the minimum time difference to prepared_sky_list.
        prepared_sky_list.append(sorted_by_closest_time[0][1])
        # Print the scienceframe, matching skyframe and time difference side by side for later comparison.
        logging.info("  "+ str(item)+ "       "+ str(sorted_by_closest_time[0][1])+ "        "+ str(abs(sciencetime - sorted_by_closest_time[0][0])))
    logging.info("\n")

    os.rename('skyFrameList', 'original_skyFrameList')

    f = open('skyFrameList', 'w')
    for image in prepared_sky_list:
        f.write(image+'\n')
    f.close()

    return prepared_sky_list

#-----------------------------------------------------------------------------#


def convertRAdec(ra, dec):
    """ converts RA from degrees to H:M:S and dec from degrees to degrees:arcmin:arcsec"""
    H = int(ra/15.)
    M = int((ra-(15*H))/.25)
    S = ((ra-(float(H)*15.))-(float(M)*.25))/(1./240.)

    ra = str(H)+'h'+str(M)+'m'+str(S)+'s'

    return ra

#-----------------------------------------------------------------------------#

def timeCalc(image):
    """Read time from .fits header. Convert to a float of seconds.
    """
    telheader = astropy.io.fits.open(image)
    UT = telheader[0].header['UT']
    secs = float(UT[6:10])
    mins = float(UT[3:5])
    hours = float(UT[0:2])
    time = secs+mins*60.+hours*(60.*60.)

    return time

#-----------------------------------------------------------------------------#

def MEFarithpy(MEF, image, op, result):

    if os.path.exists(result):
        os.remove(result)
    scimage = astropy.io.fits.open(MEF+'.fits')
    arithim = astropy.io.fits.open(image+'.fits')
    for i in range(88):
        if scimage[i].name=='SCI':
            if op=='multiply':
                scimage[i]=scimage[i]*arithim
            if op=='divide':
                scimage[i]=scimage[i]/arithim
    scimage.writeto(result, output_verify='ignore')
#-----------------------------------------------------------------------------#

def MEFarith(MEF, image, op, result):

    if os.path.exists(result):
        os.remove(result)
    iraf.fxcopy(input=MEF+'[0]', output=result)
    for i in range(1,88):
        iraf.fxinsert(input=MEF+'['+str(i)+']', output=result+'['+str(i)+']', groups='', verbose = 'no')
    for i in range(1,88):
        header = astropy.io.fits.open(result)
        extname = header[i].header['EXTNAME']
        if extname == 'SCI':
            iraf.imarith(operand1=result+'['+str(i)+']', op=op, operand2 = image, result = result+'['+str(i)+', overwrite]', divzero = 0.0)

#-----------------------------------------------------------------------------#

def downloadQueryCadc(query, directory='./rawData'):
    """
    Finds and downloads all CADC files for a particular gemini program ID to
    the current working directory.
    """

    cadc = Cadc()
    job = cadc.create_async(query)
    job.run().wait()
    job.raise_if_error()
    result = job.fetch_result().to_table()

    # Store product id's for later
    pids = list(result['productID'])

    urls = cadc.get_data_urls(result)
    cwd = os.getcwd()
    os.chdir(directory)
    for url, pid in zip(urls, pids):
        try:
            filename = getFile(url)
            logging.debug("Downloaded {}".format(filename))
        except Exception as e:
            logging.error("A frame failed to download.")
            os.chdir(cwd)
            raise e
    os.chdir(cwd)


def getFile(url):
    """
    Gets a file from the specified url and returns the filename.
    """
    r = requests.get(url, stream=True)
    # Parse out filename from header
    try:
        filename = re.findall("filename=(.+)", r.headers['Content-Disposition'])[0]
    except KeyError:
        # 'Content-Disposition' header wasn't found, so parse filename from URL
        # Typical URL looks like:
        # https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/data/pub/GEM/N20140505S0114.fits?RUNID=mf731ukqsipqpdgk
        filename = (url.split('/')[-1]).split('?')[0]
    
    # Write the fits file to the current directory, verifying the md5 hash as we go. Store partial results in a temporary file.
    writeWithTempFile(r, filename)

    return filename


def writeWithTempFile(request, filename):
    """ Write the fits file, verifying the md5 hash as we go. Store partial results in a temporary file. """
    temp_downloads_path = '.temp-downloads'
    if not os.path.exists(temp_downloads_path):
        os.mkdir(temp_downloads_path)
    try:
        server_checksum = request.headers['Content-MD5']
    except KeyError:
        # Catch case that header didn't contain a 'content-md5' header
        logging.warning("Content-MD5 header not found for file {}. Skipping checksum validation.".format(filename))
        server_checksum = None

    # Write out content (first to a temp file) optionally doing an md5 verification.
    download_checksum = hashlib.md5()
    with tempfile.TemporaryFile(mode='w+b', prefix=filename, dir=temp_downloads_path) as f:
        for chunk in request.iter_content(chunk_size=128):
            f.write(chunk)
            download_checksum.update(chunk)
        if server_checksum and (server_checksum != download_checksum.hexdigest()):
            logging.error("Problem downloading {} from {}.".format(filename, url))
            raise IOError
        f.seek(0)
        with open(filename, 'w') as out_fp:
            out_fp.write(f.read())

    return filename


def checkForMDFiles(path, query):
    """
    Gemini sometimes attaches an "-md!" prefix to their files to indicate there's
    some sort of metadata problem with them. The CADC archive doesn't currently
    contain that information, but it can be queried from the Gemini Science Archive.
    """
    try:
        # This abuses getFile because (currently) if the content-disposition header isn't found it
        # uses a hard-coded part of the url as the filename (which makes no sense here). However,
        # it doesn't seem as if that should break things.
        json_filename = getFile(query)
    except Exception:
        logging.warning("MD Report check failed. Some fits files with broken metadata may have made it through.")
        return

    with open(json_filename) as f:
        data = json.load(f)
    for item in data:
        try:
            if not item['mdready']:
                # Do some sort of removal of the actual fits file here
                fits_filename = item['filename'].rstrip('.bz2')
                if os.path.exists(os.path.join(path, fits_filename)):
                    os.remove(os.path.join(path, fits_filename))
        except Exception:
            continue

class CalibrationsError(Exception):
    """General error for tasks while getting calibrations."""
    pass

class CalibrationsNotFoundError(CalibrationsError):
    """Raised when calibrations aren't found for a particular science frame."""
    pass

class TelluricsNotFoundError(Exception):
    """Raised when tellurics aren't found for a particular science frame."""
    pass

class ScienceObservationError(Exception):
    """Raised when there's a problem with a particular science directory.
    Can be used to skip the reduction for that particular directory."""
    pass

class ObservationDirError(Exception):
    """Raised when there's a problem with a particular observation directory,
    (science OR telluric) and that directory should be excluded from the data reduction process."""

class SkyFrameError(ObservationDirError):
    """Raised when there's a problem with the skyFrameList in a particular directory."""

class WavelengthError(Exception):
    """Raised when a non-standard wavelength is detected in the pipeline."""

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
            self.date = header[0].header['DATE-OBS'].replace('-','')
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


class ProductTagger:
    def __init__(self, configFile):
        user_clobber=iraf.envget("clobber")
        iraf.reset(clobber='yes')
        self.scienceDirectories = self.parseConfig(configFile)


    def run(self):
        """
        Tags products (and updates datalabels) in science directories with an extra extension that identifies calibrations.
        """
        for scienceDir in self.scienceDirectories:  
            try:
                cal_ext = self.makeCalExt(scienceDir)
            except Exception:
                print("Failed to make calibration extension for {}. Skipping.".format(scienceDir))
                continue

            try:
                bpm_file = self.getBPMFile(scienceDir)
            except Exception:
                print("Failed to make bad pixel mask extension for {}. It will not be included.".format(scienceDir))
                bpm_file = None

            
            for productType in ["uncorrected", "telluric_corrected", "fluxcal_AND_telluric_corrected"]:
                try:
                    self.tagProducts(scienceDir, productType, cal_ext, bpm_file)
                except Exception:
                    print("Failed to tag {} products in {}. Skipping.".format(productType, scienceDir))


    def makeCalExt(self, scienceDir):
        shift = glob.glob(os.path.join(scienceDir, "calibrations", "*_shift.fits"))
        flat = glob.glob(os.path.join(scienceDir, "calibrations", "*_flat.fits"))
        arc = glob.glob(os.path.join(scienceDir, "calibrations", "*_arc.fits"))
        ronchi = glob.glob(os.path.join(scienceDir, "calibrations", "*_ronchi.fits"))

        cals = [shift, flat, arc, ronchi]

        if not all([len(x) <= 1 for x in cals]):
            logging.error("A calibration file wasn't found in {}.".format(scienceDir))
            raise IOError()

        cals = [os.path.split(x[0])[-1] if x else "" for x in cals]

        a1 = np.array(['Processed Shift File', 'Processed Flat', 'Processed Arc', 'Processed Ronchi'])
        a2 = np.array(cals)
        col1 = astropy.io.fits.Column(name='Calibration Type', format='50A', array=a1)
        col2 = astropy.io.fits.Column(name='Filename', format='50A', array=a2)

        cols = astropy.io.fits.ColDefs([col1, col2])

        hdu = astropy.io.fits.BinTableHDU.from_columns(cols, name="CAL")

        return hdu


    def getBPMFile(self, scienceDir):
        bpm = glob.glob(os.path.join(scienceDir, "calibrations", "*.pl"))
        return bpm[0]



    def tagProducts(self, scienceDir, productType, cal_ext, bpm_file):
        if productType == "uncorrected":
            prefix = "ctfbrsn"
            products = glob.glob(os.path.join(scienceDir, "products_"+productType, prefix+"N*"))
        elif productType == "telluric_corrected":
            prefix = "actfbrsn"
            products = glob.glob(os.path.join(scienceDir, "products_"+productType, prefix+"N*"))
        elif productType == "fluxcal_AND_telluric_corrected":
            prefix = "factfbrsn"
            products = glob.glob(os.path.join(scienceDir, "products_"+productType, prefix+"N*"))
        else:
            raise ValueError("Invalid product type: {}".format(productType))

        for product in products:
            try:
                with astropy.io.fits.open(product, mode='update') as hdu1:
                    hdu1['PRIMARY'].header['DATALAB'] = hdu1['PRIMARY'].header['DATALAB']+'-'+prefix.upper()
                    hdu1['PRIMARY'].header['FLATIMAG'] = hdu1['PRIMARY'].header['FLATIMAG'].split(os.path.sep)[1]+'.fits'
                    hdu1['PRIMARY'].header['BPMFILE'] = hdu1['PRIMARY'].header['BPMFILE'].split(os.path.sep)[1]+'.pl'
                    if not self.hasCalExt(hdu1):
                        hdu1.append(cal_ext)
                    hasBPMFlag = self.hasBPMExt(hdu1)
                    hdu1.flush()
                if not hasBPMFlag:
                    iraf.imcopy(bpm_file, product+"[BPM,type=mask,append]")
                    fits.delval(product, "DATALAB", extname='BPM')
                    fits.delval(product, "NSFLON1", extname='BPM')
                    fits.delval(product, "NSFLOF1", extname='BPM')

            except Exception:
                print("Problem adding cal and bpm extension to {}. Skipping.".format(product))
                continue


    def parseConfig(self, configFile):
        try:
            with open(configFile) as f:
                lines = f.readlines()
        except IOError:
            return []

        try:
            line = [line for line in lines if "scienceDirectoryList =" in line][0]
            line = line.rstrip('\n').replace("scienceDirectoryList = [", "")
            line = line.replace("]", "")
            line = line.replace("'", "")
            scienceDirectories = line.split(",")
        except Exception:
            return []

        return scienceDirectories

    def hasCalExt(self, hdulist):
        return any([x.name == "CAL" for x in hdulist])

    def hasBPMExt(self, hdulist):
        return any([x.name == "BPM" for x in hdulist])


class CalibrationTagger:
    def __init__(self, calDir):
        user_clobber=iraf.envget("clobber")
        iraf.reset(clobber='yes')
        self.calDir = calDir

        self.keywordDict = {
                        "FCOMBINE": "Num flat frames used in processing",
                        "ACOMBINE": "Num arc frames used in processing",
                        "DCOMBINE": "Num dark frames used in processing",
                        "FCMB":     "Flat frame used in processing",
                        "ACMB":     "Arc frame used in processing",
                        "DCMB":     "Dark frame used in processing",
                        "SHFILE":   "File supplying mdf shift in processing",
                        "FLTFILE":  "Processed flat used in processing"
        }

        self.extensionDescriptions = {
                        "INPUTFLAT": "Processed flat field frame.",
                        "INPUTDARK": "Raw dark frame.",
                        "MEMBERFLAT": "Raw flat field frame.",
                        "MEMBERRONCHI": "Raw ronchi spatial correction frame.",
                        "MEMBERARC": "Raw arc frame."
        }

        self.extensionColumnNames = {
                        "label": "Filename",
                        "description": "Description",
                        "type": "Type (member or input)"
        }

    def run(self):
        """
        Tags calibrations with header keywords and extensions.
        """

        class Cals:
            pass
        cals = Cals()

        try:
            cals.flats = self.parseList(os.path.join(self.calDir, 'flatlist'))
            cals.flatdarks = self.parseList(os.path.join(self.calDir, 'flatdarklist'))
            cals.arcs = self.parseList(os.path.join(self.calDir, 'arclist'))
            cals.arcdarks = self.parseList(os.path.join(self.calDir, 'arcdarklist'))
            cals.ronchis = self.parseList(os.path.join(self.calDir, 'ronchilist'))
            cals.flat_file = self.getCalFile(os.path.join(self.calDir, '*_flat.fits'))
            cals.shift_file = self.getCalFile(os.path.join(self.calDir, '*_shift.fits'))
            cals.arc_file = self.getCalFile(os.path.join(self.calDir, '*_arc.fits'))
            cals.ronchi_file = self.getCalFile(os.path.join(self.calDir, '*_ronchi.fits'))
            cals.bpm_file = glob.glob(os.path.join(self.calDir, '*_bpm.pl'))[0]
        except Exception:
            logging.error("Problem tagging cals in {}. Skipping.".format(self.calDir))
            return

        self.tagFlat(cals)
        self.tagArc(cals)
        self.tagRonchi(cals)
        self.tagShift(cals)

    def makeProvenance(self, labels, descriptions, types):
        if not (len(labels) == len(descriptions)) and (len(descriptions) == len(types)):
            logging.error("Error making provenance extension, lengths of labels/descriptions/types " \
                        "should have been equal and were {}, {}, {}.".format(len(labels), len(descriptions), len(types)))
            raise IOError()

        a1 = np.array(labels)
        a2 = np.array(descriptions)
        a3 = np.array(types)
        col1 = astropy.io.fits.Column(name=self.extensionColumnNames['label'], format='64A', array=a1)
        col2 = astropy.io.fits.Column(name=self.extensionColumnNames['description'], format='64A', array=a2)
        col3 = astropy.io.fits.Column(name=self.extensionColumnNames['type'], format='64A', array=a3)

        cols = astropy.io.fits.ColDefs([col1, col2, col3])

        hdu = astropy.io.fits.BinTableHDU.from_columns(cols, name="PROVENANCE")

        return hdu


    def tagFlat(self, cals):
        """ Flats contain:
            - raw flats as members
            - raw darks as inputs
        """

        try:
            filenames = cals.flats[:]
            filenames.extend(cals.flatdarks)
            filenames = [x.split(os.path.sep)[-1] for x in filenames] # Only want filename
            
            descriptions = [self.extensionDescriptions['MEMBERFLAT']]*len(cals.flats)
            descriptions.extend([self.extensionDescriptions['INPUTDARK']]*len(cals.flatdarks))
            
            types = ['member']*len(cals.flats)
            types.extend(['input']*len(cals.flatdarks))

            provenance_extension = self.makeProvenance(filenames, descriptions, types)
            
            with fits.open(cals.flat_file, mode="update") as hdul:
                # Put updates to processed flat headers here
                if len(cals.flats) > 1 and "FLAT" not in hdul['PRIMARY'].header['DATALAB']:
                    hdul['PRIMARY'].header['DATALAB'] += "-FLAT"

                hasBPMFlag = self.hasBPMExt(hdul)
                hdul.append(provenance_extension)
                hdul.flush()
            if not hasBPMFlag:
                iraf.imcopy(cals.bpm_file, cals.flat_file+"[BPM,type=mask,append]")
                fits.delval(cals.flat_file, "DATALAB", extname='BPM')
                fits.delval(cals.flat_file, "NSFLON1", extname='BPM')
                fits.delval(cals.flat_file, "NSFLOF1", extname='BPM')
        except Exception:
            logging.error("Problem tagging {}.".format(cals.flat_file))
            

    def tagArc(self, cals):
        """Arcs contain:
            - Raw arcs as members
            - Raw darks as inputs
            - Processed flat as input 
        """
        try:
            filenames = cals.arcs[:]
            filenames.extend(cals.arcdarks)
            filenames.append(cals.flat_file)
            filenames = [x.split(os.path.sep)[-1] for x in filenames] # Only want filename
            
            descriptions = [self.extensionDescriptions['MEMBERARC']]*len(cals.arcs)
            descriptions.extend([self.extensionDescriptions['INPUTDARK']]*len(cals.arcdarks))
            descriptions.append(self.extensionDescriptions['INPUTFLAT'])
            
            types = ['member']*len(cals.arcs)
            types.extend(['input']*len(cals.arcdarks))
            types.append('input')

            provenance_extension = self.makeProvenance(filenames, descriptions, types)
            
            with fits.open(cals.arc_file, mode="update") as hdul:
                # Put updates to processed arc headers here
                if len(cals.arcs) > 1 and "ARC" not in hdul['PRIMARY'].header['DATALAB']:
                    hdul['PRIMARY'].header['DATALAB'] += "-ARC"

                #if ".fits" not in hdul['PRIMARY'].header['FLATIMAG']:
                #    hdul['PRIMARY'].header['FLATIMAG'] = hdul['PRIMARY'].header['FLATIMAG']+".fits"

                hasBPMFlag = self.hasBPMExt(hdul)
                hdul.append(provenance_extension)
                hdul.flush()
            if not hasBPMFlag:
                iraf.imcopy(cals.bpm_file, cals.arc_file+"[BPM,type=mask,append]")
                fits.delval(cals.arc_file, "DATALAB", extname='BPM')
                fits.delval(cals.arc_file, "NSFLON1", extname='BPM')
                fits.delval(cals.arc_file, "NSFLOF1", extname='BPM')
        except Exception:
            logging.error("Problem tagging {}.".format(cals.arc_file))    

    def tagRonchi(self, cals):
        """Ronchis contain:
            - Raw ronchis as members
            - Raw darks as inputs
            - Processed flat as input
        """
        try:
            filenames = cals.ronchis[:]
            filenames.extend(cals.flatdarks)
            filenames.append(cals.flat_file)
            filenames = [x.split(os.path.sep)[-1] for x in filenames] # Only want filename
            
            descriptions = [self.extensionDescriptions['MEMBERRONCHI']]*len(cals.ronchis)
            descriptions.extend([self.extensionDescriptions['INPUTDARK']]*len(cals.flatdarks))
            descriptions.append(self.extensionDescriptions['INPUTFLAT'])
            
            types = ['member']*len(cals.ronchis)
            types.extend(['input']*len(cals.flatdarks))
            types.append('input')

            provenance_extension = self.makeProvenance(filenames, descriptions, types)
            
            with fits.open(cals.ronchi_file, mode="update") as hdul:
                # Put updates to processed ronchi headers here
                if len(cals.ronchis) > 1 and "RONCHI" not in hdul['PRIMARY'].header['DATALAB']:
                    hdul['PRIMARY'].header['DATALAB'] += "-RONCHI"

                #if ".fits" not in hdul['PRIMARY'].header['FLATIMAG']:
                #    hdul['PRIMARY'].header['FLATIMAG'] = hdul['PRIMARY'].header['FLATIMAG']+".fits"
                #if ".fits" not in hdul['PRIMARY'].header['DARKIMAG']:
                #    hdul['PRIMARY'].header['DARKIMAG'] = hdul['PRIMARY'].header['DARKIMAG']+".fits"
                
                hasBPMFlag = self.hasBPMExt(hdul)
                hdul.append(provenance_extension)
                hdul.flush()
            if not hasBPMFlag:
                iraf.imcopy(cals.bpm_file, cals.ronchi_file+"[BPM,type=mask,append]")
                fits.delval(cals.ronchi_file, "DATALAB", extname='BPM')
                fits.delval(cals.ronchi_file, "NSFLON1", extname='BPM')
                fits.delval(cals.ronchi_file, "NSFLOF1", extname='BPM')
        except Exception:
            logging.error("Problem tagging {}.".format(cals.ronchi_file))


    def tagShift(self, cals):
        pass

    def hasBPMExt(self, hdulist):
        return any([x.name == "BPM" for x in hdulist])

    def parseList(self, listpath):
        """
        Tries to open a list in listpath, asserts it's non-empty, and well-formed.
        """
        with open(listpath) as f:
            lines = f.readlines()

        lines = [x.rstrip('\n') + '.fits' for x in lines]

        if (not len(lines) > 0) or not all(["N2" in x for x in lines]):
            logging.warning("List {} was either empty or did not have an 'N2' prefix in every line.".format(listpath))

        basePath = os.path.split(listpath)[0]

        lines = [os.path.join(basePath, x) for x in lines]

        return lines

    def getCalFile(self, calpath):
        try:
            calfile = glob.glob(calpath)[0]
        except IndexError as e:
            logging.error("Calibration file {} not found.".format(calpath))
            raise IOError()
        with fits.open(calfile) as hdul:
            if not len(hdul) > 0:
                logging.error("hdulist {} appears to have no header units.".format(calfile))
        return calfile


#-----------------------------------------------------------------------------#

