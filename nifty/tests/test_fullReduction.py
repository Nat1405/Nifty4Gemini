from __future__ import print_function
import os
import sys
import logging
import json
import glob
from pathlib import Path
import requests
import tarfile
import shutil
import pkg_resources
import astropy.io.fits as fits
from pyraf import iraf, iraffunctions

import pytest

from nifty.pipeline.configobj.configobj import ConfigObj
from common import data_path, download_test_data
import nifty.pipeline.nifsPipeline
from nifty.pipeline.nifsUtils import CalibrationTagger

try:
    from mock import Mock, patch, PropertyMock
except ImportError:
    pytest.skip("Install mock for the cadc tests.", allow_module_level=True)

PROVENANCE_EXT_NAME = CalibrationTagger.provenanceExtensionName

def test_quickGN2014AQ85_uncorrected(tmpdir, monkeypatch):
    """
        - Set up working directories
        - Download raw data
    """

    tmpdir = str(tmpdir)
    os.chdir(tmpdir)

    # Copy config file over and set it up for testing.
    RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
    shutil.copy(os.path.join(RECIPES_PATH, 'defaultConfig.cfg'), os.path.join(os.getcwd(), 'config.cfg'))

    with open('config.cfg') as config_file:
        config = ConfigObj(config_file, unrepr=True)
    
    config['ManualMode'] = False
    config['over'] = False
    
    config['nifsPipelineConfig']['sort'] = True
    config['nifsPipelineConfig']['calibrationReduction'] = True
    config['nifsPipelineConfig']['telluricReduction'] = False
    config['nifsPipelineConfig']['scienceReduction'] = True
    config['nifsPipelineConfig']['telluricCorrection'] = False
    config['nifsPipelineConfig']['fluxCalibration'] = False
    config['nifsPipelineConfig']['merge'] = False

    config['sortConfig']['rawPath'] = os.path.join(os.getcwd(), 'rawData')

    with open('config.cfg', 'w') as config_file:
        config.write(config_file)

    # Copy over test data
    if not os.path.exists(data_path('GN-2014A-Q-85-all')):
        download_test_data()
    os.mkdir(os.path.join(tmpdir, 'rawData'))

    # Copy relevant frames over.
    flatlist = [
        'N20140428S0169.fits'
    ]
    
    flatdarklist = [
        'N20140428S0174.fits'
    ]

    arclist = [
        'N20140428S0085.fits'
    ]

    arcdarklist = [
        'N20140428S0179.fits'
    ]

    ronchilist = [
        'N20140428S0181.fits'
    ]

    scienceFrameList = [
        'N20140428S0086.fits'
    ]

    skyFrameList = [
        'N20140428S0087.fits'
    ]



    for frame in flatlist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in flatdarklist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in arclist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in arcdarklist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in ronchilist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in scienceFrameList:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')
    for frame in skyFrameList:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), 'rawData')

    args = ["config.cfg"]
    nifty.nifsPipeline.start(args)

    with fits.open(os.path.join(tmpdir, 'Titan', '20140428', 'K', 'obs12', 'products_uncorrected', 'ctfbrsnN20140428S0086.fits')) as hdul:
        assert hdul['PRIMARY'].header['DATALAB'] == 'GN-2014A-Q-85-12-002'
        assert hdul[PROVENANCE_EXT_NAME].data.tolist() == [
                ('N20140428S0086.fits', CalibrationTagger.extensionDescriptions['MEMBERSCIENCE'], 'member'),

                ('N20140428S0174.fits', CalibrationTagger.extensionDescriptions['INPUTFLAT'], 'input'),
                ('rnN20140428S0181_ronchi.fits', CalibrationTagger.extensionDescriptions['INPUTRONCHI'], 'input'),
                ('wrn_N20140428S0085_arc.fits', CalibrationTagger.extensionDescriptions['INPUTARC'], 'input'),
                ('N20140428S0087.fits', CalibrationTagger.extensionDescriptions['INPUTSKY'], 'input')
        ]






