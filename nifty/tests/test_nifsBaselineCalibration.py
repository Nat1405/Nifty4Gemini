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
import nifty.pipeline.steps.nifsBaselineCalibration as nifsBaselineCalibration

try:
    from mock import Mock, patch, PropertyMock
except ImportError:
    pytest.skip("Install mock for the cadc tests.", allow_module_level=True)

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

def download_test_data(path):
    r = requests.get("https://github.com/nat1405/NiftyTestData/archive/master.tar.gz")
    with open("NiftyTestData.tar.gz", 'wb') as fd:
        for chunk in r.iter_content(chunk_size=128):
            fd.write(chunk)
    tar = tarfile.open("NiftyTestData.tar.gz")
    tar.extractall()
    tar.close()
    shutil.move(os.path.join(os.getcwd(), "NiftyTestData-master", "GN-2014A-Q-85-all"), data_path("GN-2014A-Q-85-all"))
    os.mkdir(data_path("GN-2014A-Q-85_one_day"))
    for frame in glob.glob(os.path.join(data_path("GN-2014A-Q-85-all"), "N20140428*")):
        shutil.copy(frame, data_path("GN-2014A-Q-85_one_day"))

   
def test_nifsBaselineCalibration(tmpdir, monkeypatch):
    """
        - Make shift image

        - Make NIFS flat field image from 5 input flats, 5 flat darks,
        and check metadata looks okay.

        - Make processed Arc

        - Make processed Ronchi
    """

    tmpdir = str(tmpdir)
    os.mkdir(os.path.join(tmpdir, 'K'))
    os.mkdir(os.path.join(tmpdir, 'Calibrations_K'))
    
    os.chdir(os.path.join(tmpdir, 'Calibrations_K'))
    
    def mocklog(foo):
        pass
    monkeypatch.setattr(logging, "info", mocklog)
    monkeypatch.setattr(logging, "warning", mocklog)
    monkeypatch.setattr(logging, "error", mocklog)

    # Copy config file over and set it up for testing.
    RECIPES_PATH = pkg_resources.resource_filename('nifty', 'recipes/')
    shutil.copy(os.path.join(RECIPES_PATH, 'defaultConfig.cfg'), os.path.join(os.getcwd(), 'config.cfg'))

    if not os.path.exists(data_path('GN-2014A-Q-85-all')):
        download_test_data(data_path('GN-2014A-Q-85-all'))

    with open('config.cfg') as config_file:
        config = ConfigObj(config_file, unrepr=True)
    
    config['ManualMode'] = False
    config['over'] = False
    config['scienceDirectoryList'] = []
    config['telluricDirectoryList'] = []
    config['calibrationDirectoryList'] = [os.getcwd()]

    with open('config.cfg', 'w') as config_file:
        config.write(config_file)

    # Copy relevant calibrations over, and write lists.
    flatlist = [
        'N20140428S0169.fits',
        'N20140428S0170.fits',
        'N20140428S0171.fits',
        'N20140428S0172.fits',
        'N20140428S0173.fits'
    ]
    
    flatdarklist = [
        'N20140428S0174.fits',
        'N20140428S0175.fits',
        'N20140428S0176.fits',
        'N20140428S0177.fits',
        'N20140428S0178.fits'
    ]

    arclist = [
        'N20140428S0085.fits'
    ]

    arcdarklist = [
        'N20140428S0179.fits',
        'N20140428S0180.fits'
    ]

    ronchilist = [
        'N20140428S0181.fits',
        'N20140428S0182.fits'
    ]

    for frame in flatlist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), os.getcwd())
    for frame in flatdarklist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), os.getcwd())
    for frame in arclist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), os.getcwd())
    for frame in arcdarklist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), os.getcwd())
    for frame in ronchilist:
        shutil.copy(data_path(os.path.join('GN-2014A-Q-85-all', frame)), os.getcwd())

    with open('flatlist', 'w') as f:
        f.write('\n'.join([x.split('.')[0] for x in flatlist]))
    with open('flatdarklist', 'w') as f:
        f.write('\n'.join([x.split('.')[0] for x in flatdarklist]))
    with open('arclist', 'w') as f:
        f.write('\n'.join([x.split('.')[0] for x in arclist]))
    with open('arcdarklist', 'w') as f:
        f.write('\n'.join([x.split('.')[0] for x in arcdarklist]))
    with open('ronchilist', 'w') as f:
        f.write('\n'.join([x.split('.')[0] for x in ronchilist]))

    nifsBaselineCalibration.start()

    with fits.open(os.path.join(tmpdir, 'Calibrations_K', 'rgnN20140428S0169_flat.fits')) as hdul:
        assert hdul['PRIMARY'].header['DATALAB'] == 'GN-2014A-Q-85-16-001-RGN-FLAT'

    with fits.open(os.path.join(tmpdir, 'Calibrations_K', 'rgnN20140428S0181_ronchi.fits')) as hdul:
        assert hdul['PRIMARY'].header['DATALAB'] == 'GN-2014A-Q-85-16-013-RGN-RONCHI'

    with fits.open(os.path.join(tmpdir, 'Calibrations_K', 'wrgnN20140428S0085_arc.fits')) as hdul:
        assert hdul['PRIMARY'].header['DATALAB'] == 'GN-2014A-Q-85-12-001-WRGN-ARC'









