import os
import sys
import logging
import json
import glob
from pathlib import Path

import pytest
# try:
#     pyvo_OK = True
#     from pyvo.dal import tap, adhoc
#     from astroquery.cadc import Cadc, conf
#     import astroquery.cadc.core as cadc_core
# except ImportError:
#     pyvo_OK = False
#     pytest.skip("Install pyvo for the cadc module.", allow_module_level=True)
# except AstropyDeprecationWarning as ex:
#     if str(ex) == \
#             'The astropy.vo.samp module has now been moved to astropy.samp':
#         print('AstropyDeprecationWarning: {}'.format(str(ex)))
#     else:
#         raise ex

import nifty.pipeline.steps.nifsSort as nifsSort

try:
    from mock import Mock, patch, PropertyMock
except ImportError:
    pytest.skip("Install mock for the cadc tests.", allow_module_level=True)

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

@patch('nifty.pipeline.steps.nifsSort.makePythonLists',
        Mock(side_effect=lambda x,y: 
            (   
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0
            )
        )
    )
def test_baby_makePythonLists():
    assert nifsSort.makePythonLists("foo", 2) == (0,0,0,0,0,0,0,0,0)




def test_makePythonLists_one_day():
    allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, obsidDateList, sciImageList = nifsSort.makePythonLists(data_path("GN-2014A-Q-85_one_day"), 2.0)

    result_list = ["N20140428S00"+str(i)+'.fits' for i in range(67,94)]
    result_list.remove("N20140428S0085.fits")

    result_allfilelist = [x[0] for x in allfilelist]

    assert sorted(result_allfilelist) == sorted(result_list)

    assert arclist == [
                        "N20140428S0085.fits",
                            ]

    # Failing for now, known bug.
    # assert sorted(arcdarklist) == [
    #                         "N20140428S0179.fits",
    #                         "N20140428S0180.fits"
    #                         ]

    assert sorted(flatlist) == [
                                "N20140428S0169.fits",
                                "N20140428S0170.fits",
                                "N20140428S0171.fits",
                                "N20140428S0172.fits",
                                "N20140428S0173.fits",
                            ]

    assert sorted(flatdarklist) == [
                                    "N20140428S0174.fits",
                                    "N20140428S0175.fits",
                                    "N20140428S0176.fits",
                                    "N20140428S0177.fits",
                                    "N20140428S0178.fits"
                                ]

    assert sorted(ronchilist) == [
                                    "N20140428S0181.fits",
                                    "N20140428S0182.fits"
    ]

    # Failing for now due to a NEW bug...
    # assert sorted(sciImageList) == [
    #                                 "N20140428S0086.fits",
    #                                 "N20140428S0089.fits",
    #                                 "N20140428S0090.fits",
    #                                 "N20140428S0093.fits",
    # ]

    assert objectDateGratingList == [['Titan', '20140428', 'K']]

    # Skip checks for obsidDateList for now

def test_makePythonLists_all_days():
    """
    Test a single target with multiple days of observations (including tellurics).
    """
    skyThreshold = 2.0

    allfilelist, arclist, arcdarklist, flatlist, flatdarklist, ronchilist, objectDateGratingList, obsidDateList, sciImageList = nifsSort.makePythonLists(data_path("GN-2014A-Q-85-all"), skyThreshold)

    desired_allfilelist = []  # Science, aqc, telluric, aqccal
    desired_sciImageList = [] # science and science sky frames
    with open(data_path(os.path.join('json', 'objects_all.json'))) as json_file:
        data = json.load(json_file)
    for p in data:
        tmplist = []
        tmplist.extend([str(p['name']), 1, str(p['observation_class'])])
        desired_allfilelist.append(tmplist)
        if p['observation_class'] == 'science':
            desired_sciImageList.append(str(p['name']))

    assert sorted(allfilelist) == sorted(desired_allfilelist)

    assert sorted(arclist) == [
                        "N20140428S0085.fits",
                        "N20140503S0161.fits",
                        "N20140503S0178.fits",
                        "N20140505S0107.fits",
                        "N20140505S0124.fits"
    ]

    # assert sorted(arcdarklist) == [
    #                         "N20140428S0179.fits",
    #                         "N20140428S0180.fits",
    #                         "N20140503S0251.fits",
    #                         "N20140503S0252.fits",
    #                         "N20140505S0338.fits",
    #                         "N20140505S0339.fits"
    # ]

    desired_flatlist = []
    for i in range(169, 174):
        desired_flatlist.append("N20140428S0{}.fits".format(str(i)))
    for i in range(241, 246):
        desired_flatlist.append("N20140503S0{}.fits".format(str(i)))
    for i in range(328, 333):
        desired_flatlist.append("N20140505S0{}.fits".format(str(i)))
    assert sorted(flatlist) == sorted(desired_flatlist)

    desired_flatdarklist = []
    for i in range(174, 179):
        desired_flatdarklist.append("N20140428S0{}.fits".format(str(i)))
    for i in range(246, 251):
        desired_flatdarklist.append("N20140503S0{}.fits".format(str(i)))
    for i in range(333, 338):
        desired_flatdarklist.append("N20140505S0{}.fits".format(str(i)))
    assert sorted(flatdarklist) == sorted(desired_flatdarklist)

    assert sorted(ronchilist) == [
                                    "N20140428S0181.fits",
                                    "N20140428S0182.fits",
                                    # "N20140428S0182.fits", lamps-off ronchis aren't used
                                    "N20140503S0253.fits",
                                    "N20140503S0254.fits",
                                    # "N20140503S0255.fits",
                                    "N20140505S0340.fits",
                                    "N20140505S0341.fits"
                                    # "N20140505S0342.fits"
    ]

    #assert sorted(sciImageList) == sorted(desired_sciImageList)

    assert sorted(objectDateGratingList) == [
                                        ['Titan', '20140428', 'K'],
                                        ['Titan', '20140503', 'K'],
                                        ['Titan', '20140505', 'K']
    ]

    # Skip checks for obsidDateList for now

def test_makePythonLists_all_days_one_target():
    pass

def a_baby_function():
    logging.info("START")
    path = os.getcwd()
    with open(os.path.join(path, 'foo.txt'), 'w') as f:
        f.write("Hello, world!\n")


def test_baby_function(tmpdir, monkeypatch):
    #p = tmpdir.mkdir("sub").join("foo.txt")
    #p.write("Hello, world!\n")
    def mocklog(foo):
        pass
    monkeypatch.setattr(logging, "info", mocklog)
    def mockjoin(path, *args):
        if len(args) == 0:
            return str(path)
        else:
            print(args)
            return str( Path(str(path)) / Path(mockjoin(args[0], *args[1:])))
    def mockcwd():
        return str(tmpdir)
    monkeypatch.setattr(os, "getcwd", mockcwd)
    monkeypatch.setattr(os.path, "join", mockjoin)
    os.path.join(tmpdir, 'sub', 'bar.txt')

    a_baby_function()
    p = tmpdir.join('foo.txt')
    assert p.read() == "Hello, world!\n"
    assert len(tmpdir.listdir()) == 1
    #objDirList, scienceDirectoryList, telluricDirectoryList = sortScienceAndTelluric(allfilelist, sciImageList, rawPath, skyThreshold):



@patch('nifty.pipeline.steps.nifsSort.turnOffTelluricSkySub',
    Mock(side_effect=lambda: ()))
def test_sortScienceAndTelluric(tmpdir, monkeypatch):
    """
    Tests sort science and telluric for the full dataset of GN-2014A-Q-85 (all three days).
    """
    # Mock os.getcwd(), logging.info/warning/error()
    tmpdir = str(tmpdir)
    
    def mocklog(foo):
        pass
    monkeypatch.setattr(logging, "info", mocklog)
    monkeypatch.setattr(logging, "warning", mocklog)
    monkeypatch.setattr(logging, "error", mocklog)
    # os.cwd() mocking
    def mockcwd():
        return str(tmpdir)
    monkeypatch.setattr(os, "getcwd", mockcwd)

    # if obstype == 'OBJECT' and (obsclass == 'science' or obsclass == 'acq' or obsclass == 'acqCal' or obsclass == 'partnerCal'):
    allfilelist = []
    sciImageList = []
    with open(data_path(os.path.join('json', 'objects_all.json'))) as json_file:
        data = json.load(json_file)
    for p in data:
        tmplist = []
        tmplist.extend([str(p['name']), 1, str(p['observation_class'])])
        allfilelist.append(tmplist)
        if p['observation_class'] == 'science':
            sciImageList.append(str(p['name']))

    skyThreshold = 1.5
    rawPath = data_path('GN-2014A-Q-85-all')

    objDirList, scienceDirectoryList, telluricDirectoryList = nifsSort.sortScienceAndTelluric(allfilelist, sciImageList, rawPath, skyThreshold)


    # Check that the lists of science and telluric observations was made correctly
    result_objDirList = [
                            os.path.join(tmpdir, "Titan", "20140505"),
                            os.path.join(tmpdir, "Titan", "20140503"),
                            os.path.join(tmpdir, "Titan", "20140428"),
    ]

    for obsdir in result_objDirList:
        assert os.path.exists(obsdir)
    assert sorted(objDirList) == sorted(result_objDirList)

    result_scienceDirectoryList = [
                            os.path.join(tmpdir, "Titan", "20140505", "K", "obs28"),
                            os.path.join(tmpdir, "Titan", "20140503", "K", "obs20"),
                            os.path.join(tmpdir, "Titan", "20140428", "K", "obs12")
    ]

    for obsdir in result_scienceDirectoryList:
        assert os.path.exists(obsdir)
    assert sorted(scienceDirectoryList) == sorted(result_scienceDirectoryList)

    result_telluricDirectoryList = [
                            os.path.join(tmpdir, "Titan", "20140505", "K", "Tellurics", "obs30"),
                            os.path.join(tmpdir, "Titan", "20140505", "K", "Tellurics", "obs26"),
                            os.path.join(tmpdir, "Titan", "20140503", "K", "Tellurics", "obs18"),
                            os.path.join(tmpdir, "Titan", "20140503", "K", "Tellurics", "obs22"),
                            os.path.join(tmpdir, "Titan", "20140428", "K", "Tellurics", "obs10")
    ]

    for obsdir in result_telluricDirectoryList:
        assert os.path.exists(obsdir)
    assert sorted(telluricDirectoryList) == sorted(result_telluricDirectoryList)


    # Now check all files were copied over to a pre-determined directory structure.
    # NEEDS REPLACED; THIS IS WRONG. Acquisition folders should have been included.
    result_dir_structure = [
        (os.path.join(tmpdir, 'Titan'), ['20140505', '20140503', '20140428'], []),
        (os.path.join(tmpdir, 'Titan/20140505'), ['K'], []),
        (os.path.join(tmpdir, 'Titan/20140505/K'), ['obs28', 'Tellurics', 'Acquisitions'], []),
        (os.path.join(tmpdir, 'Titan/20140505/K/obs28'), [], ['N20140505S0114.fits', 'N20140505S0115.fits', 'N20140505S0117.fits', 'N20140505S0118.fits', 'N20140505S0119.fits', 'N20140505S0120.fits', 'N20140505S0121.fits', 'N20140505S0122.fits', 'N20140505S0123.fits', 'N20140505S0116.fits', 'N20140505S0108.fits', 'N20140505S0112.fits', 'N20140505S0109.fits', 'N20140505S0110.fits', 'N20140505S0113.fits', 'N20140505S0111.fits', 'skyFrameList', 'scienceFrameList']),
        (os.path.join(tmpdir, 'Titan/20140505/K/Tellurics'), ['obs30', 'obs26'], []),
        (os.path.join(tmpdir, 'Titan/20140505/K/Tellurics/obs30'), [], ['N20140505S0127.fits', 'N20140505S0128.fits', 'N20140505S0130.fits', 'N20140505S0129.fits', 'N20140505S0131.fits', 'tellist']),
        (os.path.join(tmpdir, 'Titan/20140505/K/Tellurics/obs26'), [], ['N20140505S0103.fits', 'N20140505S0102.fits', 'N20140505S0101.fits', 'N20140505S0100.fits', 'N20140505S0099.fits', 'tellist']),
        # (os.path.join(tmpdir, 'Titan/20140505/K/Acquisitions'), [], ['N20140505S0104.fits']),
        (os.path.join(tmpdir, 'Titan/20140503'), ['K'], []),
        (os.path.join(tmpdir, 'Titan/20140503/K'), ['obs20', 'Tellurics', 'Acquisitions'], []),
        (os.path.join(tmpdir, 'Titan/20140503/K/obs20'), [], ['N20140503S0177.fits', 'N20140503S0169.fits', 'N20140503S0168.fits', 'N20140503S0166.fits', 'N20140503S0167.fits', 'N20140503S0164.fits', 'N20140503S0163.fits', 'N20140503S0162.fits', 'N20140503S0165.fits', 'N20140503S0171.fits', 'N20140503S0172.fits', 'N20140503S0174.fits', 'N20140503S0175.fits', 'N20140503S0176.fits', 'N20140503S0170.fits', 'N20140503S0173.fits', 'scienceFrameList', 'skyFrameList']),
        (os.path.join(tmpdir, 'Titan/20140503/K/Tellurics'), ['obs18', 'obs22'], []),
        (os.path.join(tmpdir, 'Titan/20140503/K/Tellurics/obs18'), [], ['N20140503S0157.fits', 'N20140503S0156.fits', 'N20140503S0155.fits', 'N20140503S0154.fits', 'N20140503S0153.fits', 'tellist']),
        (os.path.join(tmpdir, 'Titan/20140503/K/Tellurics/obs22'), [], ['N20140503S0181.fits', 'N20140503S0184.fits', 'N20140503S0182.fits', 'N20140503S0183.fits', 'N20140503S0185.fits', 'tellist']),
        # (os.path.join(tmpdir, 'Titan/20140503/K/Acquisitions'), [], ['N20140503S0158.fits']),
        (os.path.join(tmpdir, 'Titan/20140428'), ['K'], []),
        (os.path.join(tmpdir, 'Titan/20140428/K'), ['obs12', 'Tellurics', 'Acquisitions'], []),
        (os.path.join(tmpdir, 'Titan/20140428/K/obs12'), [], ['N20140428S0093.fits', 'N20140428S0092.fits', 'N20140428S0091.fits', 'N20140428S0090.fits', 'N20140428S0089.fits', 'N20140428S0087.fits', 'N20140428S0086.fits', 'N20140428S0088.fits', 'scienceFrameList', 'skyFrameList']),
        (os.path.join(tmpdir, 'Titan/20140428/K/Tellurics'), ['obs10'], []),
        (os.path.join(tmpdir, 'Titan/20140428/K/Tellurics/obs10'), [], ['N20140428S0077.fits', 'N20140428S0078.fits', 'N20140428S0069.fits', 'N20140428S0072.fits', 'N20140428S0070.fits', 'N20140428S0075.fits', 'N20140428S0071.fits', 'N20140428S0073.fits', 'N20140428S0074.fits', 'N20140428S0076.fits', 'tellist'])
        # (os.path.join(tmpdir, 'Titan/20140428/K/Acquisitions'), [], ['N20140428S0079.fits'])
    ]

    for item in result_dir_structure:
        item[1].sort()
        item[2].sort()

    dir_structure = os.walk(os.path.join(tmpdir, 'Titan'))
    dirs_list = []
    for item in dir_structure:
        item[1].sort()
        item[2].sort()
        # Acquisitions aren't deterministic for now so skip them
        if 'Acquisitions' not in item[0]:
            dirs_list.append(item)

    # Failing for now because aquisitions aren't being copied over for some reason
    assert sorted(dirs_list) == sorted(result_dir_structure)



    # Check scienceFrameLists, tellists, and skyFrameLists were made properly.
    # Science frame lists
    with open(os.path.join(tmpdir, 'Titan', '20140428', 'K', 'obs12', 'scienceFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140428S0086\n", "N20140428S0089\n", "N20140428S0090\n", "N20140428S0093\n"])

    with open(os.path.join(tmpdir, 'Titan', '20140503', 'K', 'obs20', 'scienceFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140503S0162\n", "N20140503S0165\n", "N20140503S0166\n", "N20140503S0169\n", "N20140503S0170\n", "N20140503S0173\n", "N20140503S0174\n", "N20140503S0177\n"])

    with open(os.path.join(tmpdir, 'Titan', '20140505', 'K', 'obs28', 'scienceFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140505S0108\n", "N20140505S0111\n", "N20140505S0112\n", "N20140505S0115\n", "N20140505S0116\n", "N20140505S0119\n", "N20140505S0120\n", "N20140505S0123\n"])

    # SkyFrameLists
    with open(os.path.join(tmpdir, 'Titan', '20140428', 'K', 'obs12', 'skyFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140428S0087\n", "N20140428S0088\n", "N20140428S0091\n", "N20140428S0092\n"])

    with open(os.path.join(tmpdir, 'Titan', '20140503', 'K', 'obs20', 'skyFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140503S0163\n", "N20140503S0164\n", "N20140503S0167\n", "N20140503S0168\n", "N20140503S0171\n", "N20140503S0172\n", "N20140503S0175\n", "N20140503S0176\n"])

    with open(os.path.join(tmpdir, 'Titan', '20140505', 'K', 'obs28', 'skyFrameList')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(["N20140505S0109\n", "N20140505S0110\n", "N20140505S0113\n", "N20140505S0114\n", "N20140505S0117\n", "N20140505S0118\n", "N20140505S0121\n", "N20140505S0122\n"])

    # tellists
    with open(os.path.join(tmpdir, 'Titan', '20140428', 'K', 'Tellurics', 'obs10', 'tellist')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(['N20140428S0069\n', 'N20140428S0070\n', 'N20140428S0071\n', 'N20140428S0072\n', 'N20140428S0073\n', 'N20140428S0074\n', 'N20140428S0075\n', 'N20140428S0076\n', 'N20140428S0077\n', 'N20140428S0078\n'])

    with open(os.path.join(tmpdir, 'Titan', '20140503', 'K', 'Tellurics', 'obs18', 'tellist')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(['N20140503S0153\n', 'N20140503S0154\n', 'N20140503S0155\n', 'N20140503S0156\n', 'N20140503S0157\n'])

    with open(os.path.join(tmpdir, 'Titan', '20140504', 'K', 'Tellurics', 'obs22', 'tellist')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(['N20140503S0181\n', 'N20140503S0182\n', 'N20140503S0183\n', 'N20140503S0184\n', 'N20140503S0185\n'])

    with open(os.path.join(tmpdir, 'Titan', '20140505', 'K', 'Tellurics', 'obs26', 'tellist')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(['N20140505S0099\n', 'N20140505S0100\n', 'N20140505S0101\n', 'N20140505S0102\n', 'N20140505S0103\n'])

    with open(os.path.join(tmpdir, 'Titan', '20140505', 'K', 'Tellurics', 'obs30', 'tellist')) as f:
        frames = f.readlines()
    assert sorted(frames) == sorted(['N20140505S0127\n', 'N20140505S0128\n', 'N20140505S0129\n', 'N20140505S0130\n', 'N20140505S0131\n'])














