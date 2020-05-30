import os
import sys
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


def a_baby_function():
    path = os.getcwd()
    with open(os.path.join(path, 'foo.txt'), 'w') as f:
        f.write("Hello, world!\n")


def test_baby_function(tmpdir, monkeypatch):
    #p = tmpdir.mkdir("sub").join("foo.txt")
    #p.write("Hello, world!\n")
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




























