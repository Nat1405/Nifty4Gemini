import os, glob, shutil
from nifty.pipeline.nifsUtils import downloadQueryCadc


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    return os.path.join(data_dir, filename)

def download_test_data():
    if os.path.exists(data_path("GN-2014A-Q-85-all")):
        shutil.rmtree(data_path("GN-2014A-Q-85-all"))
    os.mkdir(data_path("GN-2014A-Q-85-all"))

    query = "SELECT observationID, type, publisherID, productID \
             FROM caom2.Plane AS Plane JOIN caom2.Observation AS Observation ON Plane.obsID = Observation.obsID \
             WHERE  ( Observation.instrument_name = 'NIFS' \
             AND Observation.collection = 'GEMINI' \
             AND Observation.proposal_id = 'GN-2014A-Q-85' )"

    downloadQueryCadc(query, directory=data_path("GN-2014A-Q-85-all"))

    one_day_files = glob.glob(os.path.join(data_path("GN-2014A-Q-85-all"), "*20140428*"))

    if os.path.exists(data_path("GN-2014A-Q-85_one_day")):
        shutil.rmtree(data_path("GN-2014A-Q-85_one_day"))
    os.mkdir(data_path("GN-2014A-Q-85_one_day"))

    for file in one_day_files:
        shutil.copy(file, data_path("GN-2014A-Q-85_one_day"))
