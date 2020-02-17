"""Contains functions for getting sample data from BaseSpace.
"""

import logging
from .utils import api_collection, API_SERVER
from . import bs_project
logger = logging.getLogger(__name__)


def get_files(bs_sample_id, extension='fastq.gz'):
    """Gets a list of ('fastq.gz') files of a sample.

    Args:
        bs_sample_id: BaseSpace ID of a sample
        extension: Filter the files by file extension, default is set to 'fastq.gz'

    Returns: A list of file data, each is a dictionary.

    """
    if not bs_sample_id:
        return []
    files = api_collection("v1pre3/samples/%s/files" % bs_sample_id)
    if files:
        if extension:
            files = [file for file in files if str(file.get("Name", "")).endswith(extension)]
        files = sorted(files, key=lambda f: f.get("Name"))
    return files


def get_fastq_pairs(bs_sample_id):
    """Gets the BaseSpace urls for the pair of FASTQ files of a sample.

    Args:
        bs_sample_id: BaseSpace ID of a sample

    Returns: A list of FASTQ pairs (lists).

    """
    # Returns an empty list if there is no file
    files = get_files(bs_sample_id)
    if not files:
        return []

    # Pairs stores the pairs of FASTQ files as a dictionary of dictionaries.
    # e.g. pairs = {"abc": {"R1": "xxx", "R2": "xxx}, "def": {"R1": "xxx", "R2": "xxx}}
    pairs = {}
    for file in files:
        filename = file.get("Name")
        # Raise an error if both R1 and R2 are in the filename.
        if "_R1_" in filename and "_R2_" in filename:
            raise ValueError("Unable to determine whether the file from sample ID=%s is R1 or R2: %s" % (
                bs_sample_id, filename
            ))
        # Key is used to identify the files from the same pair
        # In the filename, the string before _R1_ or _R2_ should be the same.
        if "_R1_" in filename:
            key = str(filename).rsplit("_R1_", 1)[0]
        elif "_R2_" in filename:
            key = str(filename).rsplit("_R2_", 1)[0]
        else:
            raise ValueError("R1 or R2 is not found in the filename: %s" % filename)
        # Get the URL of the file
        href = file.get("Href")
        uri = API_SERVER + href
        # Save the URI to pairs
        fastq_pair = pairs.get(key, dict())
        if "_R1_" in filename:
            fastq_pair["R1"] = uri
        if "_R2_" in filename:
            fastq_pair["R2"] = uri
        pairs[key] = fastq_pair
    # Convert the "pairs" dictionary to a list of pairs
    pair_list = []
    for fastq_pair in pairs.values():
        r1 = fastq_pair.get("R1")
        r2 = fastq_pair.get("R2")
        if not r1:
            raise ValueError("Missing R1 for file %s in sample ID=%s" % (r2, bs_sample_id))
        if not r2:
            raise ValueError("Missing R2 for file %s in sample ID=%s" % (r1, bs_sample_id))
        pair_list.append([r1, r2])
    return pair_list


def get_sample(project_name, sample_name):
    """Gets the information of a sample.

    Args:
        project_name (str): The name of the project.
        sample_name (str): The name of the sample (Sample ID).

    Returns: A dictionary containing the sample information.

    """
    samples = bs_project.get_samples(project_name)
    for sample in samples:
        if sample.get("Name") == sample_name:
            return sample

    return None


def get_files_by_name(project_name, sample_name):
    """Gets a list of files for a sample

    Args:
        project_name (str): The name of the project.
        sample_name (str): The name of the sample (Sample ID).

    Returns: A list of file information (dictionaries).

    """
    sample = get_sample(project_name, sample_name)
    if sample:
        href = sample.get("Href")
        api_url = href + "/files"
        return api_collection(api_url)
    return None
