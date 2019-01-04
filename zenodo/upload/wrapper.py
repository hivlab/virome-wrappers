__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import os
import requests
import gzip
import re
from subprocess import Popen, PIPE

if not os.environ["ZENODO_PAT"]:
      raise ValueError("Missing ZENODO_PAT environment variable with zenodo API access token!")

# Setup access token
access_token = {"access_token": os.environ["ZENODO_PAT"]}

# Upload depository id
deposition_id = snakemake.params[0]

# Input files
files = snakemake.input

# Create tar.gz file for upload
zipfile = list(set([re.sub("_\d+", "", file) for file in files]))[0] + ".tar.gz"
cmd = ["tar", "-cvzf", zipfile] + files
p = Popen(cmd, stdout = PIPE, stderr = PIPE)
stout, stderr = p.communicate()

# Compose files query and upload url
base_url = "https://zenodo.org/api"
url = os.path.join(base_url, "deposit/depositions/{}/files".format(deposition_id))

# Get info for remote files
r = requests.get(url, params = access_token)
filename = [deposit["filename"] for deposit in r.json()]

# Upload, if file is not present
if os.path.basename(zipfile) not in filename:
    with open(zipfile, "rb") as handle:
        r = requests.post(url, params = access_token,
                          data = {"filename": str(zipfile)},
                          files = {"file": handle})

        if r.status_code != 201:
                raise requests.HTTPError(f"Error in data upload, status code: {r.status_code}   {r.json()['message']}")

else:
    print("Doing nothing. File {} is already uploaded!\nPlease delete local and remote copy of the file\nif you wish to upload new version.".format(os.path.basename(zipfile)))