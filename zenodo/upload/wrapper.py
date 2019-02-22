__author__ = "Taavi Päll"
__copyright__ = "Copyright 2018, Taavi Päll"
__email__ = "tapa741@gmail.com"
__license__ = "MIT"

import os
import requests
import gzip
import re
import hashlib
from subprocess import Popen, PIPE

# Function to calculate md5sum.
def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

if not os.environ["ZENODO_PAT"]:
      raise ValueError("Missing ZENODO_PAT environment variable with zenodo API access token!")

# Setup access token.
access_token = {"access_token": os.environ["ZENODO_PAT"]}

# Upload depository id.
deposition_id = snakemake.params[0]

# Input files.
files = snakemake.input

# Compose files query and upload url.
base_url = "https://zenodo.org/api"
url = os.path.join(base_url, "deposit/depositions/{}/files".format(deposition_id))

# Create a session with a hook to raise error on bad request.
session = requests.Session()
session.hooks = {
   'response': lambda r, *args, **kwargs: r.raise_for_status()
}

# Get info for remote files.
r = session.get(url, params = access_token)
remotefiles = [{"filename": i["filename"], "checksum": i["checksum"], "id": i["id"]} for i in r.json()]

# Create tar.gz file for upload
zipfile = list(set([re.sub("_\d+", "", file) for file in files]))[0] + ".tar.gz"
cmd = ["tar", "-cvzf", zipfile] + files
p = Popen(cmd, stdout = PIPE, stderr = PIPE)
stout, stderr = p.communicate()
if p.returncode != 0:
    errmsg = "%s. Code: %s" % (stderr.strip(), p.returncode)
    raise Exception(errmsg)

# Calculate checksum
checksum = md5(zipfile)

# Check if zipfile or checksum present in remotefiles
remotefile = [i for i in remotefiles if os.path.basename(zipfile) in i["filename"] or checksum in i["checksum"]]

# Upload, if file is not present
if len(remotefile) == 0:
    with open(zipfile, "rb") as handle:
        r = session.post(url, params = access_token, 
                                data = {"filename": str(zipfile)}, 
                               files = {"file": handle})
elif remotefile[0]["checksum"] == checksum and remotefile[0]["filename"] != os.path.basename(zipfile):
    # Rename if checksum matches but not filename
    r = session.put(os.path.join(url, remotefile[0]["id"]), 
                          params = access_token, 
                            data = {"filename": str(zipfile)})
elif remotefile[0]["checksum"] != checksum and remotefile[0]["filename"] == os.path.basename(zipfile):
    # File checksum does not match, delete remote file and upload fresh one
    r = session.delete(os.path.join(url, remotefile[0]["id"]), 
                          params = access_token)
    with open(zipfile, "rb") as handle:
        r = session.post(url, params = access_token, 
                                data = {"filename": str(zipfile)}, 
                               files = {"file": handle})
else:
    print("Doing nothing."
      " File {} is already uploaded!".format(zipfile))
    pass
