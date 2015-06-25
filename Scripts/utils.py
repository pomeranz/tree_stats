import os
from os import path
import httplib2
import urllib
import json
from itertools import izip_longest

# from Bio import PDB

def check_dir(dirname):
    if not path.exists(dirname):
        os.mkdir(dirname)

# Ensembl REST API stuff
http = httplib2.Http(".cache")

server = "http://rest.ensembl.org/"
def ens_get(ext, *args, **kwargs):
    if len(args):
        ext += '&'.join([ urllib.quote(a) for a in args]) 
        
    if len(kwargs):
        ext += '&' + urllib.urlencode(kwargs)

    print ext

    resp, content = http.request(server+ext, method="GET", headers={"Content-Type":"application/json"})
    
    if not resp.status == 200:
        exc = IOError(resp)
        exc.errno = resp
        raise exc

    decoded = json.loads(content)
    return decoded


def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)
