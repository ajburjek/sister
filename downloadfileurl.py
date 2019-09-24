#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Repeating Ground Track Orbits in High-Fidelity Geopotential"""
import requests
"""Python 3.7
   Simpson Aerospace, Copyright 2019
   Christopher R. Simpson
   simpsonchristo@gmail.com """
#------------------------------------------------------------------------------
def download_file(url):
    local_filename = url.split('/')[-1]
    # NOTE the stream=True parameter below
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filename, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192): 
                if chunk: # filter out keep-alive new chunks
                    f.write(chunk)
                    # f.flush()
    return local_filename