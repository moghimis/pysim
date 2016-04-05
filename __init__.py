#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Init file for assimilation
"""
__author__    = "Saeed Moghimi"
__copyright__ = "Copyright 2015, Oregon State University"
__license__   = "GPL"
__version__   = "0.1"
__email__     = "moghimis@gmail.com"

from __init__ import *
from base_info import *
#############################################
# Saeed Moghimi; moghimis@gmail.com   
# Logs:
# 1.0  
# 3.0 
# 4.0 
#########################################################

import os
import datetime
import time
import base_info 
import numpy as np

logfile = base_info.base_dir+'/main_log'+datetime.datetime.now().isoformat()+'.txt'
if (not os.path.exists(base_info.base_dir)):
    os.system('mkdir -p '+base_info.base_dir )

def logf(txt1,txt2,logfile):
    str1=txt1+' > '+txt2+'\n'
    os.system('echo "'+ str1+ '" >> '+ logfile )
    os.system('echo "'+ str1+ '"')

ijxy = np.append(base_info.ix,base_info.jy)
local_inp = base_info.base_dir+'/inp'


print 'init assim imported ...'  



 
