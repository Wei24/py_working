import numpy as np
import matplotlib.pyplot as plt
import os, sys
# from config import get_and_create_download_dir
import shutil
from astropy.io import fits
import urllib2
from split_cli import split_cli as split
from ptclean_cli import ptclean_cli as ptclean
from suncasa.utils import helioimage2fits as hf
import sunpy.map as smap
from sunpy.net import vso
from astropy import units as u
from astropy.time import Time
from taskinit import ms, tb, qa, iatool
from clean_cli import clean_cli as clean
from sunpy import lightcurve
from sunpy.time import TimeRange
from matplotlib.dates import DateFormatter
from astropy.io import fits
from astropy.coordinates import SkyCoord
from sunpy import lightcurve as lc
from sunpy.time import TimeRange, parse_time
import pickle
import datetime
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.patches as patches
from matplotlib import gridspec
import glob
from suncasa.utils import DButil
import copy
import pdb

fig = plt.figure(figsize=(12, 16), dpi=100)
gs1 = gridspec.GridSpec(4, 4)
gs1.update(left=0.08, right=0.32, wspace=0.05)
tdir='./fitsfile123/'
tb.open('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal0')
starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
endtim = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
tb.close()
timerange = '{0}~{1}'.format(starttim.iso.replace('-', '/').replace(' ', '/'),endtim.iso.replace('-', '/').replace(' ', '/'))
cur_imagefile=[]
cur_fitsfile=[]
vis='/srg/ywei/data/eovsa/sep6_dofinal/IDB20170906T190319-195320.ms.corrected'
for sp in range(4):
    index=0
    cur_imagefile.append('slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal0.image/')
    cur_fitsfile.append(tdir+'slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal0.fits')
    print cur_imagefile[0]
    hf.imreg(vis=vis,imagefile=cur_imagefile[0],fitsfile=cur_fitsfile[0],timerange=timerange)
    cur_imagefile.append('slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal1.image/')
    cur_fitsfile.append(tdir+'slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal1.fits')
    hf.imreg(vis=vis,imagefile=cur_imagefile[1],fitsfile=cur_fitsfile[1],timerange=timerange)
    cur_imagefile.append('slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal2.image/')
    cur_fitsfile.append(tdir+'slf_192920.spw'+str(sp+index*4+1).zfill(2)+'.slfcal2.fits')
    hf.imreg(vis=vis,imagefile=cur_imagefile[2],fitsfile=cur_fitsfile[2],timerange=timerange)
    cur_fitsfile.append('./without_mask_final/slf_192920_final_s'+str(sp+index*4+1).zfill(2)+'.fits')
    cur_ax=[] 
    emap=[]
    for order in range(4):
        cur_ax.append(plt.subplot(gs1[order+sp*4]))
        emap.append(smap.Map(cur_fitsfile[order]))
        emap[order].data=emap[order].data.reshape((256,256))
        emap[order].plot_settings['cmap'] = plt.get_cmap('jet')
        emap[order].plot(axes=cur_ax[order])
