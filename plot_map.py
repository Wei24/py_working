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

def WEI_plot(dofile=True ,vis=None, timerange=None, spw='', aiafits='', imagehead='', workdir='', spwCol=3, phasecenter='J2000 11h00m48 06d14m60'):
    if vis[-1] == '/':
        vis = vis[:-1]
    ms.open(vis)
    spwInfo=ms.getspectralwindowinfo()
    ms.close()
    tb.open(vis)
    starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
    endtim = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
    tb.close()
    tb.open(vis+'/SPECTRAL_WINDOW')
    reffreqs=tb.getcol('REF_FREQUENCY')
    bdwds=tb.getcol('TOTAL_BANDWIDTH')
    cfreqs=reffreqs+bdwds/2.
    tb.close()
    sbeam=35.
    #get timerange from vis file
    if not timerange:
        timerange = '{0}~{1}'.format(starttim.iso.replace('-', '/').replace(' ', '/'),
                                     endtim.iso.replace('-', '/').replace(' ', '/'))
    nspw=len(spwInfo)
    #draw plot aia
    fig = plt.figure(figsize=(12, 7), dpi=100)
    gs1 = gridspec.GridSpec(4, 3)
    gs1.update(left=0.08, right=0.32, wspace=0.05)
    ax1=plt.subplot(gs1[11])
    aiamap=smap.Map(aiafits)
    aiamap.plot(axes=ax1)
    #do clean spwCol by spwCol
    for cur_spwCol in range(0,np.floor_divide(nspw,spwCol)):
        if ((cur_spwCol+1)*spwCol) < nspw:
            cur_spwRange=str(cur_spwCol*spwCol+1) + '~' + str((cur_spwCol+1)*spwCol)
        else:
            cur_spwRange=str(cur_spwCol*spwCol+1) + '~' + '31'
        imagename=imagehead+cur_spwRange+'SPWs'
        cur_eovsaFits = imagename+'.fits'
        if cur_spwCol < 6:
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6mask_' + str(cur_spwCol+1) + '.rgn'
        else:
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6/mask_6.rgn'
        if dofile:
            #clean(vis=vis, spw=cur_spwRange, timerange=timerange, imagename=imagename, imsize=[256,256], niter=100, cell=['2arcsec'] )
            #clean(vis=vis, spw=cur_spwRange, timerange=timerange, imagename=imagename, imsize=[512,512], niter=1000, cell=['1arcsec'],stokes='XX', gain=0.05,weighting='briggs', mode='mfs',imagermode='csclean',psfmode='clark',robust=0.0,restoringbeam = ['10.0arcsec'], mask=cur_mask,pbcor=True)
            clean(vis=vis, spw=cur_spwRange, timerange=timerange, imagename=imagename, imsize=[512,512], niter=1000, cell=['1arcsec'],stokes='XX', gain=0.05,weighting='briggs', mode='mfs',imagermode='csclean',psfmode='clark',robust=0.0,restoringbeam = ['10.0arcsec'], mask='',pbcor=True)
            print 'fits name =' + str(cur_eovsaFits)
            hf.imreg(vis=vis,imagefile=imagename+'.image',fitsfile=imagename+'.fits',timerange=timerange)
        #plot eovsa
        cur_emap=smap.Map(cur_eovsaFits)
        cur_axname=plt.subplot(gs1[cur_spwCol+1])
        (npol,nf,nx,ny)=cur_emap.data.shape
        print 'shape = ' + str(cur_emap.data.shape)
        if npol != 1:
            print 'To be determined'
        else:
            cur_emap.data=cur_emap.data.reshape((512,512))
            cur_emap.plot_settings['cmap'] = plt.get_cmap('jet')
            cur_emap.plot(axes=cur_axname)
        #cur_emap.plot(axes=cur_axname,title='') 
    
    
        
        
