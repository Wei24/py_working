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
antennas='0~12'
pol='XX'
visa=[]
visa.append('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal')
visa.append('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal0')
visa.append('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal01')
visa.append('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcaled')
tb.open('/srg/ywei/data/eovsa/sep6_dofinal/slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal')
starttim = Time(tb.getcell('TIME', 0) / 24. / 3600., format='mjd')
endtim = Time(tb.getcell('TIME', tb.nrows() - 1) / 24. / 3600., format='mjd')
tb.close()
trange = '{0}~{1}'.format(starttim.iso.replace('-', '/').replace(' ', '/'),endtim.iso.replace('-', '/').replace(' ', '/'))
workdir='/srg/ywei/data/eovsa/sep6_dofinal'
for index in range(4):
    imgprefix=workdir+'/slfcal/images/'+str(index)+'testing_slf_192920'
    img_final=imgprefix+'_final'
    spws=[str(s+1) for s in range(30)]
    slfcaledms=visa[index]
    tb.open(slfcaledms+'/SPECTRAL_WINDOW')
    reffreqs=tb.getcol('REF_FREQUENCY')
    bdwds=tb.getcol('TOTAL_BANDWIDTH')
    cfreqs=reffreqs+bdwds/2.
    tb.close()
    sbeam=35.
    from matplotlib import gridspec as gridspec
    from sunpy import map as smap
    from matplotlib import pyplot as plt
    fig = plt.figure(figsize=(12,10))
    gs = gridspec.GridSpec(6, 5)
    for s,sp in enumerate(spws):
        cfreq=cfreqs[int(sp)]
        bm=max(sbeam*cfreqs[1]/cfreq,10.)
        imname=img_final+'_s'+sp.zfill(2)
        fitsfile=imname+'.fits'
        if sp < 19:
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6/mask_' + str(np.floor_divide(sp,3)+1) + '.rgn'
        else:
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6/mask_6.rgn'
        if not os.path.exists(fitsfile):
            print 'cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm)
            if 1==1:
            #try:
                clean(vis=slfcaledms,
                        antenna=antennas,
                        imagename=imname,
                        spw=sp,
                        #mode='channel',
                        mode='mfs',
                        timerange=trange,
                        imagermode='csclean',
                        psfmode='clark',
                        imsize=[256,256],
                        cell=['2arcsec'],
                        niter=1000,
                        gain=0.05,
                        stokes=pol,
                        weighting='briggs',
                        robust=0.0,
                        restoringbeam=[str(bm)+'arcsec'],
                        #phasecenter='J2000 11h14m09 04d52m53',
                        phasecenter='J2000 11h00m48 06d14m60',
                        #mask='clean_mask.rgn',
                        #mask=cur_mask,
                        mask='',
                        pbcor=True,
                        interactive=False,
                        usescratch=False)
            #except:
                #print 'cleaning spw '+sp+' unsuccessful. Proceed to next spw'
                #continue
            junks=['.flux','.model','.psf','.residual']
            for junk in junks:
                if os.path.exists(imname+junk):
                    os.system('rm -rf '+imname+junk)
            if os.path.exists(imname+'.image'):
                hf.imreg(vis=slfcaledms,imagefile=imname+'.image',fitsfile=fitsfile,
                         timerange=trange,usephacenter=False,toTb=True,verbose=False)
        else:
            ax = fig.add_subplot(gs[s])
            eomap=smap.Map(fitsfile)
            eomap.data=eomap.data.reshape((256,256))
            eomap.plot_settings['cmap'] = plt.get_cmap('jet')
            eomap.plot()
            eomap.draw_limb()
            eomap.draw_grid()
            ax.set_title('spw '+sp)
            ax.set_xlim([850,1150])
            ax.set_ylim([-300,0])
