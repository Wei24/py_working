import numpy as np
import matplotlib.pyplot as plt
import os, sys
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
from suncasa.utils import helioimage2fits as hf
import os
from astropy.time import TimeDelta
from astropy.time import Time
#task handles
doslfcal=0
doapplyselfcal=0
dofinalclean=1
combine_spws=0
#workdir='/srg/ywei/data/eovsa/sep6/'
workdir='/srg/ywei/data/eovsa/sep6/'
#refms = workdir+'IDB20170906T190319-195320.ms.corrected/'
refms = '/srg/ywei/data/eovsa/sep6/IDB20170906T190319-195320.ms.corrected/'
init_time='2017-09-06T19:19:00.000'
time_interval=10.0
ttotle=90
spws=[str(s+1) for s in range(30)]
antennas=''
pol='XX'
slfcalms_='aa'
spwrans=['1~5','6~12','13~20','21~30']
#for tint in range(ttotle):
for tint in range(24,40):
    #define timerange:
    init_t=Time(init_time,format='isot')
    timed1=TimeDelta((tint*time_interval)*1.0, format='sec')
    timed2=TimeDelta(time_interval, format='sec')
    start_time=init_t+timed1
    end_time=start_time+timed2
    midtime_mjd = (start_time.mjd + end_time.mjd) / 2.
    eph = hf.read_horizons(t0=Time(midtime_mjd, format='mjd'))
    trange = '{0}~{1}'.format(start_time.iso.replace('-', '/').replace(' ', '/'),end_time.iso.replace('-', '/').replace(' ', '/'))
    slfcalms =workdir+'slfcal/No'+str(tint)+'IDB20170906T'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]+'.ms.corrected.xx.slfcal'
    slfcaledms =workdir+'slfcal/No'+str(tint)+'IDB20170906T'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]+'.ms.corrected.xx.slfcaled'
    calprefix=workdir+'slfcal/caltbs/No'+str(tint)+'slf_19'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]
    #imgprefix=workdir+'slfcal/images/No'+str(tint)+'slf_19'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]
    imgprefix=workdir+'slfcal/images/sbs/No'+str(tint)+'slf_19'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]
    if not os.path.exists(slfcalms):
        split(vis=refms,outputvis=slfcalms,datacolumn='data',timerange=trange,correlation='XX')
        listobs(slfcalms,listfile=slfcalms+'.listobs')
    clearcal(slfcalms)
    delmod(slfcalms)
    if doapplyselfcal:
        caltable_list=['/srg/ywei/data/eovsa/caltbs/slf_192920.G0/','/srg/ywei/data/eovsa/caltbs/slf_192920.G1/','/srg/ywei/data/eovsa/caltbs/slf_192920.G2/']
        for nc in range(len(caltable_list)):
            applycal(vis=slfcalms,gaintable=caltable_list[nc],spw=','.join(spws),selectdata=True,antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
            if nc<(len(caltable_list)-1):
                print 'nc',type(slfcalms_)
                slfcalms_=slfcalms+str(nc)
                split(slfcalms,slfcalms_,datacolumn='corrected')
                slfcalms=slfcalms_
            else:
                split(slfcalms,slfcaledms,datacolumn='corrected')
                
                
            
    if doslfcal:
       tb.open(slfcalms+'/SPECTRAL_WINDOW')
       reffreqs=tb.getcol('REF_FREQUENCY')
       bdwds=tb.getcol('TOTAL_BANDWIDTH')
       cfreqs=reffreqs+bdwds/2.
       tb.close()
       sbeam=40.
       strtmp=[t.replace(':','') for t in trange.split('~')]
       timestr='t'+strtmp[0]+'-'+strtmp[1]
       refantenna='0'
       nround=3
       niters=[100,300,300]
       robusts=[1.0,1.0,1.0]
       doapplycal=[1,1,1]
       calmodes=['p','p','a']
       uvranges=['','','']
       slftbs=[]
       spwrans=['1~5','6~12','13~20','21~30']
       rgns=['slfcal_'+trange.split(':')[3]+trange.split(':')[4].split('.')[0]+'_spw'+s.replace('~','-')+'.rgn' for s in spwrans]
       for n in range(nround):
       #for n in range(1,2):
           slfcal_tb_g= calprefix+'.G'+str(n)
           for sp in spws:
               cfreq=cfreqs[int(sp)]
               bm=max(sbeam*cfreqs[1]/cfreq,10.)
               slfcal_img = imgprefix+'.spw'+sp.zfill(2)+'.slfcal'+str(n)
               #apply masks
               if sp < 19:
                   cur_mask='/srg/ywei/data/eovsa/mask/sep_6/larger_' + str(np.floor_divide(sp,3)+1) + '.rgn'
               else:
                   cur_mask='/srg/ywei/data/eovsa/mask/sep_6/larger_6.rgn'
               if n < 2:
                   spbg=max(int(sp)-2,0)
                   sped=min(int(sp)+2,30)
                   spwran=str(spbg)+'~'+str(sped)
                   #if int(sp) < 5:
                   #    spwran='1~5'
                   #if int(sp) >= 5 and int(sp) < 10:
                   #    spwran='4~10'
                   #if int(sp) >= 10 and int(sp) < 15:
                   #    spwran='8~15'
                   #if int(sp) >= 15 and int(sp) < 20:
                   #    spwran='10~20'
                   #if int(sp) >= 20:
                   #    spwran='15~30'
                   #if int(sp) < 5:
                   #    rgn = rgns[0]
                   #if int(sp) > 5 and int(sp) <= 12:
                   #    rgn = rgns[1]
                   #if int(sp) > 12 and int(sp) <= 20:
                   #    rgn = rgns[2]
                   #if int(sp) > 20 and int(sp) <= 30:
                   #    rgn = rgns[3]
                   rgn=rgns[0]
                   interactive=False
                   print 'using spw {0:s} as model'.format(spwran)
                   print 'using rgn {0:s}'.format(rgn)
               else:
                   spwran=sp
                   interactive=False
               try:
                   clean(vis=slfcalms,
                           antenna=antennas,
                           imagename=slfcal_img,
                           uvrange=uvranges[n],
                           #spw=sp,
                           spw=spwran,
                           mode='mfs',
                           timerange=trange,
                           imagermode='csclean',
                           psfmode='clark',
                           imsize=[256,256],
                           cell=['2arcsec'],
                           niter=niters[n],
                           gain=0.05,
                           stokes=pol,
                           weighting='briggs',
                           robust=robusts[n],
                           #phasecenter='J2000 11h14m09 04d52m53',
                           phasecenter='J2000 11h00m48 06d14m60',
                           #mask='box [ [ 75pix , 90pix] , [205pix, 165pix ] ]',
                           mask=cur_mask,
                           #mask='',
                           restoringbeam=[str(bm)+'arcsec'],
                           pbcor=False,
                           interactive=interactive,
                           usescratch=True)
    
               except:
                   print 'error in cleaning spw: '+sp
                   print 'using nearby spws for initial model'
                   sp_e=int(sp)+2
                   sp_i=int(sp)-2
                   if sp_i < 0:
                       sp_i = 0
                   if sp_e > 30:
                       sp_e = 30
                   sp_=str(sp_i)+'~'+str(sp_e)
                   try:
                       tget(clean)
                       spw=sp_
                       clean()
                   except:
                       print 'still not successful. abort...'
                       break
    
               print 'processing spw: '+sp
               #gain solution, phase only
               #if os.path.exists(slfcal_tb_g):
                   # remove flags in the slfcal table
               #    tb.open(slfcal_tb_g,nomodify=False)
               #    flags=tb.getcol('FLAG')
               #    flags=np.zeros(flags.shape)
               #    tb.putcol('FLAG',flags)
               #    tb.close()
               gaincal(vis=slfcalms, refant=refantenna,antenna=antennas,caltable=slfcal_tb_g,spw=sp, uvrange='',\
                       #gaintable=slftbs,selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode='p',\
                       gaintable=[],selectdata=True,timerange=trange,solint='inf',gaintype='G',calmode=calmodes[n],\
                       combine='',minblperant=2,minsnr=2,append=True)
               if not os.path.exists(slfcal_tb_g):
                   #listcal(vis=slfcalms, caltable=slfcal_table)
                   #print 'solutions found in spw: '+sp
                   print 'No solution found in spw: '+sp
    
           if os.path.exists(slfcal_tb_g):
               #slftbs.append(slfcal_tb_g)
               slftbs=[slfcal_tb_g]
               #plotcal(caltable=slfcal_tb_g,antenna=antennas,xaxis='antenna',yaxis='phase',\
                       #subplot=421,plotrange=[0,12,-180,180],iteration='spw')
    
           if doapplycal[n]:
               clearcal(slfcalms)
               delmod(slfcalms)
               applycal(vis=slfcalms,gaintable=slftbs,spw=','.join(spws),selectdata=True,\
                        antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
           if n < nround-1:
                slfcalms_=slfcalms+str(n)
                if os.path.exists(slfcalms_):
                    os.system('rm -rf '+slfcalms_)
                split(slfcalms,slfcalms_,datacolumn='corrected')
                slfcalms=slfcalms_
               #prompt=raw_input('Continuing to selfcal?')
               #if prompt.lower() == 'n':
               #    if os.path.exists(slfcaledms):
               #        os.system('rm -rf '+slfcaledms)
               #    split(slfcalms,slfcaledms,datacolumn='corrected')
               #    print 'Final calibrated ms is {0:s}'.format(slfcaledms)
               #    break
               #if prompt.lower() == 'y':
               #    slfcalms_=slfcalms+str(n)
               #    if os.path.exists(slfcalms_):
               #        os.system('rm -rf '+slfcalms_)
               #    split(slfcalms,slfcalms_,datacolumn='corrected')
               #    slfcalms=slfcalms_
           else:
               if os.path.exists(slfcaledms):
                   os.system('rm -rf '+slfcaledms)
               split(slfcalms,slfcaledms,datacolumn='corrected')
               print 'Final calibrated ms is {0:s}'.format(slfcaledms)

    if dofinalclean:
        if combine_spws:
            spwrans=['1~5','6~12','13~20','21~30']
        else:
            spwrans=map(str,range(1,31))
        #for (rgn,spwran) in zip(rgns,spwrans):
        for spwran in spwrans:
            imname=imgprefix+'_final_'+spwran.replace('~','-')
            fitsfile=imname+'.fits'
            sbeam=40.0
            tb.open(slfcaledms+'/SPECTRAL_WINDOW')
            reffreqs=tb.getcol('REF_FREQUENCY')
            bdwds=tb.getcol('TOTAL_BANDWIDTH')
            cfreqs=reffreqs+bdwds/2.
            tb.close()
            if combine_spws:
                cfreq=(cfreqs[int(spwran.split('~')[0])]+cfreqs[int(spwran.split('~')[1])])/2
            else:
                cfreq=cfreqs[int(spwran)]
            bm=max(sbeam*cfreqs[1]/cfreq,10.)
            print 'cleaning ',str(spwran),'  with beamsize= ', str(bm)
            clean(vis=slfcaledms,
                    antenna=antennas,
                    imagename=imname,
                    spw=spwran,
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
            #hf.imreg(vis=slfcaledms,imagefile=imname+'.image',fitsfile=fitsfile,timerange=trange,usephacenter=True,toTb=True,verbose=False)
            hf.imreg(vis=slfcaledms,ephem=eph,imagefile=imname+'.image',fitsfile=fitsfile,timerange=trange,verbose=False)
            #imname_ori=imgprefix+'_ori_'+spwran.replace('~','-')
            #fitsfile_ori=imname_ori+'.fits'
            #clean(vis=refms,
            #       antenna=antennas,
            #       imagename=imname_ori,
            #       spw=sp,
            #       #mode='channel',
            #       mode='mfs',
            #       timerange=trange,
            #       imagermode='csclean',
            #       psfmode='clark',
            #       imsize=[256,256],
            #       cell=['2arcsec'],
            #       niter=1000,
            #       gain=0.05,
            #       stokes=pol,
            #       weighting='briggs',
            #       robust=0.0,
            #       restoringbeam=[str(bm)+'arcsec'],
            #       #phasecenter='J2000 11h14m09 04d52m53',
            #       phasecenter='J2000 11h00m48 06d14m60',
            #       #mask='clean_mask.rgn',
            #       #mask=cur_mask,
            #       mask='',
            #       pbcor=True,
            #       interactive=False,
            #       usescratch=False)
            #hf.imreg(vis=slfcaledms,imagefile=imname_ori+'.image',fitsfile=fitsfile_ori,timerange=trange,usephacenter=False,toTb=True,verbose=False)
