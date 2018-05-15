import numpy as np
import matplotlib.pyplot as plt
import os, sys
# from config import get_and_create_download_dir
import shutil
from astropy.io import fits
import urllib2
from split_cli import split_cli as split
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

def makelist(tdir='',keyword1='',keyword2=''):
    li=[]
    for root, dirs, files in os.walk(tdir):
        for file in files:
            if keyword1 in file and keyword2 in file:
                li.append(os.path.join(root, file))
    return li
def seed_grow(im=None,x=None,y=None,threshold=None):
    seed = (x, y)
    mask = np.zeros(im.shape[:2], dtype=np.uint8)
    mask[seed] = 1
    vect = [seed]
    area = [im[seed]]

    while vect:
        n = len(vect)
        mean = np.sum(np.array(area), axis=0) / len(area)
        for i in range(n):
            seed = vect[i]
            s0 = seed[0]
            s1 = seed[1]
            for p in [
                (s0 - 1, s1 - 1),
                (s0 - 1, s1),
                (s0 - 1, s1 + 1),
                (s0, s1 - 1),
                (s0, s1 + 1),
                (s0 + 1, s1 - 1),
                (s0 + 1, s1),
                (s0 + 1, s1 + 1)
            ]:
                if p[0] < 0 or p[0] >= im.shape[0] or p[1] < 0 or p[1] >= im.shape[1]:
                    continue
                if mask[p] == 1:
                    continue
                 #define the threshold
                #if abs(mean - im[p]) <= 5:
                if im[p] >= threshold:
                    mask[p] = 1
                    vect.append(p)
                    area.append(im[p])
        vect = vect[n:]


    #mask = (1 - mask) * 255
    #im = PIL.Image.fromarray(mask)
    #mask.show()
    return mask    

def clean_significance(workdir='',xstart=None,xend=None,ystart=None,yend=None):
    #recognize signal's area automaticlly(distance from brightest point to center point)
    #workdir='/srg/ywei/data/eovsa/sep6_dofinal/slfcal/images'
    fitslist=[]
    fitslist.append(makelist(tdir=workdir,keyword1='0testing',keyword2='fits'))
    fitslist.append(makelist(tdir=workdir,keyword1='1testing',keyword2='fits'))
    fitslist.append(makelist(tdir=workdir,keyword1='2testing',keyword2='fits'))
    fitslist.append(makelist(tdir=workdir,keyword1='3testing',keyword2='fits'))
    round_loop=[]
    for roun in range(4):
        sp_loop=[]
        background_list=[]
        briest_list=[]
        significance=[]
        mask=[]
        dist=[]
        cur_list=fitslist[roun]
        for sp in range(30):
            cur_emap=smap.Map(cur_list[sp])
            sz=cur_emap.data.shape
            backregion=cur_emap.data[0,0,xstart:xend,ystart:yend]
            background_list.append(np.sqrt(np.nanmean(np.power(backregion,2))))
            briest_list.append(np.nanmax(cur_emap.data))
            significance.append(np.nanmax(cur_emap.data)/np.sqrt(np.nanmean(np.power(backregion,2))))
            cur_mask=seed_grow(im=cur_emap.data[0,0,:,:],x=0,y=0,threshold=60000.0)
            #print 'for spw' + str(sp)+', in round ' +str(roun)+' significance is '+ str(np.nanmax(cur_emap.data[0,0,:,:])/np.sqrt(np.nanmean(np.power(backregion,2))))
            #print 'for spw' + str(sp)+', in round ' +str(roun)+' significance is '+ str(np.nanmax(cur_emap.data[0,0,:,2]))
            mask.append(np.nanmean(cur_mask)) 
            cur_im=cur_emap.data[0,0,:,:]
            loc=np.where(cur_im==np.nanmax(cur_im))
            distance=np.sqrt(abs(loc[0][0]-128)**2 + abs(loc[1][0]-128)**2)
            dist.append(distance)
            
            print ' mask now is '+str(distance)
        sp_loop.append(background_list)
        sp_loop.append(briest_list)
        sp_loop.append(significance)
        sp_loop.append(mask)
        sp_loop.append(dist)
        round_loop.append(sp_loop)
    #return round_loop    
    print 'shape now is  ' + str(round_loop[0][2][0])
    rx=[0,1,2,3]
    #plt.xlabel('nround')
    #plt.ylabel('significance')
    #plt.title('spw 10')
    #plt.plot(rx,[round_loop[0][2][9],round_loop[1][2][9],round_loop[2][2][9],round_loop[3][2][9]],'ro')
    rsp=np.empty(30)
    yax = np.empty(30)
    for nn in range(30):
        yax[nn]=abs(round_loop[0][4][nn]-round_loop[3][4][nn])
        rsp[nn]=nn
        print yax[nn]
    plt.plot(rsp,yax,'ro')
    #plt.plot(rx,[round_loop[0][2][9],round_loop[1][2][9],round_loop[2][2][9],round_loop[3][2][9]],'ro')
    plt.show()
        
def make_fov(fitsfile=None,xsize=None):
    curmap=smap.Map(fitsfile)
    row = curmap.data.shape[2]
    col = curmap.data.shape[3]
    positon = np.nanargmax(curmap.data)
    m, n = divmod(positon, col)
    length = (xsize/2) * u.arcsec
    curmap.data=curmap.data.reshape(curmap.data.shape[2],curmap.data.shape[3])
    x0 = curmap.xrange[0] + curmap.scale[1] * (n + 0.5) * u.pix
    y0 = curmap.yrange[0] + curmap.scale[0] * (m + 0.5) * u.pix
    x1 = x0 - length
    x2 = x0 + length
    y1 = y0 - length
    y2 = y0 + length
    fov = [[x1.value, x2.value], [y1.value, y2.value]]
    return fov
#def get_coraiafits(cor_emap=None)
def nearest_eovsa_aia(emap=None,alist=None):
    etime=Time(emap.date).jd
    tdlist=[]
    for aindex,afits in enumerate(alist):        
        tdlist.append(abs(etime-Time(smap.Map(afits).date).jd))
    rindex=tdlist.index(min(tdlist))
    return rindex
#def mk_aiaeovsa_movie()
#   #fig = plt.figure(figsize=(8, 8))
#   #if len(aia_list) > len(eovsa_list)
#   for eindex,efits in enumerate(e_list):
#       fig = plt.figure(figsize=(8, 8))
#       cur_emap=smap.Map(e_list[eindex])
#       aia_index=nearest_eovsa_aia(cur_emap,aia_list)
        
def plt_align(ax=None,afits=None, efits_list=None,xsize=None,fov=None):
    ax.set_title('contour level = 50%')
    amap=smap.Map(afits)
    #ax.set_title(atitle)
    bottom_left = SkyCoord(fov[0][0]*u.arcsec, fov[1][0]*u.arcsec,frame=amap.coordinate_frame)
    up_right = SkyCoord(fov[0][1]*u.arcsec, fov[1][1]*u.arcsec,frame=amap.coordinate_frame)
    #submap = amap.submap(bottom_left, up_right)
    amap=amap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec)
    amap.draw_limb()
    amap.draw_grid()
    clevels1 = np.linspace(0.2, 0.9, 5)
    #clevels1 = np.linspace(0.75, 0.9, 1)
    #clevels1 = np.linspace(0.5, 0.9, 1)
    #clevels1 = np.linspace(0.2, 0.9, 1)
    color_str=['r','y','g','b']
    label_str=['spw1~5','spw6~12','spw13~20','spw21~30']
    lines=[]
    for eindex,efits in enumerate(efits_list):
        cur_emap=smap.Map(efits)
        cur_emap.data=cur_emap.data.reshape(cur_emap.data.shape[2],cur_emap.data.shape[3])
        XX, YY = np.meshgrid(np.arange(cur_emap.data.shape[1]), np.arange(cur_emap.data.shape[0]))
        emapx, emapy = cur_emap.pixel_to_data(XX * u.pix, YY * u.pix)
        #ax.contour(emapx.value, emapy.value, cur_emap.data, levels=clevels1 * np.nanmax(cur_emap.data), cmap=cm.jet)
        ax.contour(emapx.value, emapy.value, cur_emap.data, levels=clevels1 * np.nanmax(cur_emap.data), colors=color_str[eindex])
        cs=ax.contour(emapx.value, emapy.value, cur_emap.data, levels=clevels1 * np.nanmax(cur_emap.data), colors=color_str[eindex])
        lines.append(cs.collections[0])
        #ax.text
    #ax.legend(lines, label_str)
    cur_out=amap.plot(axes=ax)
    #cur_out=amap.draw_rectangle((fov[0][0], fov[1][0]) * u.arcsec, xsize * u.arcsec, xsize * u.arcsec,axes=ax)
    return cur_out
         

def plotting_them(workdir=None,imagesize=None,K11=None,K12=None,K21=None,K22=None,fov=None):
    aiamap_list=makelist(tdir=workdir,keyword1=K11,keyword2=K12)
    eovsa_list=makelist(tdir=workdir,keyword1=K21,keyword2=K22)
    print str(eovsa_list)
    if not fov:
        fov=make_fov(fitsfile=eovsa_list[0],xsize=400)
    fst_emap=smap.Map(eovsa_list[0])
    aia_index=nearest_eovsa_aia(emap=fst_emap,alist=aiamap_list) 
    print 'cur_aiamap= '+ str(aia_index)
    fig = plt.figure(figsize=(12, 12), dpi=100)
    ax = plt.subplot(1, 1, 1)
    ax.set_title('contour level = 50%')
    plt_align(ax=ax,afits=aiamap_list[aia_index], efits_list=eovsa_list,xsize=imagesize,fov=fov)

def mk_movie(workdir=None,imagesize=None,K11=None,K12=None,K21=None,K22=None):
    default_fov = make_fov('/srg/ywei/data/eovsa/sep6/slfcal/images/No64slf_192950_final_1-5.fits',xsize=400)
    for i in range(89):
    #for i in range(20,26):
        K21='No'+str(i)+'slf'
        print 'eovsa index=  '+ str(i)
        plotting_them(workdir=workdir,imagesize=imagesize,K11=K11,K12=K12,K21=K21,K22=K22,fov=default_fov)
        imgfile_name=workdir+'movie/image_'+str(i).zfill(2)+'.png'
        plt.savefig(imgfile_name)
    cur_dir=workdir+'movie/'
    #DButil.img2html_movie(cur_dir)

#def testing_cor():
#    temap=smap.Map('/srg/ywei/data/eovsa/sep6/slfcal/images/No50slf_192730_final_6-12.fits') 
#    amaplist=makelist(tdir='/srg/ywei/data/eovsa/sep6/slfcal/images/',keyword1='aia',keyword2='fits')
#    testing_result=nearest_eovsa_aia(emap=temap,alist=amaplist)
#    print testing_result
#    print str(amaplist[testing_result])

def plot_map_onaxes(ax=None,aiafits=None,in_fits=None,fov=None,axtitle=None):
    amap=smap.Map(aiafits)
    amap=amap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec)
    cur_map=smap.Map(in_fits)
    cur_map.data=cur_map.data.reshape(cur_map.data.shape[2],cur_map.data.shape[3])
    clevels1 = np.linspace(0.2, 0.9, 5)
    XX, YY = np.meshgrid(np.arange(cur_map.data.shape[1]), np.arange(cur_map.data.shape[0]))
    mapx, mapy = cur_map.pixel_to_data(XX * u.pix, YY * u.pix)
    ax.contour(mapx.value, mapy.value, cur_map.data, levels=clevels1 * np.nanmax(cur_map.data), cmap=cm.jet)
    mapx, mapy = cur_map.pixel_to_data(XX * u.pix, YY * u.pix)
    amap.draw_grid()
    #ax.set_title('spw')
    cur_out=amap.plot(axes=ax,title=axtitle) 
    return cur_out

def show_sbs(workdir=None,No=None,K11=None,K12=None,fov=None,xsize=None):
    #fig, ((ax1, ax2,ax3,ax4,ax5,ax6),(ax7, ax8,ax9,ax10,ax11,ax12),(ax13, ax14,ax15,ax16,ax17,ax18),(ax19, ax20,ax21,ax22,ax23,ax24),(ax25, ax26,ax27,ax28,ax29,ax30)) = plt.subplots(nrows=5, ncols=6)
    #fig = plt.figure(figsize=(18, 12),dpi=100)
    #frefits=workdir+'No'+str(No)+'slf_192910_final_1.fits'
    frefits=glob.glob('/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/No'+str(No)+'slf*final_1.fits')
    plt.suptitle(str(Time(smap.Map(frefits).date).iso))
    fremap=smap.Map(frefits[0])
    aiafits_list=makelist(tdir='/srg/ywei/data/eovsa/sep6/slfcal/images/',keyword1=K11,keyword2=K12)
    aindex=nearest_eovsa_aia(emap=fremap,alist=aiafits_list)
    if not fov:
        fov=make_fov(fitsfile=frefits,xsize=xsize)
    for sp in range(1,31):
        efits=glob.glob('/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/No'+str(No)+'slf*final_'+str(sp)+'.fits')[0]
        afits=aiafits_list[aindex]
        cur_ax=fig.add_subplot(5,6,0+sp)
        cur_title='spw'+str(sp)
        plot_map_onaxes(ax=cur_ax,aiafits=afits,in_fits=efits,fov=fov,axtitle=cur_title)

def movie_sbs(xsize=None):
    default_fov=make_fov('/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/No58slf_192850_final_10.fits',xsize=xsize)
    workdir='/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/'
    fig = plt.figure(figsize=(18, 12),dpi=100)
    plt.ioff()
    for No in range(9,89):
        show_sbs(workdir=workdir,No=No,K11='aia',K12='fits',fov=default_fov,xsize=xsize)
        imgfile_name=workdir+'movie/image_'+str(No).zfill(2)+'.png'
        plt.savefig(imgfile_name) 
        
def fits2flux():
    result_array=np.zeros((89,30,256,256))
    time_array=np.zeros(89)
    for tin in range(89):
        ref_fits=glob.glob('/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/No'+str(tin)+'slf*final_1.fits')
        time_array[tin]=Time(smap.Map(ref_fits[0]).date).mjd
        for sp in range(1,31):
            cur_fits=glob.glob('/srg/ywei/data/eovsa/sep6/slfcal/images/sbs/No'+str(tin)+'slf*final_'+str(sp)+'.fits')
            cur_map=smap.Map(cur_fits[0])
            result_array[tin,sp-1,:,:]=cur_map.data[0,0,:,:]
    return (result_array,time_array)
         
def light_curve(x1=None,x2=None,y1=None,y2=None):
    edata=fits2flux()[0]
    etime=fits2flux()[1]
    rarray=np.zeros((edata.shape[0],edata.shape[1]))
    fig = plt.figure(figsize=(14, 9),dpi=100)
    for tin in range(edata.shape[0]):
        for sp in range(edata.shape[1]):
            rarray[tin,sp]=np.nanmean(edata[tin,sp,x1:x2,y1:y2])
    print etime.shape
    print rarray.shape
    for sp in range(edata.shape[1]):
    #for sp in range(13,17):
        cur_ax=fig.add_subplot(5,6,0+sp)
        #cur_ax.set_ylim([0,100000])
        single_lc(ax=cur_ax,xdata=etime,ydata=rarray[:,sp])
    #plt.xlabel('time in jd')
    #plt.ylabel('flux')
    #plt.title('sub_region')
    #plt.legend()
    #plt.show()
    #return rarray

def single_lc(ax=None,xdata=None,ydata=None):
    out=ax.plot(xdata,ydata)
    return out 

def get_freq_wlcm_bm(vis=None):
    tb.open(vis+'/SPECTRAL_WINDOW')
    tb.open(vis+'/SPECTRAL_WINDOW')
    reffreqs=tb.getcol('REF_FREQUENCY')
    bdwds=tb.getcol('TOTAL_BANDWIDTH')
    cfreqs=reffreqs+bdwds/2.
    tb.close()
    sbeam=40.
    bm=np.zeros(30)
    wlcm=np.zeros(30)
    for sp in range(30):
        cfreq=cfreqs[sp]
        wlcm[sp]=30000000000/cfreq
        bm[sp]=max(sbeam*cfreqs[1]/cfreq,10.)
    return (cfreqs,wlcm,bm)

def dign_curve(x1=None,x2=None,y1=None,y2=None):
    edata=fits2flux()[0]
    etime=fits2flux()[1]
    fnb_result=get_freq_wlcm_bm(vis='/srg/ywei/data/eovsa/sep6/slfcal/No60IDB20170906T2910.ms.corrected.xx.slfcaled')
    cfreq=fnb_result[0][0:30]
    wlcm=fnb_result[1]
    bsize=fnb_result[2]
    bt=np.zeros((edata.shape[0],edata.shape[1]))
    for tin in range(1,2):
    #for tin in range(edata.shape[0]):
        for sp in range(edata.shape[1]):
            cur_flux=np.nanmean(edata[tin,sp,x1:x2,y1:y2])
            bt[tin,sp]=1.36*float(cur_flux)*(float(wlcm[sp]/bsize[sp])**2)*1000
        plt.xlabel('frequency in HZ')
        plt.ylabel('brightness tempreature')
        plt.title('')
        plt.loglog(cfreq,bt[tin,:])
