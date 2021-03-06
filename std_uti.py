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
from astropy.time import TimeDelta
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
import pickle

def makelist(tdir='',keyword1='',keyword2=''):
    li=[]
    for root, dirs, files in os.walk(tdir):
        for file in files:
            if keyword1 in file and keyword2 in file:
                li.append(os.path.join(root, file))
    return li

def make_fov(fitsfile=None,xsize=None,eovsa=None):
    curmap=smap.Map(fitsfile)
    if eovsa==True:
        row = curmap.data.shape[2]
        col = curmap.data.shape[3]
        curmap.data=curmap.data.reshape(curmap.data.shape[2],curmap.data.shape[3])
    else:
        row = curmap.data.shape[0]
        col = curmap.data.shape[1]
    positon = np.nanargmax(curmap.data)
    m, n = divmod(positon, col)
    length = (xsize/2) * u.arcsec
    x0 = curmap.xrange[0] + curmap.scale[1] * (n + 0.5) * u.pix
    y0 = curmap.yrange[0] + curmap.scale[0] * (m + 0.5) * u.pix
    x1 = x0 - length
    x2 = x0 + length
    y1 = y0 - length
    y2 = y0 + length
    fov = [[x1.value, x2.value], [y1.value, y2.value]]
    return fov

def download_aiafits(init_time=None,time_interval=None,wave1=None,wave2=None,instrument=None,path=None):
    client = vso.VSOClient()
    init_t=Time(init_time,format='isot')
    timed=TimeDelta(time_interval, format='sec')
    end_t=init_t+timed
    qr = client.query(vso.attrs.Time(init_t.iso, end_t.iso), vso.attrs.Instrument(instrument),vso.attrs.Wave(wave1 * u.AA, wave2 * u.AA))
    res = client.get(qr, path=path)

def nearest_eovsa_aia(emap=None,alist=None,time=None):
    tdlist=[]
    if not time:
        etime=Time(emap.date).jd
    else:
        etime=time
    for aindex,afits in enumerate(alist):
        tdlist.append(abs(etime-Time(smap.Map(afits).date).jd))
    rindex=tdlist.index(min(tdlist))
    return rindex

# def make_time_dic(workdir=None,radio_dir=None,kw2=None,start_timeindex=None, end_timeindex=None):
#    keyword1_list=['94a','131a','171a','193a','211a','304a','335a']
#    spws=[0,1,2,3]
#    file_list=[]
#    dic_list=[]
#    for kindex,kw1 in enumerate(keyword1_list):
#        cur_list = makelist(tdir=workdir, keyword1=kw1,keyword2=kw2)
#        file_list.append(cur_list)
#    for i in range(start_timeindex,end_timeindex):
#        time_dic={}
#        radio_list=makelist(tdir=radio_dir,keyword1='No'+str(i)+'slf',keyword2=kw2)
#        ref_fits=radio_list[0]
#        ref_time=Time(smap.Map(ref_fits).date).jd
#        time_dic['time']=str(Time(ref_time, format='jd').iso)
#        for kindex,kw1 in enumerate(keyword1_list):
#            cur_index=nearest_eovsa_aia(alist=file_list[kindex], time=ref_time)
#            time_dic[kw1]=file_list[kindex][cur_index]
#        cur_spw_list=[]
#        for spw in spws:
#            cur_spw_list.append(radio_list[spw])
#        time_dic['radio']=cur_spw_list
#        dic_list.append(time_dic)
#    pickle.dump(dic_list, open('/srg/ywei/data/aia/time_dic.p', 'wb'))
#    return dic_list

def draw_goes(goestime=None, ax=None):
    if goestime:
        btgoes = goestime[0]
        etgoes = goestime[1]
    else:
        datstrg = datstr.replace('-', '/')
        btgoes = datstrg + ' ' + qa.time(qa.quantity(tim[0] - 1800, 's'), form='clean', prec=9)[0]
        etgoes = datstrg + ' ' + qa.time(qa.quantity(tim[tidx[-1] - 1] + 1800, 's'), form='clean', prec=9)[0]
    #if verbose:
    #    print 'Acquire GOES soft X-ray data in from ' + btgoes + ' to ' + etgoes

    #ax3 = plt.subplot(gs1[2])

    goesscript = os.path.join(workdir, 'goes.py')
    goesdatafile = os.path.join(workdir, 'goes.dat')
    os.system('rm -rf {}'.format(goesscript))
    fi = open(goesscript, 'wb')
    fi.write('import os \n')
    fi.write('from sunpy.time import TimeRange \n')
    fi.write('from sunpy import lightcurve as lc \n')
    fi.write('import pickle \n')
    fi.write('goesplottim = TimeRange("{0}", "{1}") \n'.format(btgoes, etgoes))
    fi.write('goes = lc.GOESLightCurve.create(goesplottim) \n')
    fi.write('fi2 = open("{}", "wb") \n'.format(goesdatafile))
    fi.write('pickle.dump(goes, fi2) \n')
    fi.write('fi2.close()')
    fi.close()

    try:
        os.system('python {}'.format(goesscript))
    except NameError:
        print "Bad input names"
    except ValueError:
        print "Bad input values"
    except:
        print "Unexpected error:", sys.exc_info()[0]
        print "Error in generating GOES light curves. Proceed without GOES..."

    if os.path.exists(goesdatafile):
        fi1 = file(goesdatafile, 'rb')
        goest = pickle.load(fi1)
        fi1.close()

        dates = mpl.dates.date2num(parse_time(goest.data.index))
        goesdif = np.diff(goest.data['xrsb'])
        gmax = np.nanmax(goesdif)
        gmin = np.nanmin(goesdif)
        ran = gmax - gmin
        db = 2.8 / ran
        goesdifp = goesdif * db + gmin + (-6)
        ax.plot_date(dates, np.log10(goest.data['xrsb']), '-', label='1.0--8.0 $\AA$', color='red', lw=2)
        ax.plot_date(dates[0:-1], goesdifp, '-', label='derivate', color='blue', lw=0.4)

        ax.set_ylim([-7, -3])
        ax.set_yticks([-7, -6, -5, -4, -3])
        ax.set_yticklabels([r'$10^{-7}$', r'$10^{-6}$', r'$10^{-5}$', r'$10^{-4}$', r'$10^{-3}$'])
        ax.set_title('Goes Soft X-ray', fontsize=12)
        ax.set_ylabel('Watts m$^{-2}$')
        ax.set_xlabel(datetime.datetime.isoformat(goest.data.index[0])[0:10])
        ax.axvspan(dates[899], dates[dates.size - 899], alpha=0.2)

        ax.yaxis.grid(True, 'major')
        ax.xaxis.grid(False, 'major')
        ax.legend(prop={'size': 6})

        formatter = mpl.dates.DateFormatter('%H:%M')
        ax.xaxis.set_major_formatter(formatter)
        ax.fmt_xdata = mpl.dates.DateFormatter('%H:%M')

def make_time_dic(workdir=None,radio_dir=None,kw2=None,start_timeindex=None, end_timeindex=None):
    keyword1_list=['94a','131a','171a','193a','211a','304a','335a','gst']
    spws=[0,1,2,3]
    file_list=[]
    dic_list=[]
    for kindex,kw1 in enumerate(keyword1_list):
        cur_list = makelist(tdir=workdir, keyword1=kw1,keyword2=kw2)
        file_list.append(cur_list)

    print 'make image list' 
    time_list=[]
    for kindex,flist in enumerate(file_list):
        cur_number = len(flist)
        jd_list = np.zeros((cur_number))
        for cn in range(cur_number):
            print smap.Map(flist[cn]).date
            jd_list[cn] = Time(smap.Map(flist[cn]).date).jd
        time_list.append(jd_list)


    for i in range(start_timeindex,end_timeindex):
        time_dic={}
        radio_list=makelist(tdir=radio_dir,keyword1='No'+str(i)+'slf',keyword2=kw2)
        ref_fits=radio_list[0]
        ref_time=Time(smap.Map(ref_fits).date).jd
        time_dic['time']=str(Time(ref_time, format='jd').iso)
        for tindex,tlist in enumerate(time_list):
            diff_list=[]
            for itime in tlist:
                diff_list.append(abs(itime-ref_time))
            cur_min=diff_list.index(min(diff_list))
            time_dic[keyword1_list[tindex]]=file_list[tindex][cur_min]            
        cur_spw_list=[]
        for spw in spws:
            cur_spw_list.append(radio_list[spw])
        time_dic['radio']=cur_spw_list
        dic_list.append(time_dic)
    pickle.dump(dic_list, open('/srg/ywei/data/aia/new_time_dic.p', 'wb'))
    return dic_list

def add_to_dic(save_file=None,kw1=None,kw2=None,workdir=None):
    load_dic = pickle.load(open('/srg/ywei/data/aia/update_bbso_time_dic.p', 'rb'))
    file_list=makelist(tdir=workdir, keyword1=kw1,keyword2=kw2)
    time_list=[]
    for cfile in file_list:
        time_list.append(Time(smap.Map(cfile).date).jd)
    # noinspection PyInterpreter
    for cur_dic in load_dic:
        cur_time=Time(cur_dic['time']).jd
        diff_list=[]
        for candi_time in time_list:
            diff_list.append(abs(cur_time-candi_time))
        cur_min=diff_list.index(min(diff_list))
        cur_dic[kw1]=file_list[cur_min]
    pickle.dump(load_dic, open('/srg/ywei/data/aia/update_'+kw1+'_time_dic.p', 'wb'))
    return load_dic

def change_date_format(workdir=None,kw1=None,kw2=None):
    file_list=makelist(tdir=workdir, keyword1=kw1,keyword2=kw2)
    for cfile in file_list:
        print cfile
        with fits.open(cfile,mode='update') as chdul:
            hdr = chdul[0].header
            #wrong_tstring=chdul[0].header['date-obs']
            #hdr['CUNIT1']='arcsec'
            #hdr['CUNIT2']='arcsec'
            #hdr['date_obs']='2017-09-06T'+wrong_tstring.split(" ",2)[2]
            #hdr['date-obs']='2017-09-06T'+wrong_tstring.split(" ",2)[2]
            #hdr['date-obs']=hdr['date-obs']+'T'+hdr['time-obs']
            #hdr['CRVAL1'] = float(hdr['CRVAL1'])
            #hdr['CRVAL2'] = float(hdr['CRVAL2'])
            hdr['CTYPE1'] = 'HPLN-TAN'
            hdr['CTYPE2'] = 'HPLT-TAN'
            chdul.flush()

def regenerating_gst_fits(workdir=None,kw1=None,kw2=None):
    file_list=makelist(tdir=workdir, keyword1=kw1,keyword2=kw2)
    for cfile in file_list:
        gst_file_name=cfile.split('/')[-1]
        os.system('cp /srg/ywei/data/aia/sep06/aia_lev1_171a_2017_09_06t19_21_45_35z_image_lev1.fits /srg/ywei/data/aia/sep06/gst_oriform'+cfile.split('/')[-1])
        cur_file='/srg/ywei/data/aia/sep06/gst_oriform/'+cfile.split('/')[-1]
        print cfile
        print cur_file

        with fits.open(cfile,mode='read') as chdul:
            chdr = chdul[0].header
            cdata = chdul[0].data
        with fits.open(cfile,mode='update') as cur_hdul:
            hdr = cur_hdul[1].header
            cur_hdul[0].data = cdata
            #wrong_tstring=chdul[0].header['date-obs']
            #hdr['CUNIT1']='arcsec'
            #hdr['CUNIT2']='arcsec'
            #hdr['date_obs']='2017-09-06T'+wrong_tstring.split(" ",2)[2]
            #hdr['date-obs']='2017-09-06T'+wrong_tstring.split(" ",2)[2]
            #hdr['date-obs']=hdr['date-obs']+'T'+hdr['time-obs']
            hdr['CRVAL1'] = chdr['CRVAL1']
            hdr['CRVAL2'] = chdr['CRVAL1']
            hdr['CDELT1'] = chdr['CDELT1']
            hdr['CDELT2'] = chdr['CDELT2']
            hdr['CRPIX1'] = chdr['CRPIX1']
            hdr['CRPIX2'] = chdr['CRPIX2']
            hdr['NAXIS1'] = chdr['NAXIS1']
            hdr['NAXIS2'] = chdr['NAXIS2']
            hdr['DATE-OBS'] = chdr['DATE-OBS']
            chdul.flush()



def alignment_plotting(ax=None,ax_title=None,fov=None,radio_maps=None,radio_fits=None,image_maps=None,image_fits=None,eovsamap=None):
    #make maps list with a given fov
    rmaps=[]
    if not radio_maps:
        for rfts in radio_fits:
            #eovsa's data has to be reshaped
            rmap=smap.Map(rfts)
            if eovsamap==True:
                rmap.data=rmap.data.reshape(rmap.data.shape[2],rmap.data.shape[3])
            rmaps.append(rmap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec))
    else:
        for rmap in radio_maps:
            rmaps.append(rmap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec))
    imaps=[]
    if not image_maps:
        for ifts in image_fits:
            imap=smap.Map(ifts)
            imaps.append(imap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec))
    else:
        for imap in image_maps:
            imaps.append(imap.submap(fov[0]*u.arcsec,fov[1]*u.arcsec))           
        
    #radio contour part
    color_str=['r','y','g','b']
    clevels1 = np.linspace(0.7, 0.9, 1)
    print 'there are '+ str(len(rmaps))+ ' maps in rmaps'
    for rindex,rmap in enumerate(rmaps):
        XX, YY = np.meshgrid(np.arange(rmap.data.shape[1]), np.arange(rmap.data.shape[0]))
        rmapx, rmapy = rmap.pixel_to_data(XX * u.pix, YY * u.pix)
        locals()['line'+str(rindex+1)]=ax.contour(rmapx.value, rmapy.value, rmap.data, levels=clevels1 * np.nanmax(rmap.data), colors=color_str[rindex])
    #ax.legend((line1,line2,line3,line4),('lable1','2','3','4'))
    #ax.legend()
    #setting title
    #ax.set_title(ax_title)
    ax.set_title('')
    #image part
    cur_out=imaps[0].plot(axes=ax)
    return cur_out

def one_frame(workdir=None,timeindex=None,single_plot=None):
    if single_plot==True:
        fig = plt.figure(figsize=(12, 8),dpi=100)
    #make fov
    print 'step1'
    cur_fov=make_fov(fitsfile='/srg/ywei/data/eovsa/sep6/slfcal/images/No37slf_192520_final_21-30.fits',xsize=100,eovsa=True)
    #aia wavelength list
    eovsa_fits_dir='/srg/ywei/data/eovsa/sep6/slfcal/images/'
    #aiaw_list=['94','131','171','193','211','304','335','1600','1700','4500']
    aiaw_list=['94','131','171','193','211','304']
    #make eovsa's file list for this moment
    rangelist=['1-5','6-12','13-20','21-30']
    eovsa_list=[]
    print 'step2'
    for rl in rangelist:
        tmplist=makelist(tdir=eovsa_fits_dir,keyword1='No'+str(timeindex)+'slf',keyword2=rl)
        eovsa_list.append(tmplist[0])
    ref_map=smap.Map(eovsa_list[0])
    #draw radio contour and image together
    for aiaw_index, aiawave in enumerate(aiaw_list):
        aiafts_list=makelist(tdir=workdir,keyword1=aiawave+'a',keyword2='fits')
        cur_aiafit=[]
        print 'step3'
        nearest_index=nearest_eovsa_aia(emap=ref_map,alist=aiafts_list)
        print 'step4'
        cur_aiafit.append(aiafts_list[nearest_index])
        cur_ax=fig.add_subplot(3,4,1+aiaw_index)
        cur_title='aia'+aiawave
        print 'step5'
        alignment_plotting(ax=cur_ax,ax_title=cur_title,fov=cur_fov,radio_fits=eovsa_list,image_fits=cur_aiafit,eovsamap=True)

def one_frame_dic(single_plot=None,in_dic=None,fig=None):
    key_list=['94a','131a','171a','193a','211a','304a','magnetogram','continuum','bbso']
    #key_list=['94a','131a','171a','193a','211a','304a','goes']
    if single_plot==True:
        fig = plt.figure(figsize=(12, 8),dpi=100)
    cur_fov=make_fov(fitsfile=in_dic['radio'][0],xsize=100,eovsa=True)
    for ikey_index,ikey in enumerate(key_list):
        if ikey=='goes':
            cur_ax = fig.add_subplot(3, 3, 1 + ikey_index)
            cur_title = ikey
            draw_goes(goestime=['2017/09/06 18:55:00','2017/09/06 19:55:00'], ax=cur_ax)
        image_list=[]
        image_list.append(in_dic[ikey])
        cur_ax=fig.add_subplot(3,3,1+ikey_index)
        cur_title=ikey
        alignment_plotting(ax=cur_ax,ax_title=cur_title,fov=cur_fov,radio_fits=in_dic['radio'],image_fits=image_list,eovsamap=True)
        #if ikey=='bbso':
        #    gst_fov=make_fov(fitsfile=in_dic['bbso'],xsize=120,eovsa=False)
        #    #gst_fov=[[565,605],[-245,-205]]
        #    alignment_plotting(ax=cur_ax,ax_title=cur_title,fov=gst_fov,radio_fits=in_dic['radio'],image_fits=image_list,eovsamap=True)
        #else:
        #    alignment_plotting(ax=cur_ax,ax_title=cur_title,fov=cur_fov,radio_fits=in_dic['radio'],image_fits=image_list,eovsamap=True)

def make_movie(dic_list_file=None):
    dic_list=pickle.load(open(dic_list_file,'rb'))
    #make figure
    fig = plt.figure(figsize=(12, 8),dpi=100)
    for i, cur_dic in enumerate(dic_list):
        one_frame_dic(single_plot=False, in_dic=cur_dic, fig=fig)



#just for test, now, comment it from baozi, comment again from hackintosh, now comment from macbook pro,edit on pycharm
#and more, clone from github directly,test: from pycharm to baozi
#test again from pycharm to baozi again
#now test from baozi to local
