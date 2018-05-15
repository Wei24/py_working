from suncasa.utils import helioimage2fits as hf
import os
import numpy as np
#task handles
dofullsun=0
dopartsun=0
dofullsun_slfed=0
doslfcal=1
dosbd=0
dofinalclean=0
workdir='/srg/ywei/data/eovsa/sep6_dofinal/'
refms = workdir+'IDB20170906T190319-195320.ms.corrected/'
slfcalms =workdir+'slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcal'
slfcaledms = workdir+'slfcal/IDB20170906T190319-195320.ms.corrected.xx.slfcaled'
if not os.path.exists(slfcalms):
    split(vis=refms,outputvis=slfcalms,datacolumn='data',timerange='19:29:20~19:29:30',correlation='XX')
    listobs(slfcalms,listfile=slfcalms+'.listobs')
#W24:Re-initializes the calibration for a visibility data set
clearcal(slfcalms)
#W24:Deletes model representations in the MS
delmod(slfcalms)
trange='2017/09/06/19:29:20~2017/09/06/19:29:30'
spws=[str(s+1) for s in range(30)]
antennas='0~12' 
pol='XX'
calprefix=workdir+'slfcal/caltbs/slf_192920'
imgprefix=workdir+'slfcal/images/slf_192920'

if dofullsun:
    #initial mfs clean to find out the image phase center
    img_init=imgprefix+'_init'
    os.system('rm -rf '+img_init+'*')
    clean(vis=slfcalms,
            antenna='0~12',
            imagename=img_init,
            spw='1~15',
            mode='mfs',
            timerange=trange,
            imagermode='csclean',
            psfmode='clark',
            imsize=[512,512],
            cell=['5arcsec'],
            niter=1000,
            gain=0.05,
            stokes=pol,
            restoringbeam=['15arcsec'],
            interactive=False,
            pbcor=True,
            usescratch=True)

    hf.imreg(vis=slfcalms,imagefile=img_init+'.image',fitsfile=img_init+'.fits',
             timerange=trange,usephacenter=False,verbose=True)

    from sunpy import map as smap
    from matplotlib import pyplot as plt
    eomap=smap.Map(img_init+'.fits')
    eomap.data=eomap.data.reshape((512,512))
    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
    eomap.plot()
    eomap.draw_limb()
    eomap.draw_grid()
    plt.show()
    viewer(img_init+'.image')

if dopartsun:
    img_init=imgprefix+'_init_partial'
    os.system('rm -rf '+img_init+'*')
    spwrans=['1~5','6~12','13~20','21~30']
    rgns=['slfcal_164100_spw'+s.replace('~','-')+'.rgn' for s in spwrans]
    imnames=[]
    for (rgn,spwran) in zip(rgns,spwrans):
        imname=img_init+spwran.replace('~','-')
        try:
            clean(vis=slfcalms,
                    antenna='0~12',
                    imagename=img_init+spwran.replace('~','-'),
                    spw=spwran,
                    mode='mfs',
                    timerange=trange,
                    imagermode='csclean',
                    psfmode='clark',
                    imsize=[256,256],
                    cell=['2arcsec'],
                    niter=1000,
                    gain=0.05,
                    stokes=pol,
		    #W24
                    restoringbeam=['20arcsec'],
		    #['11:00:48.861,06.14.60.451']
                    # phasecenter='J2000 11h14m09 04d52m53',
                    phasecenter='J2000 11h00m48 06d14m60',
                    weighting='briggs',
                    robust=1.0,
                    interactive=False,
                    pbcor=True,
                    usescratch=False)
            imnames.append(imname)
        except:
            print('error in cleaning spw: '+spwran)
    #viewer(imnames)

if doslfcal:
    os.system('rm -rf '+calprefix+'*')
    os.system('rm -rf '+imgprefix+'*')
    #trange='14:20:00~14:30:00'
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
    rgns=['slfcal_192920_spw'+s.replace('~','-')+'.rgn' for s in spwrans]
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
		#interactive=True
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
            plotcal(caltable=slfcal_tb_g,antenna=antennas,xaxis='antenna',yaxis='phase',\
                    subplot=421,plotrange=[0,12,-180,180],iteration='spw')

        if doapplycal[n]:
            clearcal(slfcalms)
            delmod(slfcalms)
            applycal(vis=slfcalms,gaintable=slftbs,spw=','.join(spws),selectdata=True,\
                     antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)
        if n < nround-1: 
            prompt=raw_input('Continuing to selfcal?')
            if prompt.lower() == 'n':
                if os.path.exists(slfcaledms):
                    os.system('rm -rf '+slfcaledms)
		#W24    
                split(slfcalms,slfcaledms,datacolumn='corrected')
                print 'Final calibrated ms is {0:s}'.format(slfcaledms)
                break
            if prompt.lower() == 'y':
                slfcalms_=slfcalms+str(n)
                if os.path.exists(slfcalms_):
                    os.system('rm -rf '+slfcalms_)
                split(slfcalms,slfcalms_,datacolumn='corrected')
                slfcalms=slfcalms_
        else:
            if os.path.exists(slfcaledms):
                os.system('rm -rf '+slfcaledms)
            split(slfcalms,slfcaledms,datacolumn='corrected')
            print 'Final calibrated ms is {0:s}'.format(slfcaledms)

#split(slfcalms,outputvis=slfcaledms,datacolumn='corrected')
if dosbd: #this is for obtaining single band delays (across all channels within a band), not used for this event
    slfcal_tb_k = calprefix+'.SBD'
    gaincal(vis=slfcalms, refant=refantenna,antenna=antennas,caltable=slfcal_tb_k, spw=','.join(spws), uvrange='',
            gaintable=slftbs,selectdata=True,timerange=trange,solint='inf',gaintype='K') 
    plotcal(caltable=slfcal_tb_k,antenna=antennas,xaxis='freq',yaxis='delay',\
            subplot=421,plotrange=[-1,-1,-2,2],iteration='antenna')
    slfcal_tb_k_mod=slfcal_tb_k+'.mod'
    os.system('cp -r '+slfcal_tb_k+' '+slfcal_tb_k_mod) 
    tb.open(slfcal_tb_k_mod,nomodify=False)
    subtb=tb.query('SPECTRAL_WINDOW_ID == 1')
    par=subtb.getcol('FPARAM')
    for sp in spws:
        subtb=tb.query('SPECTRAL_WINDOW_ID == '+sp)
        if sp != '1':
            subtb.putcol('FPARAM',par)
    tb.close()
    slftbs.append(slfcal_tb_k_mod)
    clearcal(slfcalms)
    delmod(slfcalms)
    applycal(vis=slfcalms,gaintable=slftbs,spw=','.join(spws),selectdata=True,\
             antenna=antennas,interp='nearest',flagbackup=False,applymode='calonly',calwt=False)

if dofinalclean:
    imgprefix=workdir+'/slfcal/images/slf_192920'
    img_final=imgprefix+'_final'
    spws=[str(s+1) for s in range(30)]
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
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6/larger_' + str(np.floor_divide(sp,3)+1) + '.rgn'
        else:
            cur_mask='/srg/ywei/data/eovsa/mask/sep_6/larger_6.rgn'
        if not os.path.exists(fitsfile):
            print 'cleaning spw {0:s} with beam size {1:.1f}"'.format(sp,bm)
            try:
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
            except:
                print 'cleaning spw '+sp+' unsuccessful. Proceed to next spw'
                continue
            junks=['.flux','.model','.psf','.residual','.mask']
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

    plt.show()

if dofullsun_slfed:
    img_slfed=imgprefix+'_slfed'
    os.system('rm -rf '+img_slfed+'*')
    clean(vis=slfcaledms,
            antenna='0~12',
            imagename=img_slfed,
            spw='1~5',
            mode='mfs',
            timerange=trange,
            imagermode='csclean',
            psfmode='clark',
            imsize=[512,512],
            cell=['5arcsec'],
            niter=50,
            gain=0.05,
            stokes=pol,
            restoringbeam=['15arcsec'],
            interactive=False,
            pbcor=True,
            usescratch=True)

    hf.imreg(vis=slfcalms,imagefile=img_slfed+'.image',fitsfile=img_slfed+'.fits',
             timerange=trange,usephacenter=False,verbose=True)

    from sunpy import map as smap
    from matplotlib import pyplot as plt
    eomap=smap.Map(img_slfed+'.fits')
    eomap.data=eomap.data.reshape((512,512))
    eomap.plot_settings['cmap'] = plt.get_cmap('jet')
    eomap.plot()
    eomap.draw_limb()
    eomap.draw_grid()
    plt.show()

