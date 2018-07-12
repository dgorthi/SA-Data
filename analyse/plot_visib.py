import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cPickle as cp
import os,sys, argparse, glob
import ephem, time
from astropy.time import Time
import aipy as ap
import plot
import yaml

with open('config_plots.yaml','r') as fp:
    conf = yaml.load(fp)

## Process your data
params = conf['mycorr']
vis,jds = None, []

for filename in params['files']:
    with open(filename,'r') as fp:
        if vis==None:  vis = cp.load(fp)
        else:
            v = cp.load(fp)
            for k in vis.keys():
                if k == 'metadata': 
                    for k2 in vis['metadata'].keys():
                        vis['metadata'][k2] = np.append(vis[k][k2], v[k][k2])
                    continue
                vis[k] = np.concatenate((vis[k],v[k]),axis=0)

nsam = [len(vis[ant]) for ant in vis.keys() if not ant=='metadata']
assert(len(set(nsam)) == 1),'What the hell.'

trange = Time(vis['unix_time'])
lst = trange.sidereal_time('apparent','21.443d')
sr = lst.radian   #sidereal radians

ants_mycorr = plot.get_ant_list(params)
vis_corr = plot.get_corr(vis,ants_mycorr)

# (2) Extract HERACorr data
if conf['heracorr']:
    heraparams = conf['heracorr']
    
    print ("Retreiving correlator data...")
    files = glob.glob(heraparams['file_dir'])
    files.sort()
    uv_chans = (np.array([804,624,584])-2)//4
    _ants_heracorr = plot.get_ant_list(heraparams)
    ants_heracorr = ['%d_%d'%(ant1,ant2) for ant1,ant2 in _ants_heracorr]
    

    _uvlst,_uvdata = plot.get_uvdata(files,','.join(ants_heracorr),\
                                     heraparams['pol'])

    # Select only data in the range to be compared
    mask = np.where((_uvlst<sr.max()) & (_uvlst>sr.min()))
    uvlst = _uvlst[mask]
    uvdata = {}
    for ant in _uvdata.keys():
        uvtemp = ma.getdata(_uvdata[ant])[mask]
        #average over 4 channels
        uvdata[ant] = np.sum(uvtemp.reshape(-1,256,4),axis=2)

## (2) Plot
fig,ax={},{}
if conf['plot']['abs']:
    print ("Plotting absolute value of visibilities")
    fig['abs']={};ax['abs']={}
    for chan in params['chan']:
        fig['abs'][chan],ax['abs'][chan] = plt.subplots(1,1)
        for ant in vis_corr.keys():
            amp = np.log10(np.abs(vis_corr[ant][:,chan]))           
            ax['abs'][chan].plot(sr,amp-amp.max(),lw=2,label=ant)

        if conf['heracorr']:
            for antpair in ants_heracorr:
                if uvdata.has_key(antpair):
                    amp = np.log10(np.abs(uvdata[antpair][:,uv_chans[chan]]))
                else:
                    antpair = '_'.join(antpair.split('_')[::-1])
                    amp = np.log10(np.abs(uvdata[antpair][:,uv_chans[chan]]))
                ax['abs'][chan].plot(uvlst,amp-amp.max(),label='HERA %s'%antpair)

        fig['abs'][chan].suptitle('Absolute value   Channel:%d'%chan)
        ax['abs'][chan].legend()

if conf['plot']['phase']:
    print ("Plotting phase of visibilities")
    fig['phase']={};ax['phase']={}
    for chan in params['chan']:
        fig['phase'][chan],ax['phase'][chan] = plt.subplots(1,1)
        for ant in vis_corr.keys():
            ang = np.angle(vis_corr[ant][:,chan])
            ax['phase'][chan].plot(sr,ang,lw=2,label=ant)

        if conf['heracorr']:
            for antpair in ants_heracorr:
                if uvdata.has_key(antpair):
                    ang = np.angle(uvdata[antpair][:,uv_chans[chan]])
                else:
                    antpair = '_'.join(antpair.split('_')[::-1])
                    ang = np.angle(uvdata[antpair][:,uv_chans[chan]])
                ax['abs'][chan].plot(uvlst,ang,label='HERA %s'%antpair)

        fig['phase'][chan].suptitle('Phase   Channel:%d'%chan)
        ax['phase'][chan].legend()
#    for chan in params['chan']:
#        plt.figure()
#        plt.title('Phase   Channel:%d'%chan)
#        if (params['nowrap']==True):
#            for ants in vis_corr[chan].keys():
#                plt.plot(jds,np.angle(vis_corr[chan][ants]),label=ants)
#        else:
#            for ants in vis_corr[chan].keys():
#                plt.plot(jds,np.unwrap(np.angle(vis_corr[chan][ants][::-1]),discont=1.5*np.pi)[::-1],label=ants)
#        plt.legend()

if conf['plot']['real'] or conf['plot']['imag']:
    print ("Plotting %s part of visibilities"%plot_type)
    for chan in params['chan']:
        plt.figure()
        plt.title('Channel:%d'%chan)
        if conf['plot']['real']:
            for ants in vis_corr[chan].keys():
                plt.plot(jds,np.real(vis_corr[chan][ants]),label=ants) 
        if conf['plot']['imag']:
            for ants in vis_corr[chan].keys():
                plt.plot(jds,np.imag(vis_corr[chan][ants]),label=ants)
        plt.legend()

plt.show()

## Gains obtained for PAMs_ant12_2017-07-29.cp by eye-balling.
#gains[-1] = {84:1,  85:1,                  86:np.exp(1j*np.pi),    87:np.exp(1j*np.pi/2),   
#             52:1,  53:np.exp(1j*np.pi/2), 54:np.exp(1j*np.pi/2),  55:-1*np.exp(1j*np.pi/4), 
#             24: 1, 25: 1,                 26: np.exp(1j*np.pi/3), 27: -1*np.exp(1j*np.pi/6) }
#gains[1] =  {84:1,  85:1,                  86:np.exp(-1j*np.pi/6), 87:np.exp(1j*np.pi),     
#             52:1,  53:np.exp(1j*np.pi/2), 54:np.exp(1j*2*np.pi/3),55: np.exp(-1j*np.pi/4),  
#             24:1,  25:1,                  26:     -1,             27:        1              }
#gains[2] = { 84:1,  85:1,                  86:np.exp(-1j*np.pi/10),87:np.exp(-1j*np.pi/10), 
#             52:1,  53:1,                  54:np.exp(-1j*np.pi/10),55: 1 ,                   
#             24:1,  25:1,                  26: 1,                  27: 1                     }
#
