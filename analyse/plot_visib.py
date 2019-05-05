import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cPickle as cp
import os,sys, argparse, glob
import ephem, time
from astropy.time import Time, TimeDelta
import aipy as ap
import plot
import yaml

parser = argparse.ArgumentParser(description='Plot the cPickle visibility files and compare with the HERA correlator data',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('config', type=str, 
                     help ='YAML configuration file for plotting instructions.')
parser.add_argument('--hera', action='store_true', default=False,
                     help ='Plot HERA Correlator data as well for comparision.')
parser.add_argument('--no_data', action='store_true', default=False,
                     help='Do not plot your own data')
parser.add_argument('-c','--calibrate',action='store_true', default=False,
                     help='Calibrate with the gains in the config file')
parser.add_argument('-a','--amplitude', action='store_true', default=False,
                     help='Plot the absolute of the visibilities.')
parser.add_argument('-p','--phase', action='store_true', default=False,
                     help='Plot the phase of the visibilities.')
parser.add_argument('-r','--real', action='store_true', default=False,
                     help='Plot the real component of the visibilities.')
parser.add_argument('-i','--imag', action='store_true', default=False,
                     help='Plot the imaginary component of the visibilities.')

args = parser.parse_args()

with open(args.config,'r') as fp:
    conf = yaml.safe_load(fp)

print "Reading the config file"

strt_time = Time(1507649941, format='unix')
dt = TimeDelta(8192*8*56/(191e6/512.), format='sec')
#V['unix_time'] = strt_time + np.arange(totnsam,totnsam+nsam)*dt

ants = ['24','25','26','27','52','53','54','55','56','84','85','86','87']
if args.calibrate:
    gains = {}
    for k,v in conf['gains'].items():
        gains[k] = [eval(x) for x in v]
else:
    gains = {a:np.ones(3) for a in ants}     

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

#trange = Time(strt_time+np.arange(nsam[0])*dt, format='unix')
trange = Time(vis['unix_time'])
lst = trange.sidereal_time('apparent','21.443d')
sr = lst.radian   #sidereal radians

ants_mycorr = plot.get_ant_list(params)
vis_corr = plot.get_corr(vis,ants_mycorr)

if not np.any([args.amplitude, args.phase, args.real, args.imag]):
    print 'You have not asked me to plot anything!'
    raise ValueError

# (2) Extract HERACorr data
if args.hera:
    heraparams = conf['heracorr']
    
    print ("Retreiving correlator data...")
    files = glob.glob(heraparams['file_dir'])
    files.sort()
    uv_chans = (np.array([804,624,584])+2)//4
    #uv_chans = (np.array([721, 549, 511])-2)//4
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
if args.amplitude:
    print ("Plotting absolute value of visibilities")
    fig['abs']={};ax['abs']={}
    for chan in params['chan']:
        fig['abs'][chan],ax['abs'][chan] = plt.subplots(1,1)
        if not args.no_data:
            for pair in vis_corr.keys():
                ant0, ant1 = pair.split('-')
                amp = np.log10(np.abs(gains[ant0][chan]*np.conj(gains[ant1][chan])*vis_corr[pair][:,chan])) 
                ax['abs'][chan].plot(sr,amp-amp.max(),lw=2,label=pair)

        if args.hera:
            for antpair in ants_heracorr:
                if uvdata.has_key(antpair):
                    amp = np.log10(np.abs(uvdata[antpair][:,uv_chans[chan]]))
                else:
                    antpair = '_'.join(antpair.split('_')[::-1])
                    amp = np.log10(np.abs(uvdata[antpair][:,uv_chans[chan]]))
                ax['abs'][chan].plot(uvlst,amp-amp.max(),label='HERA %s'%antpair)

        fig['abs'][chan].suptitle('Absolute value   Channel:%d'%chan)
        ax['abs'][chan].legend()

if args.phase:
    print ("Plotting phase of visibilities")
    fig['phase']={};ax['phase']={}
    for chan in params['chan']:
        fig['phase'][chan],ax['phase'][chan] = plt.subplots(1,1)
        if not args.no_data:
            for pair in vis_corr.keys():
                ant0, ant1 = pair.split('-')
                ang = np.angle(gains[ant0][chan]*np.conj(gains[ant1][chan])*vis_corr[pair][:,chan])
                ax['phase'][chan].plot(sr,ang,lw=2,label=pair)

        if args.hera:
            for antpair in ants_heracorr:
                if uvdata.has_key(antpair):
                    ang = np.angle(uvdata[antpair][:,uv_chans[chan]])
                else:
                    antpair = '_'.join(antpair.split('_')[::-1])
                    ang = np.angle(uvdata[antpair][:,uv_chans[chan]])
                ax['phase'][chan].plot(uvlst,ang,label='HERA %s'%antpair)

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

if args.real or args.imag:
    print ("Plotting visibilities..")
    for chan in params['chan']:
        plt.figure()
        plt.title('Channel:%d'%chan)
        if args.real:
            for pair in vis_corr.keys():
                ant0, ant1 = pair.split('-')
                plt.plot(sr, np.real(gains[ant0][chan]*np.conj(gains[ant1][chan])*vis_corr[pair][:,chan]),label=pair)
        if args.imag:
            for pair in vis_corr.keys():
                ant0, ant1 = pair.split('-')
                plt.plot(sr, np.imag(gains[ant0][chan]*np.conj(gains[ant1][chan])*vis_corr[pair][:,chan]),label=pair)
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
