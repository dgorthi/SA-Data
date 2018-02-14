import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import cPickle as cp
import os,sys, argparse, glob
import ephem, time
from astropy.time import Time
import aipy as ap
import plot

parser = argparse.ArgumentParser(description='Read an input cPickle file and plot the cross-correlation'\
                                 'between the antennas specified for all available channels (unless'\
                                 'otherwise spcified).',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('files',type=str, nargs='+',
                    help= 'Name of the input CPickle file containing the cross-correlations')
parser.add_argument('-ants_set',type=str,choices=['nrao','pams'],default='pams',
                    help= 'The set of antennas data was collected from') 
parser.add_argument('-ant',type=int,default=0,
                    help='All cross correlations for this antenna can be plotted')
parser.add_argument('-cc','--cross_corr', nargs=2, action='append', type=int, dest='cc',
                    help='Two antennas that should be cross-correlated')
parser.add_argument('-b', '--baseline', action='store_true', default=False, dest='bl',
                    help='Compute all cross correlations of the baseline specified by -cc flag')
parser.add_argument('-chan','--channel', nargs='+', type=int, default=[0,1,2], dest='chan',
                    help='Channel to plot visibility for')
parser.add_argument('-plot',type=str,choices=['abs','phase','real','imag'],nargs='+',default='phase',
                    help='Complex numbers are multifaceted!')
parser.add_argument('-gain','--gain', nargs=12, type=float, default=np.ones(12,dtype=float),
                    help='Correction gain factors (need to specify all \'em)')
parser.add_argument('-no_jd','--no_julian-date',action='store_false',default=True,dest='jd',
                    help='Plot against integrations instead of Julian day')
parser.add_argument('-nowrap','--unwrap_phase',action='store_true',default=False,dest='nowrap',
                    help='Unwrap phases when plotting phase of the visibilities')
parser.add_argument('-compare','--compare',action='store_true',default=False,
                    help='Compare with HERA correlator data of the next day')
args = parser.parse_args()

## Extract data from all files specified.

if args.ants_set.lower() == 'pams':
    ant_list = [84,85,86,87,52,53,54,55,24,25,26,27]
elif args.ants_set.lower() == 'nrao':
    ant_list = [36,51,69,70,71,56,54,55,24,25,26,27] 

vis,jds = None, []

for filename in args.files:
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

## (1) Retrive and conj if order is opposite, calibrate.
vis_corr = {}
gains = {}
for ant in ant_list:
    gains[ant] = [1,1,1]

if args.ant:
    ants = [(args.ant,a) for a in ant_list]
    ants.remove((args.ant,args.ant))
    visb_ants = plot.get_corr(vis,ants)
    vis_corr.update(visb_ants)

if args.cc:
    visb_ants = plot.get_corr(vis,args.cc)
    vis_corr.update(visb_ants)

if args.bl:
    antpos = {24: [0,2,0], 25: [1,2,0], 26: [2,2,0], 27: [3,2,0],
              52: [0,1,0], 53: [1,1,0], 54: [2,1,0], 55: [3,1,0],
              84: [0,0,0], 85: [1,0,0], 86: [2,0,0], 87: [3,0,0]}
    redants = plot.get_pos_reds(antpos)
    red_bls = [x for x in redants if tuple(args.cc[0]) in x]
    if not red_bls:
        red_bls = [[y[::-1] for y in [x for x in redants if tuple(args.cc[0][::-1]) in x][0]]]

    visb_ants = plot.get_corr(vis,red_bls[0])
    vis_corr.update(visb_ants)

if args.compare:
    print ("Retreiving correlator data...")
    files = glob.glob('/home/deepthi/Documents/HERA/Data/Oct11/oct11/2458038/zen.2458038.[1-3]*')
    files.sort()
    uv_chans = (np.array([804,624,584])-2)//4
    ants = ','.join(vis_corr[chan].keys()).replace('-','_')

    _uvlst,_uvdata = plot_uv_data(files,ants,'xx')

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
for plot_type in args.plot:
    if (plot_type == 'abs'):
        print ("Plotting absolute value of visibilities")
        fig['abs']={};ax['abs']={}
        for chan in args.chan:
            fig['abs'][chan],ax['abs'][chan] = plt.subplots(1,1)
            for ant in vis_corr.keys():
                amp = np.log10(np.abs(vis_corr[ant][:,chan]))           
                ax['abs'][chan].plot(sr,amp-amp.max(),label=ant)
                if args.compare:
                    try: amp = np.log10(np.abs(uvdata[ant][:,uv_chans[chan]]))
                    except KeyError:
                        ant = '-'.join(ant.split('-')[::-1])
                        amp = np.log10(np.abs(uvdata[ant][:,uv_chans[chan]]))
                    ax['abs'][chan].plot(uvlst,amp-amp.max(),label='HERA %s'%ant)
    
            fig['abs'][chan].suptitle('Absolute value   Channel:%d'%chan)
            ax['abs'][chan].legend()

    if (plot_type == 'phase'):
        print ("Plotting phase of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Phase   Channel:%d'%chan)
            if (args.nowrap==True):
                for ants in vis_corr[chan].keys():
                    plt.plot(jds,np.angle(vis_corr[chan][ants]),label=ants)
            else:
                for ants in vis_corr[chan].keys():
                    plt.plot(jds,np.unwrap(np.angle(vis_corr[chan][ants][::-1]),discont=1.5*np.pi)[::-1],label=ants)
            plt.legend()

    if (plot_type == 'real' or plot_type == 'imag'):
        print ("Plotting %s part of visibilities"%plot_type)
        for chan in args.chan:
            plt.figure()
            plt.title('Channel:%d'%chan)
            if 'real' in args.plot:
                for ants in vis_corr[chan].keys():
                    plt.plot(jds,np.real(vis_corr[chan][ants]),label=ants) 
            if 'imag' in args.plot:
                for ants in vis_corr[chan].keys():
                    plt.plot(jds,np.imag(vis_corr[chan][ants]),label=ants)
            plt.legend()

        try:
            args.plot.remove('real')
            args.plot.remove('imag')
        except(ValueError):
            pass
                    
plt.show()

#antenna_map = {0:'27N', 1:'27E', 2:'84N', 3:'84E'}
# antenna_map = {84: '84N', 85: '85N', 86: '86N', 87: '87N',
#                52: '52N', 53: '53N', 54: '54N', 55: '55N',
#                24: '24N', 25: '25N', 26: '26N', 27: '27N'  }



## Gains obtained for PAMs_ant12_2017-07-29.cp by eye-balling.
#gains[-1] = {84:1, 85:1, 86:np.exp(1j*np.pi),    87:np.exp(1j*np.pi/2),   52:1, 53:np.exp(1j*np.pi/2), 54:np.exp(1j*np.pi/2),   55:-1*np.exp(1j*np.pi/4), 24: 1, 25: 1, 26: np.exp(1j*np.pi/3),    27: -1*np.exp(1j*np.pi/6) }
#gains[1] = {84:1, 85:1, 86:np.exp(-1j*np.pi/6), 87:np.exp(1j*np.pi),     52:1, 53:np.exp(1j*np.pi/2), 54:np.exp(1j*2*np.pi/3), 55: np.exp(-1j*np.pi/4),  24: 1, 25: 1, 26:     -1,                27:        1              }
#gains[2] = {84:1, 85:1, 86:np.exp(-1j*np.pi/10),87:np.exp(-1j*np.pi/10), 52:1, 53:      1,            54:np.exp(-1j*np.pi/10), 55: 1 ,                   24: 1, 25: 1, 26: 1, 27: 1 }
#
