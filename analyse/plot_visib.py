import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp
import argparse
import sys

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

parser.add_argument('-jd','--julian-date',action='store_true',default=True,dest='jd',
                    help='Plot against Julian Day instead of integrations')


# Copied from redcal.py in heracal (https://github.com/HERA-Team/heracal/hercal/redcal.py)
def get_pos_reds(antpos, precisionFactor=1e6):
    """ Figure out and return list of lists of redundant baseline pairs. Ordered by length.
        All baselines have the same orientation with a preference for positive b_y and,
        when b_y==0, positive b_x where b((i,j)) = pos(i) - pos(j).

        Args:
            antpos: dictionary of antenna positions in the form {ant_index: np.array([x,y,z])}.
            precisionFactor: factor that when multiplied by different baseline vectors and rounded
                to integer values, gives unique integer tuples for unique baselines

        Returns:
            reds: list of lists of redundant tuples of antenna indices (no polarizations)
    """

    keys = antpos.keys()
    reds = {}
    array_is_2D = np.all(np.all(np.array(antpos.values())[:,2]==0))
    for i,ant1 in enumerate(keys):
        for ant2 in keys[i+1:]:
            delta = tuple((precisionFactor*2.0 * (np.array(antpos[ant1]) - np.array(antpos[ant2]))).astype(int))
            # Multiply by 2.0 because rounding errors can mimic changes below the grid spacing
            if delta[0] > 0 or (delta[0]==0 and delta[1] > 0) or (delta[0]==0 and delta[1]==0 and delta[2] > 0):
                bl_pair = (ant1,ant2)
            else:
                delta = tuple([-d for d in delta])
                bl_pair = (ant2,ant1)
            # Check to make sure reds doesn't have the key plus or minus rounding error
            p_or_m = (0,-1,1)
            if array_is_2D:
                epsilons = [[dx,dy,0] for dx in p_or_m for dy in p_or_m]
            else:
                epsilons = [[dx,dy,dz] for dx in p_or_m for dy in p_or_m for dz in p_or_m]
            for epsilon in epsilons:
                newKey = (delta[0]+epsilon[0], delta[1]+epsilon[1], delta[2]+epsilon[2])
                if reds.has_key(newKey):
                    reds[newKey].append(bl_pair)
                    break
            if not reds.has_key(newKey):
                reds[delta] = [bl_pair]
    orderedDeltas = [delta for (length,delta) in sorted(zip([np.linalg.norm(delta) for delta in reds.keys()],reds.keys()))]
    return [reds[delta] for delta in orderedDeltas]

def compute_correlation(ant_list):
    """ Computes the cross-correlation of the two antennas given and 
        corrects for the antenna gains as well.
        
        Args: 
            ant1:  int,    Antenna number (from HERA position map)
            ant2:  int,    Antenna number (cannot be same as ant1)
            gain1: float,  Gain of antenna1. default: 1
            gain2: float,  Gain of antenna2. default: 1
    
        Returns: 
            vis: complex array of cross-correlation amplitudes
    """
    tempvis = {}
    for chan in args.chan:
        tempvis[chan] = {}
        for ant1,ant2 in ant_list:
            if ant1==ant2:
                print ("Auto correlations not available in this version!")
                continue
            try:
                data = vis['chan%d'%chan]['%dN-%dN'%(ant1,ant2)]
            except(KeyError):
                data = np.conj(vis['chan%d'%chan]['%dN-%dN'%(ant2,ant1)])     

            tempvis[chan]['%d-%d'%(ant1,ant2)] = gains[chan][ant1]*np.conj(gains[chan][ant2])*data 

    return tempvis

def compute_jd(utc,day=10,month=10,year=2017,mjd=False):
    jd = 367*year -np.round(7*(year+np.round((month+9)/12))/4)\
         +np.round((275*month)/9) +day +1721013.5 +utc/24
    if mjd==False:
        return jd
    else:
        return jd-2400000.5 

args = parser.parse_args()


if args.ants_set.lower() == 'pams':
    ant_list = [84,85,86,87,52,53,54,55,24,25,26,27]
elif args.ants_set.lower() == 'nrao':
    ant_list = [36,51,69,70,71,56,54,55,24,25,26,27] 

vis,jds = None, []

for filename in args.files:
    with open(filename,'r') as fp:
        if vis==None:
            vis = cp.load(fp)
            v = vis
        else:
            v = cp.load(fp)
            for chan in vis.keys():
                for ants in vis[chan].keys():
                    vis[chan][ants] = np.append(vis[chan][ants],v[chan][ants]) 
        jdrange = []
        timerange = filename.split('/')[-1].split('_')[-1].rstrip('.cp').split('-')
        for t in timerange:
            jdrange.append(compute_jd(int(t.split(':')[0])+int(t.split(':')[1])/60. + 5))
        jds = np.append(jds,np.linspace(jdrange[0],jdrange[1],num=np.size(v['chan0']['84N-85N']),endpoint=False))

## (1) Retrive and conj if order is opposite
## (2) Calibrate
vis_corr = {}
gains = {}
for chan in args.chan:
    vis_corr[chan] = {}
    gains[chan] = {}
    for ant in ant_list:
        gains[chan][ant] = 1

if args.ant:
    ants = [(args.ant,a) for a in ant_list]
    ants.remove((args.ant,args.ant))
    visb_ants = compute_correlation(ants)
    for chan in args.chan:
        vis_corr[chan].update(visb_ants[chan])

if args.cc:
    visb_ants = compute_correlation(args.cc)        
    for chan in args.chan:
        vis_corr[chan].update(visb_ants[chan])

if args.bl:
    antpos = {24: [0,2,0], 25: [1,2,0], 26: [2,2,0], 27: [3,2,0],
              52: [0,1,0], 53: [1,1,0], 54: [2,1,0], 55: [3,1,0],
              84: [0,0,0], 85: [1,0,0], 86: [2,0,0], 87: [3,0,0]}
    redants = get_pos_reds(antpos)
    red_bls = [x for x in redants if tuple(args.cc[0]) in x]
    if not red_bls:
        red_bls = [[y[::-1] for y in [x for x in redants if tuple(args.cc[0][::-1]) in x][0]]]

    visb_ants = compute_correlation(red_bls[0])
    for chan in args.chan:
        vis_corr[chan].update(visb_ants[chan])

## (3) Plot
for plot_type in args.plot:
    if (plot_type == 'abs'):
        print ("Plotting absolute value of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Absolute value   Channel:%d'%chan)
            for ants in vis_corr[chan].keys():
                plt.plot(jds,np.log10(np.abs(vis_corr[chan][ants])),label=ants)
            plt.legend()

    if (plot_type == 'phase'):
        print ("Plotting phase of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Phase   Channel:%d'%chan)
            for ants in vis_corr[chan].keys():
                #plt.plot(jds,np.unwrap(np.angle(vis_corr[chan][ants][::-1]),discont=1.5*np.pi)[::-1],label=ants)
                plt.plot(jds,np.angle(vis_corr[chan][ants]),label=ants)
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
