import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp
import argparse

parser = argparse.ArgumentParser(description='Read an input cPickle file and plot the cross-correlation'\
                                 'between the antennas specified for all available channels (unless'\
                                 'otherwise spcified).',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('filename',type=str,
                    help= 'Name of the input CPickle file containing the cross-correlations')

parser.add_argument('-cc','--cross_corr', nargs=2, action='append', type=int, dest='cc',
                    help='Two antennas that should be cross-correlated')      

parser.add_argument('-chan','--channel', nargs='+', type=int, default=[0,1,2], dest='chan',
                    help='Channel to plot visibility for')

parser.add_argument('-plot',type=str,choices=['abs','phase','real','imag'],nargs='+',default='phase',
                    help='Complex numbers are multifaceted!')

# parser.add_argument('-gain','--gain', nargs=12, type=float, default=np.ones(12,dtype=float),
#                     help='Correction gain factors (need to specify all)')

args = parser.parse_args()

#ant_list = [84,85,86,87,52,53,54,55,24,25,26,27]
ant_list = [36,51,69,70,71,56,54,55,24,25,26,27] 
gains= {}
for chan in range(3):
    gains[chan] = {}
    for ant in range(12):
        gains[chan][ant_list[ant]] = 1
## Gains obtained for PAMs_ant12_2017-07-29.cp by eye-balling.
#gains[-1] = {84:1, 85:1, 86:np.exp(1j*np.pi),    87:np.exp(1j*np.pi/2),   52:1, 53:np.exp(1j*np.pi/2), 54:np.exp(1j*np.pi/2),   55:-1*np.exp(1j*np.pi/4), 24: 1, 25: 1, 26: np.exp(1j*np.pi/3),    27: -1*np.exp(1j*np.pi/6) }
#gains[1] = {84:1, 85:1, 86:np.exp(-1j*np.pi/6), 87:np.exp(1j*np.pi),     52:1, 53:np.exp(1j*np.pi/2), 54:np.exp(1j*2*np.pi/3), 55: np.exp(-1j*np.pi/4),  24: 1, 25: 1, 26:     -1,                27:        1              }
#gains[2] = {84:1, 85:1, 86:np.exp(-1j*np.pi/10),87:np.exp(-1j*np.pi/10), 52:1, 53:      1,            54:np.exp(-1j*np.pi/10), 55: 1 ,                   24: 1, 25: 1, 26: 1, 27: 1 }
#

with open(args.filename,'r') as fp:
    vis = cp.load(fp)

## (1) Retrive and conj if order is opposite
## (2) Calibrate
vis_corr = {}
for chan in args.chan:
    vis_corr[chan] = {}
    for ant1, ant2 in args.cc:
        try:
            data = vis['chan%d'%chan]['%dN-%dN'%(ant1,ant2)]
        except(KeyError):
            data = np.conj(vis['chan%d'%chan]['%dN-%dN'%(ant2,ant1)])
            
        vis_corr[chan]['%d-%d'%(ant1,ant2)] = gains[chan][ant1]*np.conj(gains[chan][ant2])*data
    
## (3) Plot
for plot_type in args.plot:
    if (plot_type == 'abs'):
        print ("Plotting absolute value of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Absolute value   Channel:%d'%chan)
            for ant1, ant2 in args.cc:
                plt.plot(np.abs(vis_corr[chan]['%d-%d'%(ant1,ant2)]),label='Ant%s-Ant%s'%(ant1,ant2))
            plt.legend()

    if (plot_type == 'phase'):
        print ("Plotting phase of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Phase   Channel:%d'%chan)
            for ant1,ant2 in args.cc:
                plt.plot(np.angle(vis_corr[chan]['%d-%d'%(ant1,ant2)]),label='Ant%s-Ant%s'%(ant1,ant2))
            plt.legend()

    if (plot_type == 'real' or plot_type == 'imag'):
        print ("Plotting %s part of visibilities"%plot_type)
        for chan in args.chan:
            plt.figure()
            plt.title('Channel:%d'%chan)
            if 'real' in args.plot:
                for ant1,ant2 in args.cc:
                    plt.plot(np.real(vis_corr[chan]['%d-%d'%(ant1,ant2)]),label='Ant%s-Ant%s Real'%(ant1,ant2))
                    
            if 'imag' in args.plot:
                for ant1,ant2 in args.cc:
                    plt.plot(np.imag(vis_corr[chan]['%d-%d'%(ant1,ant2)]),label='Ant%s-Ant%s Imag'%(ant1,ant2))
                    
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
