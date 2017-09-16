## Coarsely calibrate the visibilities by hand
## Why are the three channels showing different behaviour?

import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp

CHAN = 'chan0'

with open('../data/PAMs_ant12_2017-07-30.cp') as fp:
    vis = cp.load(fp)

gains = {}
gain[1] = 

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

args = parser.parse_args()

antenna_map = {0: '84N', 1: '85N', 2: '86N', 3: '87N',
               4: '52N', 5: '53N', 6: '54N', 7: '55N',
               8: '24N', 9: '25N', 10:'26N', 11:'27N'  }

with open(args.filename,'r') as fp:
    vis = cp.load(fp)

for plot_type in args.plot:
    if (plot_type == 'abs'):
        print ("Plotting absolute value of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Absolute value   Channel:%d'%chan)
            for ant1, ant2 in args.cc:
                plt.plot(np.abs(vis['chan%d'%chan]['%s-%s'%(antenna_map[ant1],antenna_map[ant2])]),label='Ant%s-Ant%s'%(antenna_map[ant1],antenna_map[ant2]))
            plt.legend()

    if (plot_type == 'phase'):
        print ("Plotting phase of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Phase   Channel:%d'%chan)
            for ant1,ant2 in args.cc:
                plt.plot(np.angle(vis['chan%d'%chan]['%s-%s'%(antenna_map[ant1],antenna_map[ant2])]),label='Ant%s-Ant%s'%(antenna_map[ant1],antenna_map[ant2]))
            plt.legend()

    if (plot_type == 'real' or plot_type == 'imag'):
        print ("Plotting %s part of visibilities")
        for chan in args.chan:
            plt.figure()
            plt.title('Channel:%d'%chan)
            if 'real' in args.plot:
                for ant1,ant2 in args.cc:
                    plt.plot(np.real(vis['chan%d'%chan]['%s-%s'%(antenna_map[ant1],antenna_map[ant2])]),label='Ant%s-Ant%s Real'%(antenna_map[ant1],antenna_map[ant2]))
                if 'imag' in args.plot:
                    plt.plot(np.imag(vis['chan%d'%chan]['%s-%s'%(antenna_map[ant1],antenna_map[ant2])]),label='Ant%s-Ant%s Real'%(antenna_map[ant1],antenna_map[ant2]))
            plt.legend()
        try:
            args.plot.remove('real')
            args.plot.remove('imag')
        except(ValueError):
            pass
                    
plt.show()



