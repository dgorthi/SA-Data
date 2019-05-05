## Coarsely calibrate the visibilities by hand
import numpy as np
import matplotlib.pyplot as plt
import cPickle as cp
import argparse

def calibrate(data, gains):
    calib_data = {}
    for (i, j, pol) in data.keys():
        calib_data[(i,j,pol)] = data[(i,j,pol)]/(gains[(i,'Jxx')]*\
                                np.conj(gains[(j,'Jxx')]))
    return calib_data


parser = argparse.ArgumentParser(description='Coarse calibrate and compare with data',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

#parser.add_argument('datafile', type=str, 
#                     help ='cPickle data file with measured visibilities')
parser.add_argument('-c', '--chan', type=int, nargs='+', default=[1],
                     help ='Channel to plot for comparision.')
parser.add_argument('-cc','--cross_correlation',type=int, nargs='+',
                     help ='Pair of antennas to cross correlate.')
parser.add_argument('-nc','--no_calibration',action='store_true', default=False,
                     help='Plot the non calibrated data for comparision with model.')
parser.add_argument('-a','--amplitude', action='store_true', default=False,
                     help='Plot the absolute of the visibilities.')
parser.add_argument('-p','--phase', action='store_true', default=False,
                     help='Plot the phase of the visibilities.')
parser.add_argument('-r','--real', action='store_true', default=False,
                     help='Plot the real component of the visibilities.')
parser.add_argument('-i','--imag', action='store_true', default=False,
                     help='Plot the imaginary component of the visibilities.')

args = parser.parse_args()

chans = [55, 100, 110]
freqs = np.linspace(200, 100, endpoint=False, num=256)[chans]

## Firstcal gets the delay corrections (say due to cables etc) from 
## an FFT across the band. Since I have only three frequency channels,
## not continuously spaced, getting an accurate firstcal result from the 
## algorithm is not reasonable. Calibrating by hand to fix this problem.

# Compare visibility solutions with data
with open('data_11:39-17:30.cp','r') as fp:
    jds, reds, data = cp.load(fp)

g = {(24, 'Jxx'): [1,                  1, 1],
     (25, 'Jxx'): [1,                  1,1],
     (26, 'Jxx'): [1,                  1,1],
     (27, 'Jxx'): [1,                  1,1],
     (52, 'Jxx'): [np.exp(-2j*np.pi/4),1,1],
     (53, 'Jxx'): [1,                  1,1],
     (54, 'Jxx'): [1,                  1,1],
     (55, 'Jxx'): [1,                  1,1],
     (84, 'Jxx'): [1,                  1,1],
     (85, 'Jxx'): [1,                  1,1],
     (86, 'Jxx'): [1,                  1,1],
     (87, 'Jxx'): [1,                  1,1] }

if len(args.cross_correlation) > 2:
    print "Can plot only one pair because I am feeling dumb"

cc = args.cross_correlation
cc.append('xx')

bls = [subreds for subreds in reds if tuple(cc) in subreds][0]

fig, ax = {}, {}

if args.amplitude:
    print "Plotting absolute value"
    fig['abs']={}; ax['abs']={}
    for chan in args.chan:
        fig['abs'][chan], ax['abs'][chan] = plt.subplots(1,1)
        if args.no_calibration:
            for bl in bls:
                ax['abs'][chan].semilogy(jds, np.abs(data[bl][:,chan]), 
                                         label=bl, alpha=0.8)
        for (a0,a1,pol) in bls:
            d = np.abs(data[(a0,a1,pol)][:,chan]/(g[(a0,'Jxx')][chan]*np.conj(g[(a1,'Jxx')][chan]))) 
            ax['abs'][chan].semilogy(jds, d, label=(a0,a1,pol), alpha=0.8)

        fig['abs'][chan].suptitle('Absolute Value   Chan:%d'%chan)
        ax['abs'][chan].legend()


if args.phase:
    print ("Plotting phase of visibilities")
    fig['phase']={};ax['phase']={}
    for chan in args.chan:
        fig['phase'][chan],ax['phase'][chan] = plt.subplots(1,1)
        if args.no_calibration:
            for bl in bls:
                ax['phase'][chan].plot(jds, np.angle(data[bl][:,chan]), 
                                       label=bl, alpha=0.8)
        for (a0,a1,pol) in bls:
            d = data[(a0,a1,pol)][:,chan]/(g[(a0,'Jxx')][chan]*np.conj(g[(a1,'Jxx')][chan]))
            ax['phase'][chan].plot(jds, np.angle(d), label=(a0,a1,pol), alpha=0.8)

        fig['phase'][chan].suptitle('Phase   Channel:%d'%chan)
        ax['phase'][chan].legend()

if args.real or args.imag:
    print ("Plotting visibilities...")
    for chan in args.chan:
        plt.figure()
        plt.title('Channel:%d'%chan)
        if args.real:
            if args.no_calibration:
                for bl in bls:
                    plt.semilogy(jds, np.real(data[bl][:,chan]), label=bl, alpha=1)
            for (a0,a1,pol) in bls:
                d = np.real(data[(a0,a1,pol)][:,chan]/(g[(a0,'Jxx')][chan]*np.conj(g[(a1,'Jxx')][chan]))) 
                plt.semilogy(jds, d, label=(a0,a1,pol), alpha=0.8)

        if args.imag:
            if args.no_calibration:
                for bl in bls:
                    plt.semilogy(jds, np.imag(data[bl][:,chan]), label=bl, alpha=1)
            for (a0,a1,pol) in bls:
                d = np.imag(data[(a0,a1,pol)][:,chan]/(g[(a0,'Jxx')][chan]*np.conj(g[(a1,'Jxx')][chan]))) 
                plt.semilogy(jds, d, label=(a0,a1,pol), alpha=0.8)
        plt.legend()

plt.show()
