## Voltage data extracted into HDF5 file format is read and the
## visibility matrix is computed in this program. The data that the
## program runs on is stored in /data0 of the gbt-1u machine in digilab, Berkeley.
## Output: cPickle file with the visibility matrices of all the frequency channels in the observation.

import numpy as np
import cPickle as cp
from astropy.time import Time
import h5py
import bitshuffle.h5
import argparse
import sys,time
import glob

parser = argparse.ArgumentParser(description='Read an input hdf5 file with antenna voltages and compute visibilities from it.'\
                                 'Folding period computed from meta data in the hdf5 header. Specify either (-t, -f, -files) '\
                                 'or (-nsam,-files) in that combination to override default folding period'\
                                 'and integration period. Visibilities are output in a cPickle format.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('filename',type=str,
                    help= 'Name of the input hdf5 file. Compression technique does not matter.')
parser.add_argument('-log', '--log_output', action='store_true', default=False,
                    help='Log all output of this code')
parser.add_argument('-o', '--output', type=str,
                    help='Specify the cPickled output file name')
parser.add_argument('-ant', '--antennas', type=int,
                    help='Number of antennas in the array')
parser.add_argument('-chan','--channels',type=int, 
                    help='Number of channels to unpack')
parser.add_argument('-f','--sampling_freq', type=float, default= 200,
                    help='Sampling frequency [MHz] used to collect data')
args = parser.parse_args()

if(args.log_output):
    logfp = '../log/log_computevisb_%s.txt'%(args.filename.rstrip('.hdf5')[-11:])
    print ('Writing output to %s'%logfp)
    sys.stdout = open(logfp,'w',0)

V = {}
V['metadata'] = {}
meta = V['metadata']

with h5py.File(args.filename,'r') as fp:
    if not (args.antennas):  args.antennas = fp.attrs['Nants']
    if not (args.channels):  args.channels = fp.attrs['Nchans']
    # copy over metadata
    meta['num_files'] = fp.attrs['N_bin_files']
    meta['NACC'] = fp.attrs['NACC']
    meta['UDP'] = fp.attrs['UDP_B']
    meta['filesize'] = fp.attrs['bin_filesize_MB']
    meta['strt_file'] = fp.attrs['strt_file']
    meta['end_file'] = fp.attrs['end_file']
    meta['date'] = fp.attrs['date']
    fp.close()

## Print all the arguments to log file and to V
for i in vars(args):
    print i,vars(args)[i]
    meta[i] = vars(args)[i]
    
## Determine integration period
#total no of samples = samples per udp pkt * udp pkts per cpy buf * cpy bufs per file * total number of files
tot_sam = meta['UDP']/(4*meta['antennas']*meta['channels'])*meta['NACC']*int((meta['filesize']*1024*1024.)/(meta['NACC']*meta['UDP'] + meta['NACC']))*meta['num_files']
nsam = meta['num_files']
sample_rate = args.sampling_freq*1e6/512. #One sample per 512 clks is output acc. to the hardware design
foldlen = tot_sam/nsam

meta['nsam'] = nsam
meta['foldlen'] = foldlen
meta['integration_time'] = foldlen/sample_rate

print ('Setting NSAM=%d and FOLDLEN=%d'%(nsam,foldlen))

antenna_map = {0: '84N', 1: '85N', 2: '86N', 3: '87N',
               4: '52N', 5: '53N', 6: '54N', 7: '55N',
               8: '24N', 9: '25N', 10:'26N', 11:'27N'  }

# antenna_map = {0: '27E', 1: '27N', 2: '84E', 3: '84N'  }

def comp_vis(a1,a2):
    vis = np.zeros([nsam,3],dtype=np.complex64)
    print a1,a2
    with h5py.File(args.filename,'r') as fp:
        ant1 = [fp['%s_chan0'%a1],fp['%s_chan1'%a1],fp['%s_chan2'%a1]]
        ant2 = [fp['%s_chan0'%a2],fp['%s_chan1'%a2],fp['%s_chan2'%a2]]
        for i in range(nsam):
            print i
            s = i*foldlen; e = (i+1)*foldlen
            vis[i] = np.sum([c1[s:e]*np.conjugate(c2[s:e]) for c1,c2 in zip(ant1,ant2)],axis=1)
    fp.close()
    return vis

for ant1 in range(meta['antennas']):
    for ant2 in range(ant1,meta['antennas']):
        print('\n\nant%s * ant%s'%(antenna_map[ant1],antenna_map[ant2]))
        V['%s-%s'%(antenna_map[ant1],antenna_map[ant2])] = comp_vis(antenna_map[ant1],antenna_map[ant2])

print('Finished! Computing time array..')

# Compute unix time
file_list = filter(lambda fn: 'Oct10_PAMs' in fn, glob.glob('/data/dgorthi/*'))
file_list.sort()
idx = file_list.index(args.filename)

totnsam = 0
for i in range(idx):
    with h5py.File(file_list[i],'r') as fp:
        totnsam += fp.attrs['N_bin_files']
    fp.close()
 
dt = meta['integration_time']
strt_time = 1507649941 + totnsam*dt
t = np.arange(strt_time,strt_time+nsam*dt,step=dt)
V['unix_time'] = Time(t,format='unix')


with open(args.output,'wb') as fp:
    cp.dump(V,fp,protocol=2)
