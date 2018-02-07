## Voltage data extracted into HDF5 file format is read and the
## visibility matrix is computed in this program. The data that the
## program runs on is stored in /data0 of the gbt-1u machine in digilab, Berkeley.
## Output: cPickle file with the visibility matrices of all the frequency channels in the observation.

import numpy as np
import cPickle as cp
import h5py
import bitshuffle.h5
import argparse
import sys,time

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
    logfp = '../log/log_computevisb_%s.txt'%(time.strftime('%d-%m-%Y',time.localtime()))
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
    vis = np.zeros([3,nsam],dtype=np.complex64)
    for i in range(nsam):
        print(i)
        vis[i] = np.sum(a1[i*foldlen:(i+1)*foldlen]*np.conjugate(a2[i*foldlen:(i+1)*foldlen]))
    return vis

with h5py.File(args.filename,'r') as fp:
    for ant1 in range(meta['antennas']):
        for ant2 in range(ant1,meta['antennas']):
            print('\n\nant%s * ant%s'%(antenna_map[ant1],antenna_map[ant2]))
            V['%s-%s'%(antenna_map[ant1],antenna_map[ant2])] = comp_vis(fp['%s_chan%d'%(antenna_map[ant1],chan)],fp['%s_chan%d'%(antenna_map[ant2],chan)])                
            for chan in range(meta['channels']):
                v.append()
    fp.close()

print('Finished! Writing output to file..')

with open(args.output,'wb') as fp:
    cp.dump(V,fp,protocol=2)
