## Voltage data extracted into HDF5 file format is read and the
## visibility matrix is computed in this program. The data that the
## program runs on is stored in /data0 of the gbt-1u machine in digilab, Berkeley.
## Output: cPickle file with the visibility matrices of all the frequency channels in the observation.

import numpy as np
import cPickle as cp
import h5py
import argparse
import sys,time

parser = argparse.ArgumentParser(description='Read an input hdf5 file with antenna voltages and compute visibilities from it.'\
                                 'Specify either (-t, -f, -files) or (-nsam,-files) in that combination to determine folding period'\
                                 'and integration period. Visibilities are output in a cPickle format.',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('filename',type=str,
                    help= 'Name of the input hdf5 file. Compression technique does not matter.')
parser.add_argument('-log', '--log_output', action='store_true', default=False,
                    help='Log all output of this code')
parser.add_argument('-o', '--output', type=str, default= 'Output_%s.cp'%(time.strftime('%d-%m-%Y_%H:%M:%S',time.localtime())),
                    help='Specify the cPickled output file name')
parser.add_argument('-ant', '--antennas', type=int, default= 12,
                    help='Number of antennas in the array')
parser.add_argument('-chan','--channels',type=int, default=3,
                    help='Number of channels to unpack')
parser.add_argument('-N','--nsam',type=int, default=0,
                    help='Number of time samples user wants in the output. This sets the integration time.')                    
parser.add_argument('-t','--integration_time', type=float, default=0,
                    help='Integration time in seconds')
parser.add_argument('-f','--sampling_freq', type=float, default=0,
                    help='Sampling frequency [MHz] used to collect data')
parser.add_argument('-files','--num_files',type=int,default=0,
                    help='Total number of binary files unpacked to create input hdf5 file')
parser.add_argument('--NACC',type=int, default= 8192,
                    help='Number of UDP packets per buffer')
parser.add_argument('--UDP', type=int, default= 8064,
                    help='Size (in bytes) of each UDP packet')
parser.add_argument('--filesize', type=int, default= 512,
                    help='Size of each input file in MB')
args = parser.parse_args()

if(args.log_output):
    sys.stdout = open('log_computevisb_%s.txt'%(time.strftime('%d-%m-%Y',time.localtime())),'w',0)

## Print all the arguments to log file
for i in vars(args):
    print i,vars(args)['%s'%i]
    
## Determine integration period
#total no of samples = samples per udp pkt * udp pkts per cpy buf * cpy bufs per file * total number of files
tot_sam = args.UDP/(4*args.antennas*args.channels) * args.NACC * int((args.filesize*1024*1024.)/(args.NACC*args.UDP + args.NACC)) * args.num_files

if(args.integration_time):
    if (args.nsam):
        print('Overriding the nsam=%d input. Set sampling frequency to zero to use that instead'%args.nsam)
    if not(args.sampling_freq) or not(args.num_files):
        sys.exit('Required input: -t <integration time in sec> -f <sampling frequency in MHz> -files <Number of files used to create the input hdf5 file>')
    
    sample_rate = args.sampling_freq*1e6/512. #One sample per 512 clks is output acc. to the hardware design
    foldlen = int(args.integration_time*sample_rate)
    nsam = int(tot_sam/foldlen)
else:
    if not(args.nsam) or not(args.num_files):
        sys.exit('Specify either (-t, -f, -files) or (-nsam,-files) in that combination to determine folding length and integration period.')
    nsam = args.nsam
    foldlen = int(tot_sam/nsam)
print ('Setting NSAM=%d and FOLDLEN=%d'%(nsam,foldlen))


V = {}
antenna_map = {0: '84N', 1: '85N', 2: '86N', 3: '87N',
               4: '52N', 5: '53N', 6: '54N', 7: '55N',
               8: '24N', 9: '25N', 10:'26N', 11:'27N'  }

# antenna_map = {0: '27E', 1: '27N', 2: '84E', 3: '84N'  }

def comp_vis(a1,a2):
    vis = np.zeros([nsam],dtype=np.complex64)
    for i in range(nsam):
        print(i)
        vis[i] = np.sum(a1[i*foldlen:(i+1)*foldlen]*np.conjugate(a2[i*foldlen:(i+1)*foldlen]))
    return vis

with h5py.File(args.filename,'r') as fp:
    for chan in range(args.channels):
        V['chan%d'%chan] = {}
        for ant1 in range(args.antennas):
            for ant2 in range(ant1+1,args.antennas):
                print('chan=%d\tant%s * ant%s'%(chan,antenna_map[ant1],antenna_map[ant2]))
                V['chan%d'%chan]['%s-%s'%(antenna_map[ant1],antenna_map[ant2])] = comp_vis(fp['%s_chan%d'%(antenna_map[ant1],chan)],fp['%s_chan%d'%(antenna_map[ant2],chan)])
  
print('Finished! Writing data to file..')

with open(args.output,'wb') as fp:
    cp.dump(V,fp,protocol=2)
