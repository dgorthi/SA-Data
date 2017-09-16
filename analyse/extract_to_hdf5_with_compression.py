import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py,bitshuffle.h5
import struct,glob,os,time,sys

parser = argparse.ArgumentParser(description='Extract binary data from multiple files into a hdf5 file with bit shuffle compression',
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('strt_filename',type=str,
                    help= 'Name of the first file to analyse')
parser.add_argument('stop_filename',type=str,
                    help= 'Name of the last file to analyse')

parser.add_argument('-o', '--output', type=str, default= 'Output_%s'%(time.strftime('%d-%m-%Y_%H:%M:%S',time.localtime())),
                    help='Specify the output file name')
parser.add_argument('-ant', '--antennas', type=int, default= 12,
                    help='Number of antennas in the array')
parser.add_argument('-chan','--channels',type=int, default=3,
                    help='Number of channels to unpack')
parser.add_argument('--NACC',type=int, default= 8192,
                    help='Number of UDP packets per buffer')
parser.add_argument('--UDP', type=int, default= 8064,
                    help='Size (in bytes) of each UDP packet')
parser.add_argument('--filesize', type=int, default= 512,
                    help='Size of each input file in MB')
parser.add_argument('-log', '--log_output', action='store_true', default=False,
                    help='Log all output of this code')

args = parser.parse_args()

if(args.log_output):
    sys.stdout = open('log_unpacking_%s.txt'%(time.strftime('%d-%m-%Y',time.localtime())),'w',0)

## Check if all the files are in the same directory
assert(args.strt_filename.split('/')[1] == args.stop_filename.split('/')[1])
directory = args.strt_filename.split('/')[1]
date = args.strt_filename.split('_')[1]

files = sorted(glob.glob('/'+directory+'/data_'+date[:8]+'*'), key=os.path.getmtime)
strt = files.index(args.strt_filename)
stop = files.index(args.stop_filename)
files = files[strt:stop]

## Starting points of data for each antenna- based on the Simulink
## design data structure.
file_map = {0: 2,   1: 0, 2: 38,  3: 36,
            4: 6,   5: 4, 6: 42,  7: 40,
            8: 10,  9: 8, 10: 46, 11: 44 }

## Labels for antennas- numbers correspond to ADCs on SNAP
antenna_map = {0: '84N', 1: '85N', 2: '86N', 3: '87N',
               4: '52N', 5: '53N', 6: '54N', 7: '55N',
               8: '24N', 9: '25N', 10:'26N', 11:'27N'  }

CPY_BUF  = int((args.filesize*1024*1024.)/(args.NACC*args.UDP + args.NACC))  # copy buffers per file

length_cpy = int(args.NACC*args.UDP*0.5)      # one cpy buffer size, 0.5 to account for short ints
struct_fmt = '<%dh%dB'%(length_cpy,args.NACC)
struct_len = struct.calcsize(struct_fmt)
nsam = int(length_cpy/(args.antennas*args.channels*2))

print ("Creating a global file for %d antennas and %d channels"%(args.antennas,args.channels))
print ("Data from %s"%date)
print ("Unpacking data from file: %s\n To file: %s"%(args.strt_filename,args.stop_filename))
print ("Output hdf5 file: %s"%('/'+directory+'/'+args.output))
print ("Compressed using the bitshuffle algorithm")

volts = {}

with h5py.File('/'+directory+'/'+args.output,'a') as fp: 
    ## Create datasets
    block_size=0
    for ant in range(args.antennas):
        volts[ant] = {}
        for chan in range(args.channels):
            volts[ant][chan] = fp.create_dataset('%s_chan%d'%(antenna_map[ant],chan),(0,),dtype='complex64',maxshape=(None,),chunks=(26091520,),compression=bitshuffle.h5.H5FILTER,compression_opts=(block_size, bitshuffle.h5.H5_COMPRESS_LZ4))

    #TODO: Delete the missing packets 8064 bytes (4032 h) at that site
    # missed = np.count_nonzero(cbuf[-1024:])
    # udp = 1024-missed
    for filename in files:
        if (os.path.getsize(filename)!=528547840):
            print("\nSkipping file %s\n"%filename)
            continue
        print ("Reading " +filename + " ...")
        with open(filename,'r') as f:
        
            for x in xrange(CPY_BUF):
                print(x)
                cbuf = struct.unpack(struct_fmt,f.read(struct_len))

                for chan in range(args.channels):
                    for ant in range(args.antennas):
                        volts[ant][chan].resize(volts[ant][chan].shape[0]+nsam,axis=0)
                        volts[ant][chan][-nsam:] = 1j*np.asarray(cbuf[(file_map[ant]+12*chan):length_cpy:72],dtype='float32') +np.asarray(cbuf[(file_map[ant]+12*chan+1):length_cpy:72],dtype='float32')
                        
print("Finished unpacking all the data.. You may read the file now")
