# Observing script to be run on the SNAP board after calibrating the
# ADCs with voltread_

import numpy as np
import corr.katcp_wrapper as k
import matplotlib.pyplot as plt
import struct
import math
import time
import sys

#Global variables
MAC_SNAP = 0xb827eba34523
IP_SNAP = 0x0a000ac7#((10<<24)+(0<<16)+(10<<8)+(199))
MASK = ((255<<24)+(255<<16)+(255<<8)+0)
MAC_COMP = 0x0002c94f0660
Fs       = 200                  # Sampling frequency- ADC clock

FREQ = np.linspace(Fs,Fs/2.,num=256)
ANT_LABELS = ['84','85','86','87','52','53','54','55','24','25','26','27']

if (len(sys.argv) == 4):
    CHAN1 = int(sys.argv[1])
    CHAN2 = int(sys.argv[2])
    CHAN3 = int(sys.argv[3])
else:
    print ('Input format: python obs_script.py <chan1> <chan2> <chan3>')
    sys.exit()

# Don't change this! Sync comes out to 2**17, 3**4, 5, 7
SYNC = math.factorial(9)*2**10
ACC_LEN = 2**16

def reset_10gbe():
    """Resets the 10Gb ethernet FIFO"""
    snap.write_int('valid_en',0)
    snap.write_int('rst',1)
    time.sleep(1)
    snap.write_int('rst',0)
    snap.write_int('valid_en',3)
    

def plot_chans(freq=True):
    """Plot all the inputs (you can turn off frequency axis with
    freq=False)"""
    f,ax = plt.subplots(4,3)
    for ant in range(12):
        snap.write_int('rst',1)
        snap.write_int('antenna',ant)
        snap.write_int('rst',0)

        time.sleep(ACC_LEN/(512*200e6)*1e3)
        arr = struct.unpack('>256Q',snap.read('spectrum',8*256))
        
        ax[ant%4][int(ant/4)].semilogy(FREQ,arr,'.-',lw=1)
        ax[ant%4][int(ant/4)].set_xlim(FREQ.max(), FREQ.min())
        ax[ant%4][int(ant/4)].set_title('Antenna %s'%ANT_LABELS[ant])
    plt.show()
    
def plotspec(ant,freq=True,test=False):
    """Plot spectrum of antenna.  

    Input: antenna number: Between 1-12
           test:           Bool. Plot test with signal

    Output: run plt.legend();plt.show() to see the graph.

    """
    snap.write_int('rst',1)
    snap.write_int('antenna',ant);
    snap.write_int('rst',0)

    ## Test
    if (test):
        if ant in [2,3,6,7,10,11]:
            test = (np.linspace(256,511,num=256,dtype='uint64')**2)*(ACC_LEN/256)
        else:
            test = (np.linspace(0,255,num=256,dtype='uint64')**2)*(ACC_LEN/256)
        plt.plot(test,'.',label='test')
    
    ## Wait for vector to get accumulated
    time.sleep(ACC_LEN/(512*200e6)*1e3)
    arr = struct.unpack('>256Q',snap.read('spectrum',8*256))
    if (freq):
        plt.plot(FREQ,arr,'.-',lw=1,label='%s'%ANT_LABELS[ant])
        plt.xlim(FREQ.max(), FREQ.min())
    else:
        plt.plot(arr,lw=2,label='%s'%ANT_LABELS[ant])

def write_ramp(bram):
    """Writes a frequency ramp to the test pfb 
    Input: bram: The bram to write the test vector to

    """
    snap.write_int('test',1)
    ramp = np.arange(2**9,dtype='uint64')
    snap.write(bram,struct.pack('>512L',*ramp))


print('Connecting to SNAP and programming..')
print('Make sure the ADC is calibrated in DEMUX 1 MODE, with a CLOCK FREQUENCY of 200MHz')
snap = k.FpgaClient('10.10.0.233')

time.sleep(0.1)

print('Connected!\n\nSetting registers..\nACC_LEN:\t%d\nSHIFT:\t2047\nDEST_IP:\t10.0.10.10\nDEST_PORT:\t10000\nSYNC PERIOD: %d'%(ACC_LEN,SYNC))

snap.write_int('acc_len',ACC_LEN)

snap.write_int('shift',2047)
snap.write_int('rst_of',1)
time.sleep(0.01)
snap.write_int('rst_of',0)

snap.write_int('dest_ip',((10<<24)+(0<<16)+(10<<8)+(10)))
snap.write_int('dest_port',10000)
snap.config_10gbe_core('xmit_sw', MAC_SNAP, IP_SNAP, 10001,[MAC_COMP for i in range(256)], gateway=0)

print('\nObserving frequency channels:\nCHANNEL\tFREQUENCY\n%d\t%f\n%d\t%f\n%d\t%f\n'%(CHAN1,FREQ[CHAN1],CHAN2,FREQ[CHAN2],CHAN3,FREQ[CHAN3]))
snap.write_int('chan1',CHAN1)
snap.write_int('chan2',CHAN2)
snap.write_int('chan3',CHAN3)

snap.write_int('test',0)

reset_10gbe()

snap.write_int('sync',SYNC)

