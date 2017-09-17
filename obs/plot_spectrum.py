# Program to read a few files and plot the spectrum
import numpy as np
import argparse
import matplotlib.pyplot as plt
import struct
import h5py

UDP_DATA = 8064       # bytes
NACC     = 8192       # UDP pkts
CPY_BUF  = int((512*1024*1024.)/(NACC*UDP_DATA*2 + NACC)) #copy buffers per file (65)

parser = argparse.ArgumentParser(description='Plot visibility of two antennas for a single channel')

parser.add_argument('filename', help='Name of the file with data')
parser.add_argument('ant0', type=int, help='SNAP ADC Channel of first antenna (0-11)')
parser.add_argument('ant1', type=int, help='SNAP ADC Channel of second antenna (0-11)')
parser.add_argument('chan', type=int, help='Channel to analyse (0-2)')

args = parser.parse_args()

# filename = raw_input('Enter filename:')
# ant0 = int(raw_input('Enter antenna 1 (0-11): '))
# ant1 = int(raw_input('Enter antenna 2 (0-11): '))

# chan = int(raw_input('Channel num (0-2): '))

length_cpy = int(NACC*UDP_DATA*0.5)
struct_fmt = '<%dh%dB'%(length_cpy,NACC)           #one cpy buffer size
struct_len = struct.calcsize(struct_fmt)

## Antennas don't come out in order!!! Nor do channels!!!
## ant1, ant0, ant5, ant4, ant9, ant8 (all channels)
## ant3, ant2, ant7, ant6, ant11, ant10 (all channels)

map = {0: 2, 1: 0, 2: 38, 3: 36,
       4: 6, 5: 4, 6: 42, 7: 40,
       8:10, 9: 8,10: 46,11: 44 } 

a0 = np.array([]);a0_strt = map[args.ant0]+args.chan*12
a1 = np.array([]);a1_strt = map[args.ant1]+args.chan*12

with open(args.filename,'r') as f:
    for x in xrange(CPY_BUF):
        cbuf = struct.unpack(struct_fmt,f.read(struct_len))
        a0 = np.append(a0,1j*np.asarray(cbuf[a0_strt:length_cpy:72],dtype='float32') +np.asarray(cbuf[(a0_strt+1):length_cpy:72],dtype='float32'))
        a1 = np.append(a1,1j*np.asarray(cbuf[a1_strt:length_cpy:72],dtype='float32') +np.asarray(cbuf[(a1_strt+1):length_cpy:72],dtype='float32'))

corr = a0*a1.conj()

plt.plot(np.real(corr),label='real')
plt.plot(np.imag(corr),label='imag')
plt.legend()
plt.show()

#4128768

