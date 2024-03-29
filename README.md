# SA-Data
Code, Data and analysis software written for taking data from HERA antennas.                                                 

### Antenna array                                               
|               |               |               |               |
|:-------------:|:-------------:|:-------------:|:-------------:|
|       84      |       85      |       86      |       87      |                                              
|       52      |       53      |       54      |       55      |
|       24      |       25      |       26      |       27      |                                                             


### Observations:

Log of all observing runs while at South Africa between 17/7/2017 and 31/7/2017

|Date   |Receivers |Antennas  |Time UTC-5 (UTC+2) |Channels used (Fs)  |Name of hdf5 file  |Visibility file    |
|:------|:---------|:--------:|:------------------|:-------------------|:------------------|:------------------|
|7/24   |NRAO      |84(NF),85,86,87,55   |3.25-4:06:58 (10.25am-11.05am)   |0, 150, 170(?) (200MHz)  |test_fringes.hdf5 (wc) |testvis_chan0-3_ant0-5.cp |
|7/27   |PAMs      |84, 27 |11.28-12.13 (6.30pm to 7.15pm) |20, 50, 80 (250 MHz)  |test_fringes_withPAMs.hdf5 (wc) | test_fringes_withPAMs_7_28.cp |
|7/28   |PAMs      |84, 27   |03:51:14- 04:05:29 (~10am)  |35, 60, 81 (220 MHz)    | |  |
|7/28   |PAMs      |84, 27   |04:07:03- 04:20:10 (~10am)  |41, 92, 125 (220 MHz)   | |  |
|7/28   |PAMs      |84, 27   |04:22:49- 07:02:04 (~10am)  |46, 59, 76 (220 MHz)    |test_fringes_withPAMs_7_28.hdf5  |test_fringes_withPAMs_7_28.cp |
|7/28   |PAMs      |84-87,52-55,24-27 |28_13:33:31- 28_21:58:54 (8pm-5am) |46, 59, 76 (220 MHz)   |     |  |
|7/29   |NRAO      |36,51,69,70,71,56 |12:35:26- 18:30:46 (7pm-1am)       |51, 102, 204 (200 MHz) |NRAO_ant12_2017-07-29_overnight.hdf5  |NRAO_ant12_2017-07-29.cp |
|7/30   |NRAO      |36,51,69,70,71,56 |03:43:26- 05:39:20 |51, 102, 204 (200 MHz)  | | |
|7/30   |NRAO      |36,51,69,70,71,56 |05:40:22- 07:25:48 |51, 102, 210 (200 MHz)  | | |
|7/30   |PAMs      |84-87,52-55,24-27 |30_19:02:35- 31_01:43:32 (2am-8am) |46, 59, 76 (220 MHz) |PAMs_ant12_2017-07-30.hdf5 |PAMs_ant12_2017-07-30.cp |

### Finding redundant baselines with PAMS: 

#### 1B:                                  
These groups don't match each other.
84-85, 86-87
52-53, 53-54, 85-86         
54-55, 24-25   
25-26          
26-27
                
#### 1L(ong)B:                         
85-53, 86-54                                   
                                       
#### 2B:                                                                                                         
84-86, 85-87, 52-54 
