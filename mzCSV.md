# Introduction #
This is a very simple format for storing mass spec data. The conversion
to this format should be trivial from pretty much any other format.


# Data Fields #
# Scan Number
# Retention Time  in seconds
# Fragment m/z
# Intensity .. ie ion count
# Ms Level ( 1=full scan)
# Parent m/z ( precursorMz)
# Polarity
# SRM Identifier

```
scannum,rt,mz,intensity,mslevel,precursorMz,polarity,srmid
1,0.359,42.148384,4.529040,2,,-,-111.050 [41.600-42.600]
2,0.443,74.326462,4.264498,2,,-,-116.001 [73.800-74.800]
3,0.527,80.036446,4.307005,2,,-,-124.000 [79.500-80.500]
4,0.611,92.037384,4.119885,2,,-,-135.000 [91.500-92.500]
5,0.695,108.032822,4.102715,2,,-,-151.000 [107.500-108.500]
6,0.779,114.029968,4.125493,2,,-,-157.050 [113.500-114.500]
7,0.863,116.043579,4.153846,2,,-,-160.001 [115.500-116.500]

```

## About SRM Identifiers ##
On thermo instrument id **-111.050 [41.600-42.600]**  is decoded
to be  Polarity=negative, Parent m/z = 111.050, Fragment m/z = 41.600-42.600

