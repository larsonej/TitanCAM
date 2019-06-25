Initialization for Ames Aerosol Model (Version 1.01)

Doing a cold start initialization for a new simulation

Particle grid structure (setupaer):
ngroups:       2
nelem:       3
nbins:      40
Massmin:     0.420E-20   0.420E-20
Mrat:          2.00        2.00
nelemg:       1     2
ntype:       0     1     2
ngpart:       1     2
ngelem:       1     2     2
 
Coagulation bin mapping arrays
igrp, ibin =   1  1
igrp, ibin =   1  2
low:np,ig,jg,i,j   1  1  1  1  1
low:np,ig,jg,i,j   2  1  1  1  2
igrp, ibin =   1  3
low:np,ig,jg,i,j   1  1  1  1  3
low:np,ig,jg,i,j   2  1  1  2  2
low:np,ig,jg,i,j   3  1  1  2  3
up:np,ig,jg,i,j   1  1  1  1  1
up:np,ig,jg,i,j   2  1  1  1  2
up:np,ig,jg,i,j   3  1  1  2  1
igrp, ibin =   2  1
igrp, ibin =   2  2
low:np,ig,jg,i,j   1  1  2  1  1
low:np,ig,jg,i,j   2  1  2  1  2
low:np,ig,jg,i,j   3  1  2  2  1
low:np,ig,jg,i,j   4  2  1  1  1
low:np,ig,jg,i,j   5  2  1  1  2
low:np,ig,jg,i,j   6  2  2  1  1
low:np,ig,jg,i,j   7  2  2  1  2
igrp, ibin =   2  3
low:np,ig,jg,i,j   1  1  2  1  3
low:np,ig,jg,i,j   2  1  2  2  2
low:np,ig,jg,i,j   3  1  2  2  3
low:np,ig,jg,i,j   4  1  2  3  1
low:np,ig,jg,i,j   5  1  2  3  2
low:np,ig,jg,i,j   6  2  1  1  3
low:np,ig,jg,i,j   7  2  1  2  2
low:np,ig,jg,i,j   8  2  1  2  3
low:np,ig,jg,i,j   9  2  2  1  3
low:np,ig,jg,i,j   10  2  2  2  2
low:np,ig,jg,i,j   11  2  2  2  3
up:np,ig,jg,i,j   1  1  2  1  1
up:np,ig,jg,i,j   2  1  2  1  2
up:np,ig,jg,i,j   3  1  2  2  1
up:np,ig,jg,i,j   4  2  1  1  1
up:np,ig,jg,i,j   5  2  1  1  2
up:np,ig,jg,i,j   6  2  1  2  1
up:np,ig,jg,i,j   7  2  2  1  1
up:np,ig,jg,i,j   8  2  2  1  2
up:np,ig,jg,i,j   9  2  2  2  1

Initial size distributions (initaer):

   Ielem     ibin        r            N(r)
     1        1       0.10E-02       0.29E-02
     1        2       0.13E-02       0.14E-01
     1        3       0.16E-02       0.62E-01
     1        4       0.20E-02       0.24E+00
     1        5       0.25E-02       0.84E+00
     1        6       0.32E-02       0.26E+01
     1        7       0.40E-02       0.73E+01
     1        8       0.50E-02       0.18E+02
     1        9       0.64E-02       0.40E+02
     1       10       0.80E-02       0.80E+02
     1       11       0.10E-01       0.14E+03
     1       12       0.13E-01       0.23E+03
     1       13       0.16E-01       0.33E+03
     1       14       0.20E-01       0.42E+03
     1       15       0.25E-01       0.48E+03
     1       16       0.32E-01       0.49E+03
     1       17       0.40E-01       0.45E+03
     1       18       0.51E-01       0.37E+03
     1       19       0.64E-01       0.27E+03
     1       20       0.81E-01       0.18E+03
     1       21       0.10E+00       0.10E+03
     1       22       0.13E+00       0.55E+02
     1       23       0.16E+00       0.26E+02
     1       24       0.20E+00       0.11E+02
     1       25       0.26E+00       0.41E+01
     1       26       0.32E+00       0.14E+01
     1       27       0.41E+00       0.42E+00
     1       28       0.51E+00       0.11E+00
     1       29       0.65E+00       0.27E-01
     1       30       0.81E+00       0.59E-02
     1       31       0.10E+01       0.11E-02
     1       32       0.13E+01       0.20E-03
     1       33       0.16E+01       0.31E-04
     1       34       0.20E+01       0.43E-05
     1       35       0.26E+01       0.53E-06
     1       36       0.33E+01       0.59E-07
     1       37       0.41E+01       0.58E-08
     1       38       0.52E+01       0.52E-09
     1       39       0.65E+01       0.41E-10
     1       40       0.82E+01       0.29E-11
     2        1       0.87E-03       0.41E-05
     2        2       0.11E-02       0.30E-04
     2        3       0.14E-02       0.19E-03
     2        4       0.17E-02       0.11E-02
     2        5       0.22E-02       0.57E-02
     2        6       0.28E-02       0.27E-01
     2        7       0.35E-02       0.11E+00
     2        8       0.44E-02       0.41E+00
     2        9       0.56E-02       0.14E+01
     2       10       0.70E-02       0.40E+01
     2       11       0.88E-02       0.11E+02
     2       12       0.11E-01       0.26E+02
     2       13       0.14E-01       0.54E+02
     2       14       0.18E-01       0.10E+03
     2       15       0.22E-01       0.18E+03
     2       16       0.28E-01       0.27E+03
     2       17       0.35E-01       0.37E+03
     2       18       0.44E-01       0.45E+03
     2       19       0.56E-01       0.49E+03
     2       20       0.71E-01       0.48E+03
     2       21       0.89E-01       0.42E+03
     2       22       0.11E+00       0.33E+03
     2       23       0.14E+00       0.23E+03
     2       24       0.18E+00       0.14E+03
     2       25       0.22E+00       0.81E+02
     2       26       0.28E+00       0.41E+02
     2       27       0.36E+00       0.18E+02
     2       28       0.45E+00       0.74E+01
     2       29       0.56E+00       0.27E+01
     2       30       0.71E+00       0.86E+00
     2       31       0.90E+00       0.25E+00
     2       32       0.11E+01       0.64E-01
     2       33       0.14E+01       0.15E-01
     2       34       0.18E+01       0.30E-02
     2       35       0.23E+01       0.56E-03
     2       36       0.28E+01       0.93E-04
     2       37       0.36E+01       0.14E-04
     2       38       0.45E+01       0.18E-05
     2       39       0.57E+01       0.21E-06
     2       40       0.72E+01       0.23E-07
     3        1       0.87E-03       0.17E-26
     3        2       0.11E-02       0.25E-25
     3        3       0.14E-02       0.32E-24
     3        4       0.17E-02       0.37E-23
     3        5       0.22E-02       0.39E-22
     3        6       0.28E-02       0.36E-21
     3        7       0.35E-02       0.30E-20
     3        8       0.44E-02       0.22E-19
     3        9       0.56E-02       0.15E-18
     3       10       0.70E-02       0.87E-18
     3       11       0.88E-02       0.46E-17
     3       12       0.11E-01       0.22E-16
     3       13       0.14E-01       0.94E-16
     3       14       0.18E-01       0.36E-15
     3       15       0.22E-01       0.12E-14
     3       16       0.28E-01       0.37E-14
     3       17       0.35E-01       0.10E-13
     3       18       0.44E-01       0.25E-13
     3       19       0.56E-01       0.54E-13
     3       20       0.71E-01       0.11E-12
     3       21       0.89E-01       0.19E-12
     3       22       0.11E+00       0.29E-12
     3       23       0.14E+00       0.41E-12
     3       24       0.18E+00       0.51E-12
     3       25       0.22E+00       0.57E-12
     3       26       0.28E+00       0.58E-12
     3       27       0.36E+00       0.52E-12
     3       28       0.45E+00       0.42E-12
     3       29       0.56E+00       0.30E-12
     3       30       0.71E+00       0.19E-12
     3       31       0.90E+00       0.11E-12
     3       32       0.11E+01       0.57E-13
     3       33       0.14E+01       0.26E-13
     3       34       0.18E+01       0.11E-13
     3       35       0.23E+01       0.40E-14
     3       36       0.28E+01       0.13E-14
     3       37       0.36E+01       0.40E-15
     3       38       0.45E+01       0.10E-15
     3       39       0.57E+01       0.25E-16
     3       40       0.72E+01       0.53E-17

Model will run with the following values:
ibtime:       0
ietime:     100
nprint:       1
nprint:       1
nhist:       1
nrest:    5000
nx:       1
ny:       1
nz:       1
nbins:      40
nelem:       3
time:          0.00
itime:       0
dtime:          0.00
simtitle:  A SIMPLE COAGULATION SIMULATION

Timestep itime:      0   time:         0.00

Total mass, CN mass:   0.42000000E-20  0.43802621E-11

History output #      1 at itime:      0   time:         0.00

Timestep itime:      1   time:         0.00

Total mass, CN mass:   0.42000000E-20  0.43802619E-11

History output #      2 at itime:      1   time:         0.00

Timestep itime:      2   time:         0.01

Total mass, CN mass:   0.42000000E-20  0.43802617E-11

History output #      3 at itime:      2   time:         0.01

Timestep itime:      3   time:         0.01

Total mass, CN mass:   0.42000000E-20  0.43802616E-11

History output #      4 at itime:      3   time:         0.01

Timestep itime:      4   time:         0.02

Total mass, CN mass:   0.42000000E-20  0.43802614E-11

History output #      5 at itime:      4   time:         0.02

Timestep itime:      5   time:         0.02

Total mass, CN mass:   0.42000000E-20  0.43802613E-11

History output #      6 at itime:      5   time:         0.02

Timestep itime:      6   time:         0.02

Total mass, CN mass:   0.42000000E-20  0.43802612E-11

History output #      7 at itime:      6   time:         0.02

Timestep itime:      7   time:         0.03

Total mass, CN mass:   0.42000000E-20  0.43802611E-11

History output #      8 at itime:      7   time:         0.03

Timestep itime:      8   time:         0.03

Total mass, CN mass:   0.42000000E-20  0.43802609E-11

History output #      9 at itime:      8   time:         0.03

Timestep itime:      9   time:         0.04

Total mass, CN mass:   0.42000000E-20  0.43802608E-11

History output #     10 at itime:      9   time:         0.04

Timestep itime:     10   time:         0.04

Total mass, CN mass:   0.42000000E-20  0.43802607E-11

History output #     11 at itime:     10   time:         0.04

Timestep itime:     11   time:         0.04

Total mass, CN mass:   0.42000000E-20  0.43802606E-11

History output #     12 at itime:     11   time:         0.04

Timestep itime:     12   time:         0.05

Total mass, CN mass:   0.42000000E-20  0.43802605E-11

History output #     13 at itime:     12   time:         0.05

Timestep itime:     13   time:         0.05

Total mass, CN mass:   0.42000000E-20  0.43802604E-11

History output #     14 at itime:     13   time:         0.05

Timestep itime:     14   time:         0.06

Total mass, CN mass:   0.42000000E-20  0.43802603E-11

History output #     15 at itime:     14   time:         0.06

Timestep itime:     15   time:         0.06

Total mass, CN mass:   0.42000000E-20  0.43802602E-11

History output #     16 at itime:     15   time:         0.06

Timestep itime:     16   time:         0.06

Total mass, CN mass:   0.42000000E-20  0.43802601E-11

History output #     17 at itime:     16   time:         0.06

Timestep itime:     17   time:         0.07

Total mass, CN mass:   0.42000000E-20  0.43802600E-11

History output #     18 at itime:     17   time:         0.07

Timestep itime:     18   time:         0.07

Total mass, CN mass:   0.42000000E-20  0.43802600E-11

History output #     19 at itime:     18   time:         0.07

Timestep itime:     19   time:         0.08

Total mass, CN mass:   0.42000000E-20  0.43802599E-11

History output #     20 at itime:     19   time:         0.08

Timestep itime:     20   time:         0.08

Total mass, CN mass:   0.42000000E-20  0.43802598E-11

History output #     21 at itime:     20   time:         0.08

Timestep itime:     21   time:         0.08

Total mass, CN mass:   0.42000000E-20  0.43802597E-11

History output #     22 at itime:     21   time:         0.08

Timestep itime:     22   time:         0.09

Total mass, CN mass:   0.42000000E-20  0.43802596E-11

History output #     23 at itime:     22   time:         0.09

Timestep itime:     23   time:         0.09

Total mass, CN mass:   0.42000000E-20  0.43802596E-11

History output #     24 at itime:     23   time:         0.09

Timestep itime:     24   time:         0.10

Total mass, CN mass:   0.42000000E-20  0.43802595E-11

History output #     25 at itime:     24   time:         0.10

Timestep itime:     25   time:         0.10

Total mass, CN mass:   0.42000000E-20  0.43802594E-11

History output #     26 at itime:     25   time:         0.10

Timestep itime:     26   time:         0.10

Total mass, CN mass:   0.42000000E-20  0.43802593E-11

History output #     27 at itime:     26   time:         0.10

Timestep itime:     27   time:         0.11

Total mass, CN mass:   0.42000000E-20  0.43802593E-11

History output #     28 at itime:     27   time:         0.11

Timestep itime:     28   time:         0.11

Total mass, CN mass:   0.42000000E-20  0.43802592E-11

History output #     29 at itime:     28   time:         0.11

Timestep itime:     29   time:         0.12

Total mass, CN mass:   0.42000000E-20  0.43802591E-11

History output #     30 at itime:     29   time:         0.12

Timestep itime:     30   time:         0.12

Total mass, CN mass:   0.42000000E-20  0.43802591E-11

History output #     31 at itime:     30   time:         0.12

Timestep itime:     31   time:         0.12

Total mass, CN mass:   0.42000000E-20  0.43802590E-11

History output #     32 at itime:     31   time:         0.12

Timestep itime:     32   time:         0.13

Total mass, CN mass:   0.42000000E-20  0.43802590E-11

History output #     33 at itime:     32   time:         0.13

Timestep itime:     33   time:         0.13

Total mass, CN mass:   0.42000000E-20  0.43802589E-11

History output #     34 at itime:     33   time:         0.13

Timestep itime:     34   time:         0.14

Total mass, CN mass:   0.42000000E-20  0.43802588E-11

History output #     35 at itime:     34   time:         0.14

Timestep itime:     35   time:         0.14

Total mass, CN mass:   0.42000000E-20  0.43802588E-11

History output #     36 at itime:     35   time:         0.14

Timestep itime:     36   time:         0.14

Total mass, CN mass:   0.42000000E-20  0.43802587E-11

History output #     37 at itime:     36   time:         0.14

Timestep itime:     37   time:         0.15

Total mass, CN mass:   0.42000000E-20  0.43802587E-11

History output #     38 at itime:     37   time:         0.15

Timestep itime:     38   time:         0.15

Total mass, CN mass:   0.42000000E-20  0.43802586E-11

History output #     39 at itime:     38   time:         0.15

Timestep itime:     39   time:         0.16

Total mass, CN mass:   0.42000000E-20  0.43802585E-11

History output #     40 at itime:     39   time:         0.16

Timestep itime:     40   time:         0.16

Total mass, CN mass:   0.42000000E-20  0.43802585E-11

History output #     41 at itime:     40   time:         0.16

Timestep itime:     41   time:         0.16

Total mass, CN mass:   0.42000000E-20  0.43802584E-11

History output #     42 at itime:     41   time:         0.16

Timestep itime:     42   time:         0.17

Total mass, CN mass:   0.42000000E-20  0.43802584E-11

History output #     43 at itime:     42   time:         0.17

Timestep itime:     43   time:         0.17

Total mass, CN mass:   0.42000000E-20  0.43802583E-11

History output #     44 at itime:     43   time:         0.17

Timestep itime:     44   time:         0.18

Total mass, CN mass:   0.42000000E-20  0.43802583E-11

History output #     45 at itime:     44   time:         0.18

Timestep itime:     45   time:         0.18

Total mass, CN mass:   0.42000000E-20  0.43802582E-11

History output #     46 at itime:     45   time:         0.18

Timestep itime:     46   time:         0.18

Total mass, CN mass:   0.42000000E-20  0.43802582E-11

History output #     47 at itime:     46   time:         0.18

Timestep itime:     47   time:         0.19

Total mass, CN mass:   0.42000000E-20  0.43802581E-11

History output #     48 at itime:     47   time:         0.19

Timestep itime:     48   time:         0.19

Total mass, CN mass:   0.42000000E-20  0.43802581E-11

History output #     49 at itime:     48   time:         0.19

Timestep itime:     49   time:         0.19

Total mass, CN mass:   0.42000000E-20  0.43802580E-11

History output #     50 at itime:     49   time:         0.19

Timestep itime:     50   time:         0.20

Total mass, CN mass:   0.42000000E-20  0.43802580E-11

History output #     51 at itime:     50   time:         0.20

Timestep itime:     51   time:         0.20

Total mass, CN mass:   0.42000000E-20  0.43802579E-11

History output #     52 at itime:     51   time:         0.20

Timestep itime:     52   time:         0.21

Total mass, CN mass:   0.42000000E-20  0.43802579E-11

History output #     53 at itime:     52   time:         0.21

Timestep itime:     53   time:         0.21

Total mass, CN mass:   0.42000000E-20  0.43802578E-11

History output #     54 at itime:     53   time:         0.21

Timestep itime:     54   time:         0.21

Total mass, CN mass:   0.42000000E-20  0.43802578E-11

History output #     55 at itime:     54   time:         0.21

Timestep itime:     55   time:         0.22

Total mass, CN mass:   0.42000000E-20  0.43802577E-11

History output #     56 at itime:     55   time:         0.22

Timestep itime:     56   time:         0.22

Total mass, CN mass:   0.42000000E-20  0.43802577E-11

History output #     57 at itime:     56   time:         0.22

Timestep itime:     57   time:         0.23

Total mass, CN mass:   0.42000000E-20  0.43802576E-11

History output #     58 at itime:     57   time:         0.23

Timestep itime:     58   time:         0.23

Total mass, CN mass:   0.42000000E-20  0.43802576E-11

History output #     59 at itime:     58   time:         0.23

Timestep itime:     59   time:         0.23

Total mass, CN mass:   0.42000000E-20  0.43802576E-11

History output #     60 at itime:     59   time:         0.23

Timestep itime:     60   time:         0.24

Total mass, CN mass:   0.42000000E-20  0.43802575E-11

History output #     61 at itime:     60   time:         0.24

Timestep itime:     61   time:         0.24

Total mass, CN mass:   0.42000000E-20  0.43802575E-11

History output #     62 at itime:     61   time:         0.24

Timestep itime:     62   time:         0.25

Total mass, CN mass:   0.42000000E-20  0.43802574E-11

History output #     63 at itime:     62   time:         0.25

Timestep itime:     63   time:         0.25

Total mass, CN mass:   0.42000000E-20  0.43802574E-11

History output #     64 at itime:     63   time:         0.25

Timestep itime:     64   time:         0.25

Total mass, CN mass:   0.42000000E-20  0.43802573E-11

History output #     65 at itime:     64   time:         0.25

Timestep itime:     65   time:         0.26

Total mass, CN mass:   0.42000000E-20  0.43802573E-11

History output #     66 at itime:     65   time:         0.26

Timestep itime:     66   time:         0.26

Total mass, CN mass:   0.42000000E-20  0.43802573E-11

History output #     67 at itime:     66   time:         0.26

Timestep itime:     67   time:         0.27

Total mass, CN mass:   0.42000000E-20  0.43802572E-11

History output #     68 at itime:     67   time:         0.27

Timestep itime:     68   time:         0.27

Total mass, CN mass:   0.42000000E-20  0.43802572E-11

History output #     69 at itime:     68   time:         0.27

Timestep itime:     69   time:         0.27

Total mass, CN mass:   0.42000000E-20  0.43802571E-11

History output #     70 at itime:     69   time:         0.27

Timestep itime:     70   time:         0.28

Total mass, CN mass:   0.42000000E-20  0.43802571E-11

History output #     71 at itime:     70   time:         0.28

Timestep itime:     71   time:         0.28

Total mass, CN mass:   0.42000000E-20  0.43802571E-11

History output #     72 at itime:     71   time:         0.28

Timestep itime:     72   time:         0.29

Total mass, CN mass:   0.42000000E-20  0.43802570E-11

History output #     73 at itime:     72   time:         0.29

Timestep itime:     73   time:         0.29

Total mass, CN mass:   0.42000000E-20  0.43802570E-11

History output #     74 at itime:     73   time:         0.29

Timestep itime:     74   time:         0.29

Total mass, CN mass:   0.42000000E-20  0.43802569E-11

History output #     75 at itime:     74   time:         0.29

Timestep itime:     75   time:         0.30

Total mass, CN mass:   0.42000000E-20  0.43802569E-11

History output #     76 at itime:     75   time:         0.30

Timestep itime:     76   time:         0.30

Total mass, CN mass:   0.42000000E-20  0.43802569E-11

History output #     77 at itime:     76   time:         0.30

Timestep itime:     77   time:         0.31

Total mass, CN mass:   0.42000000E-20  0.43802568E-11

History output #     78 at itime:     77   time:         0.31

Timestep itime:     78   time:         0.31

Total mass, CN mass:   0.42000000E-20  0.43802568E-11

History output #     79 at itime:     78   time:         0.31

Timestep itime:     79   time:         0.31

Total mass, CN mass:   0.42000000E-20  0.43802568E-11

History output #     80 at itime:     79   time:         0.31

Timestep itime:     80   time:         0.32

Total mass, CN mass:   0.42000000E-20  0.43802567E-11

History output #     81 at itime:     80   time:         0.32

Timestep itime:     81   time:         0.32

Total mass, CN mass:   0.42000000E-20  0.43802567E-11

History output #     82 at itime:     81   time:         0.32

Timestep itime:     82   time:         0.33

Total mass, CN mass:   0.42000000E-20  0.43802566E-11

History output #     83 at itime:     82   time:         0.33

Timestep itime:     83   time:         0.33

Total mass, CN mass:   0.42000000E-20  0.43802566E-11

History output #     84 at itime:     83   time:         0.33

Timestep itime:     84   time:         0.33

Total mass, CN mass:   0.42000000E-20  0.43802566E-11

History output #     85 at itime:     84   time:         0.33

Timestep itime:     85   time:         0.34

Total mass, CN mass:   0.42000000E-20  0.43802565E-11

History output #     86 at itime:     85   time:         0.34

Timestep itime:     86   time:         0.34

Total mass, CN mass:   0.42000000E-20  0.43802565E-11

History output #     87 at itime:     86   time:         0.34

Timestep itime:     87   time:         0.35

Total mass, CN mass:   0.42000000E-20  0.43802565E-11

History output #     88 at itime:     87   time:         0.35

Timestep itime:     88   time:         0.35

Total mass, CN mass:   0.42000000E-20  0.43802564E-11

History output #     89 at itime:     88   time:         0.35

Timestep itime:     89   time:         0.35

Total mass, CN mass:   0.42000000E-20  0.43802564E-11

History output #     90 at itime:     89   time:         0.35

Timestep itime:     90   time:         0.36

Total mass, CN mass:   0.42000000E-20  0.43802564E-11

History output #     91 at itime:     90   time:         0.36

Timestep itime:     91   time:         0.36

Total mass, CN mass:   0.42000000E-20  0.43802563E-11

History output #     92 at itime:     91   time:         0.36

Timestep itime:     92   time:         0.37

Total mass, CN mass:   0.42000000E-20  0.43802563E-11

History output #     93 at itime:     92   time:         0.37

Timestep itime:     93   time:         0.37

Total mass, CN mass:   0.42000000E-20  0.43802563E-11

History output #     94 at itime:     93   time:         0.37

Timestep itime:     94   time:         0.37

Total mass, CN mass:   0.42000000E-20  0.43802562E-11

History output #     95 at itime:     94   time:         0.37

Timestep itime:     95   time:         0.38

Total mass, CN mass:   0.42000000E-20  0.43802562E-11

History output #     96 at itime:     95   time:         0.38

Timestep itime:     96   time:         0.38

Total mass, CN mass:   0.42000000E-20  0.43802562E-11

History output #     97 at itime:     96   time:         0.38

Timestep itime:     97   time:         0.39

Total mass, CN mass:   0.42000000E-20  0.43802561E-11

History output #     98 at itime:     97   time:         0.39

Timestep itime:     98   time:         0.39

Total mass, CN mass:   0.42000000E-20  0.43802561E-11

History output #     99 at itime:     98   time:         0.39

Timestep itime:     99   time:         0.39

Total mass, CN mass:   0.42000000E-20  0.43802561E-11

History output #    100 at itime:     99   time:         0.39

Timestep itime:    100   time:         0.40

Total mass, CN mass:   0.42000000E-20  0.43802560E-11

History output #    101 at itime:    100   time:         0.40

Last timestep index:    100
Last time:         0.40
History output was to file model_his.out
Restart output was to file model_res.out
