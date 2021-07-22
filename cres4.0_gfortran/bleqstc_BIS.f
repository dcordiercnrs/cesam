c
 
c******************************************************************************
c
        blockdata bleqstc
c  initialize data for equation of state
c
c  Modified 22/5/90 to include all levels of ionization of Ar and Fe.
c  Previously only the first 15 levels were included, the remainder
c  being forced to be unionized
c
      implicit double precision (a-h,o-z)
        implicit integer(i-n)
      character name*5
      common/eqscntc/ anh0,anhe0,ihvz,iprrad,ihmin
      common/hvabndc/ ab(10),iab
      common/eqdpco/ frhi,bdcoh,bdcoz,idpco
      common/eqphcs/ c(48),ic
      common/potetcc/ chi(125),am(10),iz(10)
      common/hvname/ name(10)
      common /hvcntl/ icnthv,iwrthv,dptst0,dptst1
      common/hvomeg/ iom(26),iom1(20)
      common/hvomcl/ iomfll
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
c
      data anh0,anhe0,ihvz,iprrad,ihmin /0.5d0,6.d0,1,1,0/
      data ab,iab /0.2254d0,0.0549d0,0.4987d0,0.0335d0,0.00197d0,0.0436d0,
     1  0.00403d0,0.0565d0,0.00180d0,0.0795d0,10/
      data frhi,bdcoh,bdcoz,idpco /1.d0,1.d0,1.d0,0/
c  coefficients for s/r phderc.
        data ic,c/4,2.315472d0,7.128660d0,7.504998d0,2.665350d0,7.837752d0,
     1  23.507934d0,23.311317d0,7.987465d0,9.215560d0,26.834068d0,
     2  25.082745d0,8.020509d0,3.693280d0,10.333176d0,9.168960d0,
     3  2.668248d0,2.315472d0,6.748104d0,6.564912d0,2.132280d0,
     4  7.837752d0,21.439740d0,19.080088d0,5.478100d0,9.215560d0,
     5  23.551504d0,19.015888d0,4.679944d0,3.693280d0,8.859868d0,
     6  6.500712d0,1.334124d0,1.157736d0,3.770676d0,4.015224d0,
     7  1.402284d0,8.283420d0,26.184486d0,28.211372d0,10.310306d0,
     8  14.755480d0,45.031658d0,46.909420d0,16.633242d0,7.386560d0,
     9  22.159680d0,22.438048d0,7.664928d0/
      data iz/6,7,8,10,11,12,13,14,18,26/
      data name/'   C','   N','   O','  Ne','  Na',
     .  '  Mg','  Al','  Si','  Ar','  Fe'/
        data am/12.00d0,14.01d0,16.00d0,20.17d0,22.99d0,24.31d0,26.98d0,
     1  28.08d0,39.94d0,55.84d0/
        data (chi(i),i=1,36)
     1  /11.26d0,       24.38d0,        47.86d0,        64.48d0,
     2  391.99d0,       489.84d0,       14.54d0,        29.60d0,
     3  47.43d0,        77.45d0,        97.86d0,        551.92d0,
     4  666.83d0,       13.61d0,        35.15d0,        54.93d0,
     5  77.39d0,        113.87d0,       138.08d0,       739.11d0,
     6  871.12d0,       21.56d0,        41.07d0,        63.5d0,
     7  97.16d0,        126.4d0,        157.91d0,       207.3d0,
     8  239.d0,         1196.d0,        1360.d0,        5.14d0,
     9  47.29d0,        71.65d0,        98.88d0,        138.60d0/
 
        data (chi(i),i=37,72)
     1  /172.36d0,      208.44d0,       264.15d0,       299.78d0,
     2  1465.d0,        1646.d0,        7.64d0,         15.03d0,
     3  80.12d0,        109.29d0,       141.23d0,       186.86d0,
     4  225.31d0,       265.96d0,       327.90d0,       367.36d0,
     5  1761.2d0,       2085.d0,        5.98d0,         18.82d0,
     6  28.44d0,        119.96d0,       153.77d0,       190.42d0,
     7  241.93d0,       285.13d0,       330.1d0,        398.5d0,
     8  441.9d0,        2085.5d0,       2299.d0,        8.15d0,
     9  16.34d0,        33.46d0,        45.13d0,        166.73d0/
 
        data (chi(i),i=73,108)
     1  /205.11d0,      246.41d0,       303.87d0,       351.83d0,
     2  401.3d0,        476.0d0,        523.2d0,        2436.d0,
     3  2666.d0,        15.75d0,        27.62d0,        40.90d0,
     4  59.79d0,        75.0d0,         91.3d0,         124.d0,
     5  143.46d0,       421.d0,         480.d0,         539.5d0,
     6  21.1d0,         688.5d0,        755.5d0,        854.4d0,
     7  918.d0,         4121.d0,        4426.d0,        7.90d0,
     8  16.18d0,        30.64d0,        56.d0,          79.d0,
     9  105.d0,         133.d0,         151.d0,         235.d0/
 
        data (chi(i),i=109,125)
     1  /262.d0,        290.d0,         321.d0,         355.d0,
     2  390.d0,         457.d0,         489.d0,         1266.d0,
     3  1358.d0,        1456.d0,        1582.d0,        1689.d0,
     4  1799.d0,        1950.d0,        2045.d0,        8828.d0,
     5  9278.d0/
      data icnthv,iwrthv /0,0/
c
c  limiting arguments in exponential in test for full or no
c  ionization in s/r hvionac.
c  dptst0 is used for first level and element to ionize
c  dptst1 is used for remaining levels
c
      data dptst0,dptst1 /85.d0,19.d0/
c
c  quantities for calculating statistical weights in
c  fonction omegac.
c
      data iom/2,1,2,1,6,9,4,9,6,1,2,1,6,9,4,9,6,1,10,21,28,25,6,25,
     .  28,21/
      data iom1/2,1,1,2,10,15,21,28,28,25,7,6,6,7,25,30,28,21,21,10/
c  flag for statistical weight omegac for fully ionized atom.
c  before 5/1/84 omegac was always set to 15, for some bizarre reason.
c  for transition introduce flag iomfll so that iomfll = 0 corresponds
c  to old situation, and iomfll .ne. 0 gives correct value omegac = 1.
      data iomfll /1/
c
c  controls for Coulomb effect.
c
      data epsdmu, icoulm, iclmit, iclmsd, epssdr
     *    / 1.d-12,   0,       1,      1,  1.d-3  /
      end
