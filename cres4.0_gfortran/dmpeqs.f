c
c*******************************************************************************
c
	subroutine dmpeqs
c  dumps commons from s/r eqstfc
      implicit double precision (a-h, o-z)
	implicit integer(i-n)
      common/eqstdc/ c1(94)
      common/eqsoutc/ c2(230)
	common/dmudec/ c3(10),idmu
 
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr
c
      write(istdpr,105)
      write(istdpr,110) c1
      write(istdpr,115) c2
      write(istdpr,120) c3
      return
  105 format(///' output from dmpeqs:'/1x,20(1h*))
  110 format(/' common/eqstdc/:'/1p4e13.5/(10e13.5))
  115 format(/' common/eqsoutc/:'/(1p10e13.5))
  120 format(/' common/dmudec/:'/1p10e13.5)
      end
