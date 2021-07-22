 
c************************************************************************
 
	function sdot(n,a,i1,b,i2)
c
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 10/3/90.
c
c  Note: this double precision version of the routine has same name
c  as single precision version, but is distinguished by its file name
c
      implicit double precision(a-h,o-z)
	implicit integer(i-n)
      dimension a(1),b(1)
      sdot=0.
      j1=1
      j2=1
      do 1 i=1,n
      sdot=sdot+a(j1)*b(j2)
      j1=j1+i1
      j2=j2+i2
    1 continue
      return
      end
