c
c******************************************************************************
c
	subroutine scopy(n,a,na,b,nb)
      implicit double precision(a-h,o-z)
	implicit integer(i-n)
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=a(ia)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end
