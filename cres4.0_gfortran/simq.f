 
c*************************************************************
 
	subroutine simq(a,b,n)
 
c	resolution d'un systeme lineaire -ssp 360
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 n
	integer*4 jj,j,jy,it,i,ij,imax,i1,i2,iqs,ix,ixj,jx,ixjx,
     1jjx,ny,ia,ib,ic,k
 
	real*8 a(1),b(1)
	real*8 tol,biga,save
 
c        forward solution
      tol=0.
      jj=-n
      do 65 j=1,n
      jy=j+1
      jj=jj+n+1
      biga=0.
      it=jj-j
      do 30 i=j,n
c        search for maximum coefficient in column
      ij=it+i
      if(abs(biga)-abs(a(ij)))20,30,30
   20 biga=a(ij)
      imax=i
   30 continue
c        test for pivot less than tolerance (singular matrix)
      if(abs(biga)-tol) 35,35,40
 35   print 36,biga,tol
 36   format(2x,'dans simq biga=',1pe10.3,' < tol=',1pe10.3)
      stop
c        interchange rows if necessary
   40 i1=j+n*(j-2)
      it=imax-j
      do 50 k=j,n
      i1=i1+n
      i2=i1+it
      save=a(i1)
      a(i1)=a(i2)
      a(i2)=save
c        divide equation by leading coefficient
   50 a(i1)=a(i1)/biga
      save=b(imax)
      b(imax)=b(j)
      b(j)=save/biga
c        eliminate next variable
      if(j-n) 55,70,55
   55 iqs=n*(j-1)
      do 65 ix=jy,n
      ixj=iqs+ix
      it=j-ix
      do 60 jx=jy,n
      ixjx=n*(jx-1)+ix
      jjx=ixjx+it
   60 a(ixjx)=a(ixjx)-(a(ixj)*a(jjx))
   65 b(ix)=b(ix)-(b(j)*a(ixj))
c        back solution
   70 ny=n-1
      it=n*n
      do 80 j=1,ny
      ia=it-j
      ib=n-j
      ic=n
      do 80 k=1,j
      b(ib)=b(ib)-a(ia)*b(ic)
      ia=ia-n
   80 ic=ic-1
      return
      end
