 
c********************************************************************
 
	subroutine bvald(x,y,m,l,r,d)
 
c	d(m*(j-1)+r):=d(r,j):= derivee d'ordre r-1,  1 .le. r .le. m, de
c	la l-m+j-ieme B-splines, 1 .le. j .le. m, d'ordre m non nulle
c	au point x, y(l) .le.  x .lt.  y(l+1).
c	les r derivees de 0 a r-1 sont calculees, et on a d(m,m)
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
	implicit none
 
	integer*4 m,l,r
	integer*4 j,i,mj,k,kml
 
	real*8 x,y(1),d(1)
	real*8 cd(10,0:10),q(11),n(10,10),a1,a2,denom
 
c	tests de dimension
 
	if(m .ge. 11)then
	 write(6,*)'dans bvald, modifier les dimensions de cd(',
     1	m,'0:',m,'), de q(',m+1,') et de n(',m,',',m,')'
	 stop
	endif
 
	if(r .gt. m)then
	 write(6,*)'dans bvald, r=',r,'>m=',m
	 stop
	endif
 
c	calcul des B-splines d'ordre 1 a m non nulles sur [y(l),y(l+1)[
c	adapte de schumaker 5-5
 
	do i=1,10
 	 do j=1,10
	  n(i,j)=0.
	 enddo
	enddo
 
	do i=1,m+1
	 q(i)=0.
	enddo
	q(m)=1./(y(l+1)-y(l))
	do i=1,m
	 n(1,i)=q(i)*(y(l-m+i+1)-y(l-m+i))
	enddo
 
	do j=2,m
	 do i=m-j+1,m
	  denom=y(i+l-m+j)-y(i+l-m)
	  a1=(x-y(i+l-m))/denom
	  a2=1.-a1
	  q(i)=a1*q(i)+a2*q(i+1)
	  n(j,i)=q(i)*(y(l-m+i+j)-y(l-m+i))
	 enddo
	enddo
 
c	c(j,i) pour chacune des B-splines non nulle sur [y(l),y(l+1)[
c	adapte de schumaker 5-10
 
	do k=l-m+1,l
	 kml=m*(k-l+m-1)
	 do i=1,10
	  do j=0,10
	   cd(i,j)=0.
	  enddo	!j
	 enddo	!i
	 cd(1,1)=1.
 
	 do j=2,r
	  mj=m-j+1
	  do i=1,j
	   denom=y(k+i-1+mj)-y(k+i-1)
	   if(denom.eq.0.)then
	    cd(j,i)=0.
	 	 	  else
	    cd(j,i)=float(mj)/denom*(cd(j-1,i)-cd(j-1,i-1))
	   endif
	  enddo	!i
	 enddo		!j
 
c	do i=1,r
c	 write(6,1000)(cd(i,j),j=1,r)
1000	 format(1x,1p6d10.3)
c	enddo
c	write(6,*)' '
 
	 do j=1,r
	  d(kml+j)=0.
	  do i=1,j
	   d(kml+j)=d(kml+j)+cd(j,i)*n(m-j+1,k-l+m+i-1)
	  enddo	!i
	 enddo		!j
	enddo		!k
 
	return
 
	end
