 
c***************************************************
 
	function neuville(x,xi,fi,n)
 
c	algorithme de Neuville
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 03 06 91
 
c entrees:
c	x: point d'interpolation
c	xi: abscisses
c	fi: fonction
c	n:degre du polynome
 
	implicit none
 
	integer pt
	parameter (pt=20)
 
	integer n,i,j
 
	real*8 x,xi(0:*),fi(0:*),t(0:pt),neuville
 
	if(n .gt. pt)then
	 write(6,*)'dans neuville mettre pt=',n
	 stop
	endif
 
c	write(6,2000)(xi(i),i=0,n)
c	write(6,2000)(fi(i),i=0,n)
2000	format((1x,1p8d10.3))
 
	do i=0,n
	 t(i)=fi(i)
	 do j=i-1,0,-1
c	  write(6,2000)xi(i),xi(j),float(i),float(j)
	  t(j)=t(j+1)+(t(j+1)-t(j))*(x-xi(i))/(xi(i)-xi(j))
	 enddo	!j
	enddo	!i
	neuville=t(0)
 
	return
 
	end
