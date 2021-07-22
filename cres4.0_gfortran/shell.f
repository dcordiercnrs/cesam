 
c************************************************************************
 
	subroutine shell(n,arr)
	
c	tri d'un tableau 50 < n < 100
 
c entree
c	n : nombre de points
 
c entree/sortie
c	arr : tableau
 
c	d'apres Numerical receicipes, p 229	
c	P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
 
	implicit none
	
	integer n,lnb2,m,nn,k,i,j,l
	
	real*8 arr(1),t
	
	lnb2=log(float(n))/.69314718+1.e-5
	m=n
	do nn=1,lnb2
	 m=m/2
	 k=n-m
	 do j=1,k
	  i=j
3	  l=i+m
	  if(arr(l) .lt. arr(i))then
	   t=arr(i)
	   arr(i)=arr(l)
	   arr(l)=t
	   i=i-m
	   goto3
	  endif
	 enddo
	enddo
	
	return
	
	end  	
