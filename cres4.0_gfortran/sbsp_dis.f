 
c******************************************************************
 
	subroutine sbsp_dis(n,x,f,nd,id,fd,nx,m,xt,knot)
 
c	interpolation SBP1DN avec discontinuites
c	formation de la base nodale
c	calcul des coefficients des splines
c	s'exploite avec sbsp1dn
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 15 09 92
 
c entree
c	n : nombre de tables
c	x : abscisses
c	f(n,nx) : fonction a interpoler, MODIFIEE
c	fd: valeurs des fonctions a droite des discontinuites, MODIFIEE
c	nd : nombre de discontinuites
c	id(0:nd+1) : commence en 0 finit en nd+1!!!
c		indices des discontinuites, la premiere en 1, la derniere en nd
c	m : ordre des B-splines
c	nx : nombre de points
 
c sortie
c	xt : points de table
c	knot : nombre de points de table
c	f(n,nx) : fonction a interpoler, MODIFIEE
 
	implicit none
 
	integer	pnx,pm,pn,pd
	parameter (pnx=4000, pm=4, pd=20, pn=20)
 
	integer nx,m,lx,knot,indpc(pnx+pd),i,j,n,nd,id(0:*),
     1	kk,nc,ij
 
	real*8 x(1),xt(1),f(1),ax((pnx+pd)*pm),fd(1),s(pn*(pnx+pd)),
     1	qx(pm*pm),eps,xc(pnx+pd)
	data eps/1.d-14/
 
2000	format((1x,1p8d10.3))
 
c	tests de dimension
 
	if(nx .gt. pnx)then
	 write(6,*)'dans sbsp_dis mettre le parametre pnx a',nx
	 stop
	elseif(n .gt. pn)then
	 write(6,*)'dans sbsp_dis mettre le parametre pn a',n
	 stop
	elseif(nd .gt. pd)then
	 write(6,*)'dans sbsp_dis mettre le parametre pd a',nd
	 stop
	elseif(m .gt. pm)then
	 write(6,*)'dans sbsp_dis mettre le parametre pm a',m
	 stop
	endif
 
c	calcul des B-splines
 
c	write(6,2000)(x(i),i=1,nx)
c	write(6,2000)(fd(n*(i-1)+1),i=1,nd)
 
	call noeu_dis(x,nx,nd,id,m,xt,knot)
c	write(6,*)nx,nd,knot,(id(i),i=0,nd+1)
c	write(6,2000)(xt(i),i=1,knot)
c	print*,'eps',eps
c	pause
 
	nc=0		!nc: nombre de points de donnee nx+nd
	ij=1		!indice de la discontinuite
	i=1		!indice de la couche
	do while(i .le. nx)
	 nc=nc+1
	 xc(nc)=x(i)
	 do kk=1,n			!s : VT
	  s(n*(nc-1)+kk)=f(n*(i-1)+kk)	!hors discontinuite
	 enddo
	 if(i .eq. id(ij))then		!sur discontinuite
	  nc=nc+1
	  xc(nc)=x(id(ij))+eps
	  do kk=1,n			!s : VT
	   s(n*(nc-1)+kk)=fd(n*(ij-1)+kk)	!sur discontinuites
	  enddo
c	  write(6,*)ij,id(ij),fd(n*(ij-1)+1)
	  ij=min(ij+1,nd)
	 endif
	 i=i+1
	enddo		!i
	
	do i=1,nc		!retour dans f
	 do kk=1,n
	  f(n*(i-1)+kk)=s(n*(i-1)+kk)
	 enddo
	enddo
 
c	write(6,*)'les f, nc=',nc
c	do i=1,nc
c	 write(6,2000)xc(i),(f(n*(i-1)+k),k=1,n)
c	 write(6,*)xc(i),(f(n*(i-1)+k),k=1,n)
c	enddo	!i
c	pause
 
	lx=m
	do i=1,nc
c	 write(6,2000)(xc(j),j=1,nc)
	 call slinf(xc(i),xt,knot,lx)
	 call bval(xc(i),xt,m,lx,qx)
c	 write(6,*)'i,lx,nc,knot/xc(i),xt(lx),xt(lx+1)',i,lx,nc,knot
c	 write(6,2000)xc(i),xt(lx),xt(lx+1)
c	 write(6,*)i,lx,xc(i),xt(lx),xt(lx+1)
	 do j=1,m
	  ax(nc*(j-1)+i)=qx(j)
	 enddo	!j
	 indpc(i)=lx-m+1
c	 write(6,*)'indpc/ax',indpc(i),lx,xc(i),(ax(nc*(j-1)+i),j=1,m)
c	 write(6,2000)(ax(nc*(j-1)+i),j=1,m)
	enddo		!i
 
	call gausdn(ax,f,indpc,nc,m,n)
 
	return
 
	end
