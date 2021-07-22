 
c*************************************************************************
 
	subroutine diff_cte_3(melange,p,t,r,l,m,ro,drox,kap,dkapx,
     1	gradrad,dgradradx,gradad,dgradadx,
     2	xchim,mstar,d,ddx,v,dvdx)
 
c	calcul du coefficient de diffusion
c	diffusion de Z et du Moment Angulaire
 
c	Ici cas de D=cte
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	CESAM, version 3.2.0.0
 
c entrees
c	melange=.true.: on est dans une ZC
c	p,t,r,l,m,ro: donnees au point de calcul
c	xchim: composition chimique, par gramme
c	kap: opacite
c	gradad, gradrad: gradients
c	mstar: masse avec perte de masse
 
c sorties
c	d(1,i) : coeff de d X1 pour element i
c	d(2,i) : coeff de d Xi pour element i
c	ddx(1,i) : derivee de d(1,i) / X1
c	ddx(2,i) : derivee de d(2,i) / X1
c	ddx(3,i) : derivee de d(1,i) / Xi (nul, sauf pour He4)
c	v(i) : coeff de Xi
c	dvdx(i) : derivee de v / X1
c	v: vitesse de diffusion
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer i,j
 
	real*8	p,t,xchim(1),m,l,r,ro,d_rad,d(2,*),kap,v(1),vi,gradad,z_diff,
     1	gradrad,mstar,w_diff,drox,dkapx,dgradradx,dgradadx,d_conv,
     2	ddx(3,*),dvdx(1)
 
	logical melange
	
	logical init
	data init/.true./
 
	save init,d_rad,vi,w_diff,z_diff,d_conv
 
2000	format(1x,1p8d10.3)
 
	if(init)then
	 init=.false.
	 d_rad=100.
	 vi=0.
c	 write(6,*)'d_rad, vi?'
c	 read(5,*)d_rad,vi
	 write(6,*)
	 write(2,*)
	 d_conv=1.d13
	 write(6,10)d_rad,d_conv,vi
	 write(2,10)d_rad,d_conv,vi
10	 format(1x,'Coeff. de diff. constant, dans ZR :', 1pd10.3,
     1	', dans ZC:',1pd10.3,' vitesse de diff. constante :',1pd10.3)
	 write(6,*)
	 write(2,*)
	 if(iw .gt. 0)then
	  w_diff=100.
	  write(6,11)w_diff
	  write(2,11)w_diff	
11	  format(t2,' avec diffusion du moment angulaire w_diff =',1pd10.3,//)
	 endif	
	 if(iz .gt. 0)then
	  z_diff=100.
	  write(6,12)z_diff
	  write(2,12)z_diff	
12	  format(t2,' avec diffusion des elements lourds z_diff =',1pd10.3,//)
	 endif	
	endif
	
	do i=1,nbelem	!mises a 0
	 v(i)=0.
	 dvdx(i)=0.	
	 do j=1,2
	  d(j,i)=0.
	  ddx(j,i)=0.
	 enddo	
	 ddx(3,i)=0.
	enddo
	
	if(melange)then
	 do i=1,nbelem
	  d(2,1)=d_conv
	 enddo
	else
	 do i=1,nchim
	  d(2,i)=d_rad
	  v(i)=vi	
	 enddo
 
	 if(iw .gt. 0)d(2,iw)=w_diff	!moment angulaire
	 if(iz .gt. 0)d(2,iz)=z_diff	!elements lourds
	
	endif
	
	return
 
	end
 
