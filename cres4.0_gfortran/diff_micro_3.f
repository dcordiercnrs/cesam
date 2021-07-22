 
c*************************************************************************
 
	subroutine diff_micro_3(melange,p,t,r,l,m,ro,drox,kap,dkapx,
     1	gradrad,dgradradx,gradad,dgradadx,
     2	xchim,mstar,d,ddx,v,dvdx)
 
c	calcul des coefficients de diffusion microscopique
c	d'apres Michaux-Proffit, Vienne, Inside the stars IAU 137, p.250
c	on tient compte de Z: Vy=-xVx/(1-x-z)
c	avec derivees
 
c	les coefficients de la matice de diffusion sont non nuls sur la
c	diagonale et sur la premiere colonne. On ne considere que deux vecteurs
c	d(1,i) : coeff de d X1 pour element i
c	d(2,i) : coeff de d Xi pour element i
c	pour les algorithmes on posera d(1,1)=0, au lieu de d(1,1)=d(2,1)
c	pour les derivees on aura
c	ddx(1,i) : derivee de d(1,i) / X1
c	ddx(2,i) : derivee de d(2,i) / X1
c	ddx(3,i) : derivee de d(2,i) / Xi (nul, sauf pour He4)
c	on aura de meme ddx(1,1)=ddx(3,1)=0 au lieu de
c	ddx(1,1)=ddx(2,1)=ddx(3,1)
 
c	v(i) : coefficient de Xi
c	dvdx(i): derivee de v(i)/X1
 
c	les derivees / Z et W ne sont pas prise en compte dans cette routine
 
c	Auteur: P.Morel, OCA
c	Conseils: J. Matias, DASGAL, G. Alecian, Evry
c	adaptation a CESAM3: P. Morel, Departement J.D. Cassini, O.C.A.,
c	Observatoire de Nice
 
c	CESAM, Version 3
 
c	19-08-96 : remplacement ln lambda par Cij (c.a.d. cs)
 
c entrees
c	melange=.true.: on est dans une ZC
c	p, t, r, l, m, ro: donnees au point de calcul
c	xchim: composition chimique, par gramme
c	kap: opacite
c	gradad, gradrad: gradients
c	terminaisons x : derivees/ X1 (ie H)
c	mstar: masse avec perte de masse
 
c sorties
c	d(1,i) : coeff de d X1 pour element i
c	d(2,i) : coeff de d Xi pour element i
c	ddx(1,i) : derivee de d(1,i) / X1
c	ddx(2,i) : derivee de d(2,i) / X1
c	ddx(3,i) : derivee de d(2,i) / Xi (nul, sauf pour He4)
c	v(i) : coeff de Xi
c	dvdx(i) : derivee de v / X1
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer i,j,k
 
	real*8	p,t,xchim(1),m,l,r,ro,kap,cs(pnelem,2),as(pnelem,2),
     1	num,den,q,cte1,mstar,thetae,b,c,z,
     2	lnlamb(pnelem,2),grav,gradrad,gradad,d_turb,ai(pnelem),
     3	zi(pnelem),sa(pnelem,2),vth,unmz,drox,dkapx,d_conv,
     4	dgradradx,dgradadx,d(2,pnelem),ddx(3,pnelem),v(1),dvdx(1),
     5	ds(2,pnelem),vs(pnelem),stor,unpdx,dstor,z_diff,w_diff
 
	logical melange
	
	logical init
	data init/.true./
	
	save init,ai,zi,as,sa,d_turb,thetae,cte1,b,d_conv,z_diff,w_diff
 
2000	format(1x,1p8d10.3)
2001	format(1x,1p10d8.1)
 
	if(init)then
	 init=.false.
	 unpdx=1.000001d0
	 write(6,*)
	 write(2,*)
	 d_turb=0.
	 d_conv=1.d13
	 write(6,10)d_conv,d_turb
	 write(2,10)d_conv,d_turb
10	 format(1x,'Diff. micro. de Proffit et Michaux, dans ZR. ',/,
     1	1x,'H et He4 sont diffuses, les autres sont element test.',/,
     2	1x,'Dans ZC, Dconv=',1pd10.3,' Coeff. diff. turb.=',1pd10.3)
	 write(2,*)'pour He4: X Vx = - Y Vy ==> Vy = - X / (1 - X - Z) Vx'
	 write(6,*)'pour He4: X Vx = - Y Vy ==> Vy = - X / (1 - X - Z) Vx'	
	 if(iw .gt. 0)then
	  w_diff=0.
c	  print*,'entrer w_diff'
c	  read*,w_diff
	  write(6,11)w_diff
	  write(2,11)w_diff	
11	  format(t2,' avec diffusion du moment angulaire w_diff =',1pd10.3,//)
	 endif	
	 if(iz .gt. 0)then	
	  z_diff=0.
c	  print*,'entrer z_diff'
c	  read*,w_diff
	  write(6,12)z_diff
	  write(2,12)z_diff	
12	  format(t2,' avec diffusion des elements lourds z_diff =',1pd10.3,//)
	 endif	
	 write(6,*)
	 write(2,*)
	
c	 Case mixture of H et He & thetae=1.   :   zi=1  zj=2
 
	 thetae = 1.
	 b = 15./16. * sqrt(2.*amu/5./pi) * sqrt(kbol)**5 / echarg**4
	
	 do i=1,nchim
	  if(nom_elem(i) .eq. ' H ')then
	   zi(i)=1.	!charge
	   ai(i)=ah	!masse atm.
	  elseif(nom_elem(i) .eq. ' H2')then
	   zi(i)=1.
	   ai(i)=ah2	
	  elseif(nom_elem(i) .eq. 'He3')then
	   zi(i)=2.
	   ai(i)=ahe3	  	  	
	  elseif(nom_elem(i) .eq. 'He4')then	
c	   ihe4=i		!deja connu
	   zi(i)=2.
	   ai(i)=ahe4	
	  elseif(nom_elem(i) .eq. 'Li7')then
	   zi(i)=3.
	   ai(i)=ali7	
	  elseif(nom_elem(i) .eq. 'Be7')then
	   zi(i)=4.
	   ai(i)=abe7	  	  	
	  elseif(nom_elem(i) .eq. 'C12')then
	   zi(i)=6.
	   ai(i)=ac12	  	
	  elseif(nom_elem(i) .eq. 'C13')then
	   zi(i)=6.
	   ai(i)=ac13	  	
	  elseif(nom_elem(i) .eq. 'N14')then
	   zi(i)=7.
	   ai(i)=an14	  	
	  elseif(nom_elem(i) .eq. 'N15')then
	   zi(i)=7.
	   ai(i)=an15	  	
	  elseif(nom_elem(i) .eq. 'O16')then
	   zi(i)=8.
	   ai(i)=ao16	  		
	  elseif(nom_elem(i) .eq. 'O17')then
	   zi(i)=8.
	   ai(i)=ao17	  	
	  endif
	 enddo
	
c	 masses reduites
 
	 do i=1,nchim
	  as(i,1)=ai(i)*ai(1)/(ai(i)+ai(1))
	  sa(i,1)=sqrt(as(i,1))
	  as(i,2)=ai(i)*ai(ihe4)/(ai(i)+ai(ihe4))
	  sa(i,2)=sqrt(as(i,2))
	 enddo
	
	 cte1=g*msol/rsol**2
	endif		!initialisation
	
c	le Z
	
	if(iz .gt. 1)then
	 z=xchim(iz)
	else
	 z=z0
	endif
	unmz=1.-z
	
	do i=1,nbelem	!mises a 0
	 v(i)=0.
	 vs(i)=0.
	 dvdx(i)=0.	
	 do j=1,2
	  d(j,i)=0.
	  ds(j,i)=0.	  	
	  ddx(j,i)=0.
	 enddo	
	 ddx(3,i)=0.
	enddo
 
	if(melange) then
	 do i=1,nbelem
	  d(2,i)=d_conv
	 enddo
	else
	
c	 les equations et les derivees numeriques en Xi,
 
	 do k=0,nchim+1		!0 pour la fonction
	  if(k .gt. 0 .and. k .le. nchim)then	!pour les elements chimiques
	   stor=xchim(k)
	   xchim(k)=stor*unpdx
	   dstor=xchim(k)-stor
	  elseif(k .eq. nchim+1)then		!derivee / ro
	   stor=ro
	   ro=stor*unpdx
	   dstor=ro-stor
	  endif	
 
c	  logarithme de Coulomb et integrale cij pour H=X et He4=Y
 
	  do i=1,nchim
	   call coulombln(zi(i),zi(1)   ,thetae,ro,xchim(1),
     1	t,lnlamb(i,1),cs(i,1))
	   call coulombln(zi(i),zi(ihe4),thetae,ro,xchim(1),
     1	t,lnlamb(i,2),cs(i,2))
	  enddo
	
	  c = b * sqrt(t)**5 / ro / cs(1,2) / (.7+.3*xchim(1))
	  q= 2./sqrt(5.)/ro*b*sqrt(t)**5
	  if(r .gt. 0.)then
	   grav=cte1*ro*m/p/r**2
	  else
	   grav=0.
	  endif
	  vth=0.54*b/ro*(4.75*xchim(1)+2.25)*sqrt(t)**5/(5.+cs(1,2))*
     1	gradrad*grav		!formule 19
c	  write(6,2000)c,q,grav,vth,m,r
 
	  do i=1,nchim	
	
	   if(i .eq. 1)then	!pour l'hydrogene
	    vs(1)= c * (1.25 + 1.125*gradrad)*grav*(1.-xchim(1))
	    ds(2,1)= c*(3.+xchim(1))/(1.+xchim(1))/(3.+5.*xchim(1))
	    ds(1,1)=ds(2,1)
	    	
	   elseif(i .eq. ihe4)then	!pour He4: XVx=-Y Vy ==> Vy=-X/(1-X-Z)Vx
	    ds(2,ihe4)=ds(2,1)	!coeff de d Y
	    vs(ihe4)=-vs(1)*xchim(1)/(unmz-xchim(1))
	
	   else		!pour les autres elements: elements test
	    num=sa(i,1)*cs(i,1)-sa(i,2)*cs(i,2)
	    den=xchim(1)*num+sa(i,2)*cs(i,2)
	    ds(2,i)=q/zi(i)**2/den
	
	    ds(1,i)=(ds(2,i)*
     1	(1./(1.+xchim(1))-10./(5.*xchim(1)+3.))+
     2	xchim(1)*(num/den-.23)*ds(2,1))*ah/ai(i)*xchim(i)  !sur dX1/dm
	    vs(i)=ds(2,i)*(1.+zi(i)-ai(i)*(5.*xchim(1)+3.)/4.)*grav+
     1	xchim(1)*(num/den-.23)*vs(1)*ah-vth	
	   endif
	  enddo
	
	  if(k .eq. 0)then		!les fonctions
	   do i=1,nchim	
	    v(i)=vs(i)
	    do j=1,2
	     d(j,i)=ds(j,i)
	    enddo
	   enddo
	   if(iw .gt. 0)d(2,iw)=w_diff	!moment angulaire
	   if(iz .gt. 0)d(2,iz)=z_diff	!elements lourds	
	
c	   print*,'dans diff_micro, v/d1/d2'
c	   write(6,2000)(v(i),i=1,nchim)	
c	   write(6,2000)(d(1,i),i=1,nchim)
c	   write(6,2000)(d(2,i),i=1,nchim)
 
c	   d(1,i) : coeff de d X1 pour element i
c	   d(2,i) : coeff de d Xi pour element i
c	   ddx(1,i) : derivee de d(1,i) / X1
c	   ddx(2,i) : derivee de d(2,i) / X1
c	   ddx(3,i) : derivee de d(2,i) / Xi (nul, sauf pour He4)
c	   v(i) : coeff de Xi
c	   dvdx(i) : derivee de v / X1
	
	  else					!derivees numerique / Xk
	   if(k .eq. 1)then
	    do i=1,nchim	
	     ddx(1,i)=(d(1,i)-ds(1,i))/dstor  !derivee de d(1,i) / X1		
	     ddx(2,i)=(d(2,i)-ds(2,i))/dstor  !derivee de d(2,i) / X1	
	     dvdx(i)=(v(i)-vs(i))/dstor  !derivee de vi(i) / X1
	    enddo	     	
	    xchim(1)=stor
	   elseif(k .le. nchim)then	
	    ddx(3,k)=(d(2,k)-ds(2,k))/dstor  !derivee de d(1,i) / Xi
	    xchim(k)=stor
	   elseif(k .eq. nchim+1)then
	    do i=1,nchim
	     ddx(1,i)=ddx(1,i)+(d(1,i)-ds(1,i))/dstor*drox	!derivee /ro
	     ddx(2,i)=ddx(2,i)+(d(2,i)-ds(2,i))/dstor*drox	!derivee /ro
	     dvdx(i)=dvdx(i)+(v(i)-vs(i))/dstor*drox	!derivee /ro
	    enddo	      	
	    ro=stor
	   endif
c	   print*,k,' ddx1i,ddx2i,ddx3i,dvdx1i'
c	   do i=1,nchim
c	    write(6,2000)ddx(1,i),ddx(2,i),ddx(3,i),dvdx(i)
c	   enddo
c	   pause	   	
	  endif
	 enddo
	
c	 diffusion turbulente
	
	 do i=1,nchim
	  d(2,i)=d(2,i)+d_turb
	 enddo
	
	 d(1,1)=0.
	 ddx(1,1)=0.
	 ddx(3,1)=0.
	
	endif		!melange
 
	return
	
	end
 
