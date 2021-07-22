 
c****************************************************************
 
	subroutine kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
 
c	ajout de l'opacite conductive
 
c	adaptation de la partie concernee du code de Geneve procuree par
c	Y. Lebreton et selon ses instructions
 
c	expression of Iben(1975  AP.j. 196,525 Appendix A) and numerical
c	evaluations of the derivatives with respect to density and
c	temperature assuming dlro=0.001 and dlt=0.0001. Derivatives
c	are calculated assuming vmye=2./(1.+x)=constante=nb. nucleons par
c	electron libre
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version: 15 02 92
 
c entrees:
c	xh(1)=X : comp. chim. en fraction de masse
c	t : temperature K
c	ro : densite cgs
 
c entree/sortie :
c	kapa : opacite gr / cm2)
c	dkapdt : kappa / d t
c	dkapdr : kappa / d densite
c       dkapdx : kappa / d xchim(1)
 
	implicit none
 
	include 'cesam_3.parametres'
 
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
 
	integer j,j12,j14,j16,j20,j22,j24,i,k
c	data j12,j14,j16,j20,j22,j24/6*0/
 
	real*8 xh(1),t,ro,kappa,dkapdt,dkapdr,dkapdx,t1,ro1
	real*8 vmye,vmye_log,za,za_log,zb,xz2,xhh,z,
     1	rol6,tl6,dellg,del,eta0,eta2,
     2	a1,a2,b1,b2,rnedne,flg,blamr2,alfa,zc,f,thxlg,thx,penktl,
     3	ef,efm,gam,glg,vkcc,thy,thc,thylg,thclg,
     4	vcond,vkchr,kappa_c,dkapdr_c,dkapdt_c,dkapdx_c,
     5	xx(8),zi(8),a(8),unpdx,store,store0,dstore
        logical init
        data init/.true./
	data unpdx/1.00001d0/
 
	data j12,j14,j16,j20,j22,j24/6*0/
	data zi/1.d0,2.d0,6.d0,7.d0,8.d0,10.d0,10.d0,12.d0/
	data a/1.d0,4.d0,12.d0,14.d0,16.d0,20.d0,22.d0,24.d0/
 
c	logical init /.true./
 
2000	format((1x,1p8d10.3))
 
c	write(6,2000)(xh(i),i=1,nchim),t,ro
c	write(6,*)(nom_elem(i),i=1,nchim)
 
	if(init)then			!indices des elements
	 if(nom_elem(1) .ne. ' H ')return
	 init=.false.
	 do j=3,nchim
	  if(nom_elem(j) .eq. 'C12')then
	   j12=j
	  elseif(nom_elem(j) .eq. 'N14')then
	   j14=j
	  elseif(nom_elem(j) .eq. 'O16')then
	   j16=j
	  elseif(nom_elem(j) .eq. 'N20')then
	   j20=j
	  elseif(nom_elem(j) .eq. 'N22')then
	   j22=j
	  elseif(nom_elem(j) .eq. 'M24')then
	   j24=j
	  endif
	 enddo
	 write(6,*)'opacites conductives Iben Apj 196,515,1975, avec Y=1-X-Z'
	 write(6,*)' '
	 write(2,*)'opacites conductives Iben Apj 196,515,1975, avec Y=1-X-Z'
	 write(2,*)' '
	endif
 
c	test whether electron conduction is likely to be important
c	calculation of conductive opacity according to the analytic
 
	if(ro .le. 1.d-4 .or. t .le. 2.d4)return
	
	if(iz .gt. 1)then	!Z peut etre diffuse
	 z=xh(iz)
	else
	 z=z0
	endif
 
c	expression of Iben(1975  AP.j. 196,525 Appendix A) and numerical
c	evaluations of the derivatives with respect to density and
c	temperature and chemical composition X.
 
	ro1=ro
	t1=t
	do k=1,4
	 if(k .eq. 2)then
	  store0=t1
	  store=store0*unpdx
	  dstore=store-store0
	  t1=store
	 elseif(k .eq. 3)then
	  store0=ro1
	  store=store0*unpdx
	  dstore=store-store0
	  ro1=store
	 elseif(k .eq. 4)then
	  store0=xh(1)
	  store=max(1.d-5,store0*unpdx)
	  dstore=store-store0
	  xhh=store
	 endif
 
c 	 valeurs de T6 et ro6 en Log10
 
	 tl6=log10(t1)-6.
	 rol6=log10(ro1)-6.
 
c	 abondances en fraction de masse de H, He4, C12,
c	 N14, O16, Ne20, Ne22 et Mg24
 
	 do i=1,8
	  xx(i)=0.
	 enddo
	 xx(1)=xhh
	 xx(2)=1.-xhh-z
	 if(j12 .ge. 3)xx(3)=xh(j12)
	 if(j14 .ge. 3)xx(4)=xh(j14)
	 if(j16 .ge. 3)xx(5)=xh(j16)
	 if(j20 .ge. 3)xx(6)=xh(j20)
	 if(j22 .ge. 3)xx(7)=xh(j22)
	 if(j24 .ge. 3)xx(8)=xh(j24)
 
c	 l'opacite conductive
 
	 vmye=2./(1.+xx(1))
	 vmye_log=log10(vmye)
	 za=0.
	 zb=0.
	 do i=1,8
          xz2=xx(i)*zi(i)**2
          za=za+xz2/(a(i)**(1./3.))
	  zb=zb+xz2/a(i)
	 enddo
	 za_log=log10(za)
 
	 if(rol6 .le. 0.3) then		!Hubbart pour ro<2.d6
	  dellg=rol6+6.-1.5*tl6-vmye_log		!A1
	  del=10.**(dellg)
	  eta0=10.**(-0.52255+2.*dellg/3.)		!A6
	  eta2=eta0**2
 
c	  logique de A2 a A5
 
	  a1=-3.29243+log10(del*(1.+0.02804*del))	!A3
	  b1=-4.80946+log10(del*del*(1.+9.376/eta2))	!A4
	  if(dellg .le. 0.645) then
	   flg=-3.2862+log10(del*(1.+0.024417*del))	!A2	
	  elseif(dellg .le. 2.)then
	   flg=a1					!A3
	  elseif(dellg .le. 2.5) then
	   flg=2.*a1*(2.5-dellg)+2.*b1*(dellg-2.)	!A5
	  else
	   flg=b1					!A4
	  endif
 
c	  logique de A8 a A10
 
	  a2=log10(1.+0.021876*del)			!A8
	  b2=log10(0.4*eta0*(1.+4.1124/eta2))		!A9
	  if(dellg . le. 1.5)then	
	   penktl=a2
	  elseif(dellg .le. 2.) then
	   penktl=2.*a2*(2.-dellg)+2.*b2*(dellg-1.5)	!A10
	  else
	   flg=b1					!A9
          endif
 
c	  logique de A11, A12
 
	  if(del .le. 40.)then
	   rnedne=1.0-0.01*del*(2.8966-0.034838*del)	!A11
	  else
	   rnedne=(1.5/eta0)*(1.-0.8225/eta2)		!A12
	  endif
 
c	  logique de A13 a A21
 
	  blamr2=9.24735E-3*10.**(dellg-0.5*tl6-penktl)	!A7
     1	*(vmye*zb+rnedne)
	  alfa=0.5*log10(blamr2)
	  if(alfa .le. -3.)then		!pour H
	   thxlg=1.048-0.124*alfa				!A13
	  elseif(alfa .le. -1.)then
	   thxlg=0.13-alfa*(0.745+0.105*alfa)			!A14
	  else
	   thxlg=0.185-0.558*alfa				!A15
	  endif
	  if(alfa .le. -3.)then		!pour He4
	   thylg=0.937-0.111*alfa				!A16
	   elseif(alfa .le. 0.)then
	    thylg=0.24-alfa*(0.55+0.0689*alfa)		!A17
	   else
	    thylg=0.24-0.6*alfa				!A18
	  endif
	  if(alfa .le. -2.5)then		!pour C
	   thclg=1.27-0.1*alfa					!A19
	  elseif(alfa .le. .5)then
	   thclg=0.727-alfa*(0.511+0.0778*alfa)		!A20
	  else
	   thclg=0.843-0.785*alfa				!A21
	  endif
	  thx=10.**(thxlg)
	  thy=10.**(thylg)
	  thc=10.**(thclg)
 
c	  zc is an effective abundance of C,N,O,Ne weighted by the
c	  squares of the charges and normalised with respect to C since
c	  electron conduction has been calculated by Hubbard+Lampe only
c	  for pure H,pure He and pure C
 
	  zc=(3.*xx(3)+3.5*xx(4)+4.*xx(5)+5.*xx(6)+4.54*xx(7)+6.*xx(8))/3. !A23
	  vkchr=(xx(1)*thx+xx(2)*thy+zc*thc)*10.**(-tl6-flg)	   !A22
 
	  if(rol6 .le. 0.)then	 !Hubbart en dessous de ro=1.d6
	   vcond=vkchr
	  else	!Hubbart + Canuto pour 1.d6< ro < 2.d6
	   ef=sqrt(1.+10.**(2.*(rol6-vmye_log)/3.)) - 1.	!A24
	   gam=22.76*10.**(rol6/3.-tl6)*za		!A25
	   efm=min(1.d0,0.5d0+log10(ef))				!A30
	   glg=(0.873-0.298*za_log+(0.333-0.168*za_log)*efm)*		!A29
     1	(1.-(1.+gam)**(-0.85))
	   vkcc=6.753E-8*10.**(2.*tl6-glg)*zb/(ef**2.5*(1.+ef))	!A28
	   f=0.5*(1.-cos(pi*rol6/0.3))				!A32
	   vcond=10.**((1.-f)*log10(vkchr)+f*log10(vkcc))		!A31
	  endif
	 else	!Canuto pour ro >2.d6
	  ef=sqrt(1.+10.**(2.*(rol6-vmye_log)/3.)) - 1.	!A24
	  gam=22.76*10.**(rol6/3.-tl6)*za		!A25
	  efm=min(1.d0,0.5d0+log10(ef))				!A30
	  glg=(0.873-0.298*za_log+(0.333-0.168*za_log)*efm)*		!A29
     1	(1.-(1.+gam)**(-0.85))
	  vkcc=6.753E-8*10.**(2.*tl6-glg)*zb/(ef**2.5*(1.+ef))	!A28
	  vcond=vkcc
	 endif
	 if(k .eq. 1)then
	  kappa_c=vcond
	 elseif(k .eq. 2)then
	  t1=t
	  dkapdt_c=(vcond-kappa_c)/dstore
	 elseif(k .eq. 3)then
	  ro1=ro
	  dkapdr_c=(vcond-kappa_c)/dstore
	 elseif(k .eq. 4)then
	  xhh=xh(1)
	  dkapdx_c=(vcond-kappa_c)/dstore
	 endif
	enddo		!k
 
c	write(6,2000)kappa,kappa_c,t,ro,xh(1)
 
	kappa=kappa*kappa_c/(kappa+kappa_c)
	dkapdr=(kappa**2*dkapdr_c+kappa_c**2*dkapdr)/(kappa+kappa_c)**2
	dkapdt=(kappa**2*dkapdt_c+kappa_c**2*dkapdt)/(kappa+kappa_c)**2
	dkapdx=(kappa**2*dkapdx_c+kappa_c**2*dkapdx)/(kappa+kappa_c)**2
 
c	kappa=kappa_c
c	dkapdr=dkapdr_c
c	dkapdt=dkapdt_c
c	dkapdx=dkapdx_c
 
	return
 
	end
