c******************************************************************
 
	subroutine intertot_3(atm,m_ou_r,dxchim,q_int,nd,jlim,lim,m_zc,
     1	p,t,xchim,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,rstar,ltot,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,mstar,w,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,alfa,grad_mj,
     4	delta,cp,hp,gradad,gradrad,vaissala,beta,age,n,ns,
     5	nh1,nhe1,nhe2,etat,opa,conv,nuc,cte)
 
c	interpolation en m ou en r puis calcul des quantites thermodynamiques
 
c	l'extension du modele n'etant valable qu'avec une precision suffisante
c	on fait passer une spline naturelle en q avec les valeurs exactes
c	des derivees premieres en 1 et n
 
c	on tient compte des discon. des derivees premieres aux limites ZR/ZC
c	12 09 96, modif dans les derivees 1-ieres de l'atmosphre
 
c	23 09 96 introduction de d He4/ dr et dHe3 / dr dans vaissala
c	09 10 96 introduction du vecteur rotation
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
c	version 3
 
c entree :
c	le nom du fichier sur lequel se trouve le modele
c	atm=.true.: il y a une atmosphere a prendre en compte
c	m_ou_r='masse' ou 'rayon' suivant la variable d'interpolation
c	le point d'interpolation
c	le nom de la routine thermodynamique du modele (en external)
 
c sortie :
c	les variables au point considere
 
c	les xchim sont /mole
 
	implicit none
	
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'atmosphere_3.common'
	include 'ctephy.common'
	
	integer pnf,pjac1
	parameter (pnf=(pne-1)*(pn+pn_atm+4+6*pnzc),
     1	pjac1=pne_atm*(pm_qs*(pn_atm-1)+1))
	
	integer n,knot,nc,knotc,lim,jlim(2*pnzc),n_t,knot_t,nc_t,knotc_t,
     1	lim_t,jlim_t(2*pnzc),n_,knot_,nc_,knotc_,nes,ns,long,
     2	i,lx,j,ndeb,nd,knot_atm,k,i_dis(6*pnzc+2),ntour,knots

	external long

	real*8	age,bp(pbp),q(pn),qt(pqt),chim(pchim),mc(pnch),bp_atm(pjac1),
     1	mct(pchimt),r2(pn),m23(pn),m_zc(2*pnzc),r_zc(2*pnzc),
     2	r_ov(2*pnzc),mstar,bp_t(pbp),q_t(pn),qt_t(pqt),chim_t(pchim),
     3	mc_t(pnch),mct_t(pchimt),r2_t(pn),m23_t(pn),dt,m_zc_t(2*pnzc),
     4	r_zc_t(2*pnzc),r_ov_t(2*pnzc),mstar_t,bp_(pbp),q_(pn),qt_(pqt),
     5	chim_(pchim),mc_(pnch),mct_(pchimt),r2_(pn),m23_(pn),dt_
	
	real*8	qs(pn+pn_atm),qst(pn+6+pn_atm+2*pnzc),rs(pn+pn_atm),
     1	ms(pn+pn_atm),fx(pnelem),dfxdx(pnelem),f(pnf),frn,fmn,
     2	dfdx1(pne),dfdxn(pne),rstar,ltot,aradias3,dx,x,
     3	chim_atm(pnelem),xchim(1),dxchim(1),q_int,q_atm(pn_atm),
     4	q_atmt(pn_atm*pm_qs+2),teff,grav,l_rac,df_g(pne*(6*pnzc+2)),
     5	df_d(pne*(6*pnzc+2)),x_g(6*pnzc+2),x_d(6*pnzc+2),mi,ri,
     6	lamb,deltap,deltat,deltax,dcpp,dcpt,dcpx,d_grad,gradconv
	
	real*8	p,t,m,l,r,ro,drop,drot,drox,u,dup,dut,dux,dfdx(pnelem),
     1	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon(1),depsp,depst,depsx(1),kap,dkapp,dkapt,dkapx,
     4	alfa,grad_mj,delta,cp,hp,gradad,gradrad,vaissala,beta,
     5	nh1,nhe1,nhe2,corr,depsro,hh,be7,b8,n13,o15,f17,
     6	chimg(pnelem),dchimg(pnelem),nuc_m(pnelem),cte3,w,dw
	
	character*5  m_ou_r
	character*70 modele_atm
 
	logical init,atm,lconv(2*pnzc),lconv_t(2*pnzc),tot_conv,tot_rad,no_lim
c	data init/.false./
 
c	save init
 
	external etat,opa,conv,nuc,cte
 
	data init/.false./
        save init
 
2000	format((1x,1p8d10.3))
 
	if(.not.init)then	!initialisation des constantes
	 init=.true.
	 call cte
	
c	 appel d'initialisation de la routine de reactions nucleaires		
	
	 call nuc(t,ro,xchim,dxchim,qs,.false.,1,
     1	epsilon,depst,depsro,depsx,hh,be7,b8,n13,o15,f17)
	
	 write(6,*)'intertot_3, lecture du modele ',
     1	nom_fich2(:long(nom_fich2))
	 open(unit=4,form='unformatted',status='old',
     1	file=nom_fich2(:long(nom_fich2))//'_B.dat')
	 read(4)age,nbelem,nchim,mtot,alpha,modele,x0,y0,ne,nreac,m_qs,m_ch,
     1	nom_elem,rot_solid,z_cte,w_rot,diffusion,iz,iw,
     2	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     3	r2,m23,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,
     4	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     5	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     6	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     7	r2_,m23_,dt_
	 close (unit=4)
	
	 z0=1.d0-x0-y0
	 mpr=m_qs+1
	 if(iz .gt. 1)nucleo(iz)=1.d0
	 if(iw .gt. 1)nucleo(iw)=1.d0
	 	
	 print*,'Reprise de la solution et exploitation de la superconvergence'
	 print*,'par spline d''ordre 4 avec les der. du modele aux extremites'
	 if(m_qs .lt. 2)then
	  print*,'ordre des splines: ',m_qs,' insufisant'
	  print*,'extension du modele peu fiable'
	  pause'on continue?'
	 endif
 
	 aradias3=aradia/3.
	
	 cte3=4.*pi*rsol**3/msol
	 	
	 nes=ne-1		!on enleve psi
	
	 dx=.1	!pour localiser les derivees premieres
	
c	 si m_ch=2, la base est duale, on passe a l'ordre 3 en rappelant
c	 simplement sbsp1dn	
	
	 if(m_ch .eq. 2)then
	  print*
	  print*,'on passe xchim a l''ordre 3'
	  print*	
	  m_ch=3
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.false.,mc(1),lx,xchim,dxchim)
	 endif
	
c	 extraction des fonctions en tous les q a droite
 
	 do i=1,n
	  x=i
	  qs(i)=x
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x,lx,fx,dfxdx)
	  do j=1,nes
	   f(nes*(i-1)+j)=fx(j)
	  enddo
c	  write(6,2000)(fx(j),j=1,nes)
 
c	  do j=1,lim
c	   if(i .eq. jlim(j))then
c	    print*,'limite ',i,fx(5)
c	   endif
c	  enddo
 
c	  if(.true.)then	!test
	  if(.false.)then
	   p=exp(fx(1))
	   t=exp(fx(2))
	   r=sqrt(abs(fx(3)))
	   l=sqrt(abs(fx(4)))**3
	   m=sqrt(abs(fx(5)))**3
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,fx(5),lx,xchim,dxchim)
	   do j=1,nbelem
	    if(fx(5) .gt. 0.d0)then
	     dxchim(j)=dxchim(j)*2.d0/3.d0/sqrt(fx(5))
	    else
	     dxchim(j)=0.
	    endif
	   enddo
	
	   if(iw .gt. 1 .and. fx(3) .ne. 0.)then		!rotation
	    w=xchim(iw)/fx(3)
	   else
	    w=w_rot
	   endif
	
	   call thermo_3(p,t,xchim,m,l,r,.false.,dxchim,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,0.d0,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)
	   write(6,2000)p,t,r,l,m,ro,xchim(1)
	   pause'apres thermo test 1'
	  endif	!test
	 enddo	!extraction sur i
	
	 frn=fx(3)	!rayon**2 et masse**2/3 sur la derniere couche
	 fmn=fx(5)
	 ns=n
	 x=1			!au centre derivees premieres
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x,lx,fx,dfxdx)
	 do j=1,nes
	  dfdx1(j)=dfxdx(j)
	 enddo
	
	 if(.not.atm)then		!a l'exterieur sans atmosphere
	  x=n
	  call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x,lx,fx,dfxdx)
	  do j=1,nes
	   dfdxn(j)=dfxdx(j)		!derivees premieres
	  enddo
	  rstar=sqrt(fx(3))
	  ltot=sqrt(fx(4))**3
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(nc),lx,chim_atm,dxchim)
	
	 else		!avec atmosphere
	  open(unit=4,form='unformatted',status='old',
     1	file=nom_fich2(:long(nom_fich2))//'_B.atm')
	  read(4)bp_atm,q_atm,q_atmt,teff,grav,ltot,chim_atm,nea,
     1	n_atm,knot_atm,n23,modele_atm
	  close(unit=4)
	  	
c	  do j=1,nbelem
c	   chim_atm(j)=chim_atm(j)/nucleo(j)	!xchim est /gr dans l'atm.
c	  enddo
 
c	  on prend plutot la derniere couche de l'enveloppe
 
	  call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,mc(nc),lx,chim_atm,dxchim)
	
	  rstar=bp_atm(4)		!bp_atm(q,4) est rstar pour tous q
	  if(iz .gt. 1)chim_atm(iz)=chim_atm(iz)/rstar
	  l_rac=ltot**(2.d0/3.d0)
c	  pause'atm debut'
c	  write(6,*)'ns,n_atm,frn,fmn',ns,n_atm,frn,fmn
	  ndeb=1
	  do i=2,n_atm
	   x=i
	   call sbsp1dn(nea,bp_atm,q_atm,q_atmt,n_atm,mpr,knot_atm,.true.,
     1	x,lx,fx,dfxdx)
	   fx(3)=fx(3)**2
	   fx(4)=l_rac
	   fx(5)=fx(5)**(2.d0/3.d0)
c	   print*,i,fx(3),frn,fx(5),fmn
c	   if(fx(3) .gt. frn .and. fx(5) .gt. fmn)then	!elimination des points
	   if(fx(3) .gt. frn)then	!elimination des points en R seulement
	    ns=ns+1		!d'atmosphere dans l'enveloppe
	    qs(ns)=ns
	    do j=1,5		!pour P, T, R, L, M
	     f(nes*(ns-1)+j)=fx(j)
	    enddo
c	    write(6,2000)(fx(j),j=1,nes)
	   else
	    ndeb=i		!le point suivant sera retenu
	   endif		!points de l'atmosphere dans l'enveloppe
	  enddo
	
c	  derivees premieres a l'exterieur
 
	  call sbsp1dn(nea,bp_atm,q_atm,q_atmt,n_atm,mpr,knot_atm,.true.,
     1	dfloat(n_atm),lx,fx,dfxdx)	
	  do j=1,2		!pour P, T
	   dfdxn(j)=dfxdx(j)
	  enddo
	  dfdxn(3)=2.d0*fx(3)*dfxdx(3)
	  dfdxn(5)=2.d0/3.d0*fx(5)**(-1.d0/3.d0)*dfxdx(5)
	  dfdxn(4)=0
 
c	  write(6,2000)(dfdxn(j),j=1,nes)
c	  write(6,*)'ns',ns
c	  pause'atm fin'
	 endif		!si atmosphere
 
c	 on garde les r et les m
 
	 do i=1,ns
	  rs(i)=sqrt(f(nes*(i-1)+3))
	  ms(i)=sqrt(f(nes*(i-1)+5))**3
	 enddo
 
c	 aux limites ZR/ZC: calcul des derivee a gauche et a droite
c	 modele totalement convectif: lim=1,jlim(1)=n,lconv(1)=.false.
c	 modele totalement radiatif: lim=0,jlim(i)=-100,lconv(i)=.false.
 
	 tot_conv=lim .eq. 1 .and. jlim(1) .eq. n .and. .not.lconv(1)
	 tot_rad=lim .eq. 0
	 no_lim=tot_conv .or. tot_rad	
	
c	 no_lim=.true.	!pour test	
c	 print*,'no_lim',no_lim,tot_rad,tot_conv,lim,(jlim(i),i=1,lim)
c	 pause'apres nolim'
	
	 if(no_lim)then
	  nd=0
	
	 else
	  print*,'On tient compte des discontinuites des derivees premieres,'
	  print*
	  print*,'Nombre de limites ZR/ZC:',lim,'. Indices des limites:',
     1	(jlim(i),i=1,lim),'.'
c	  print*,'lim',lim,(jlim(i),i=1,nd)
c	  print*,'dx=?'
c	  read*,dx
 
	  nd=0
	  do i=1,lim
	   do k=-2,2
	    nd=nd+1
	    i_dis(nd)=jlim(i)+k
	    x_g(nd)=i_dis(nd)-dx
	
c	    print*	
c	    print*,'x a gauche',i,jlim(i),nd,k,i_dis(nd),x_g(nd)
	    call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x_g(nd),lx,fx,dfxdx)
	    do j=1,nes
	     df_g(nes*(nd-1)+j)=dfxdx(j)	!a gauche de la discontinuite
	    enddo
c	    write(6,2000)(dfxdx(j),j=1,ne)
 
	    x_d(nd)=i_dis(nd)+dx
c	    print*,'x droite  ',i,jlim(i),nd,k,i_dis(nd),x_d(nd)
	    call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x_d(nd),lx,fx,dfxdx)
	    do j=1,nes
	     df_d(nes*(nd-1)+j)=dfxdx(j)	!a droite de la discontinuite
	    enddo	
c	    write(6,2000)(dfxdx(j),j=1,ne)
	   enddo	!k	
	  enddo		!i sur lim
	 endif	!no_lim
c	 pause'fin des discons'
	
c	 discontinuite de la derivee a la limite enveloppe/atmosphere	
 
c	 if(.false.)then	!pour eliminer eventuellement
	 if(.true.)then
	  if(atm)then
	   nd=nd+1
	   i_dis(nd)=n
	   print*
	   print*,'discontinuite de la derivee au fond de atm.',n,nd
	   print*
	
c	   a gauche cote enveloppe
	
	   x_g(nd)=n-dx	
c	   print*,'x a gauche',nd,i_dis(nd)
c	   write(6,2000)x_g(nd)
	   call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,x_g(nd),lx,
     1	fx,dfxdx)	!cote enveloppe	   	
	   do j=1,nes
	    df_g(nes*(nd-1)+j)=dfxdx(j)	!a gauche de la discontinuite
	   enddo
c	   write(6,2000)(df_g(nes*(nd-1)+i),i=1,nes)
 
c	   a droite cote atmosphere
 
	   x_d(nd)=n+dx
c	   print*,'x a droite  ',nd,i_dis(nd)
c	   write(6,2000)x_d(nd)	   	
	   call sbsp1dn(nea,bp_atm,q_atm,q_atmt,n_atm,mpr,knot_atm,.true.,
     1	1.d0+dx,lx,fx,dfxdx)	!cote atmosphere
	   do j=1,2	!pour ln p, ln t
	    df_d(nes*(nd-1)+j)=dfxdx(j)	!a droite de la discontinuite
	   enddo
	   df_d(nes*(nd-1)+3)=2.*fx(3)*dfxdx(3)	   	!pour r
	   df_d(nes*(nd-1)+4)=0		!pour l
	   df_d(nes*(nd-1)+5)=2.d0/3.d0*fx(5)**(-1.d0/3.d0)*dfxdx(5)	!pour m	   	   	
c	   write(6,2000)(df_d(nes*(nd-1)+i),i=1,nes)
c	   pause'fond atmosphere'
 
c	   discontinuite en n23
 
	   nd=nd+1
	   i_dis(nd)=n+n23-ndeb	!ndeb avant dernier point elimine
	
c	   print*,'discontinuite en n23',n,n23,ndeb,i_dis(nd)
	
	   print*
	   print*,'discontinuite de la derivee/q en n23',nd,i_dis(nd),n23,ndeb
	   print*
	
c	   a gauche
 
	   x_g(nd)=i_dis(nd)-dx
c	   print*,'x a gauche',nd,i_dis(nd)
c	   write(6,2000)x_g(nd)
	   call sbsp1dn(nea,bp_atm,q_atm,q_atmt,n_atm,mpr,knot_atm,.true.,
     1	dble(n23-dx),lx,fx,dfxdx)	!cote gauche
	   do j=1,2	!pour ln p, ln t
	    df_g(nes*(nd-1)+j)=dfxdx(j)	!a droite de la discontinuite
	   enddo
	   df_g(nes*(nd-1)+3)=2.*fx(3)*dfxdx(3)	   	!pour r
	   df_g(nes*(nd-1)+4)=0		!pour l
	   df_g(nes*(nd-1)+5)=2.d0/3.d0*fx(5)**(-1.d0/3.d0)*dfxdx(5)	!pour m
c	   write(6,2000)(df_g(nes*(nd-1)+i),i=1,nes)
 
c	   a droite
 
	   x_d(nd)=i_dis(nd)+dx
c	   print*,'x a droite  ',nd,i_dis(nd)
c	   write(6,2000)x_d(nd)
	   call sbsp1dn(nea,bp_atm,q_atm,q_atmt,n_atm,mpr,knot_atm,.true.,
     1	dble(n23+dx),lx,fx,dfxdx)	!cote droit
	   do j=1,2	!pour ln p, ln t
	    df_d(nes*(nd-1)+j)=dfxdx(j)	!a droite de la discontinuite
	   enddo
	   df_d(nes*(nd-1)+3)=2.*fx(3)*dfxdx(3)	   	!pour r
	   df_d(nes*(nd-1)+4)=0		!pour l
	   df_d(nes*(nd-1)+5)=2.d0/3.d0*fx(5)**(-1.d0/3.d0)*dfxdx(5)	!pour m	
c	   write(6,2000)(df_d(nes*(nd-1)+i),i=1,nes)	
c	   pause'en n23'
 
	  endif		!si atmosphere
	 endif		!.true./.false.	
c	 pause'initial'
 
c	 initialisation de l'interpolation par spline naturelle
c	 avec les derivees premieres connues aux extremites
c	 et discontinuite des derivees premieres aux limites ZR/ZC
 
c	 print*,nes,ns,nd,(i_dis(i),i=1,nd)
c	 write(6,2000)(x_g(i),x_d(i),i=1,nd)
c	 do j=1,nes
c	  write(6,2000)(df_d(nes*(i-1)+j),df_g(nes*(i-1)+j),i=1,nd)
c	 enddo
c	 pause'avant sn1dn_dis'
 
	 call sn1dn_dis(nes,f,qs,qst,ns,knots,dfdx1,dfdxn,nd,
     1	i_dis,x_g,df_g,x_d,df_d)
c	 pause'sn1dn_dis'
 
c	 if(.true.)then		!test
	 if(.false.)then
	
	  do i=1,ns		!test
	   x=qs(i)	
	   call sbsp1dn(nes,f,qs,qst,ns,4,knots,.true.,x,lx,fx,dfdx)
c	   write(6,*)i
c	   write(6,2000)(fx(j),j=1,nes)
c	   write(6,2000)(dfdx(j),j=1,nes)
	   p=exp(fx(1))
	   t=exp(fx(2))
	   r=sqrt(abs(fx(3)))
	   l=sqrt(abs(fx(4)))**3
	   m=sqrt(abs(fx(5)))**3
	   call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,fx(5),lx,xchim,dxchim)
	   do j=1,nbelem
	    if(fx(5) .gt. 0.d0)then
	     dxchim(j)=dxchim(j)*2.d0/3.d0/sqrt(fx(5))
	    else
	     dxchim(j)=0.
	    endif
	   enddo
	   call thermo_3(p,t,xchim,m,l,r,.true.,dxchim,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,0.d0,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)
	   write(6,2000)p,t,r,l,m,ro,xchim(1)
c	   pause'apres thermo test 2'
	  enddo
	  pause'test 3'
	 endif
	
	endif
 
c	cas ou on donne le rayon ou la masse
 
	if(m_ou_r .eq. 'masse')then
	 mi=m**(2.d0/3.d0)
	 call slinf(m,ms,ns,lx)		!entre lx et lx+1
	 q_int=(m-ms(lx+1))/(ms(lx+1)-ms(lx))+qs(lx+1)	!approx. lin.
	 corr=1.d30
	 ntour=0
c	 print*,mi,q_int
	 do while(ntour .lt. 30 .and. abs(corr) .gt. 1.d-12)
	
	  call sbsp1dn(nes,f,qs,qst,ns,4,knots,.true.,q_int,lx,fx,dfdx)
	  	
	  corr=(fx(5)-mi)/dfdx(5)		!mq=mu(qn) et corr=f(Xn)/f'(Xn)
	  q_int=q_int-corr			!NR: Xn+1=Xn-f(Xn)/f'(Xn)
	  q_int=max(1.d0,min(q_int,dfloat(ns)))
	  ntour=ntour+1
	 enddo
c	 write(6,*)ntour,corr
	 if(ntour .gt. 30)then
	  call slinf(m,ms,ns,lx)		!entre lx et lx+1
	  write(6,*)'intertot_3 en m: pas conv. nwt; m,lx,q_int,mi,fx(5),corr',
     1	 m,lx,q_int, mi, fx(5),corr
	  q_int=(m-ms(lx+1))/(ms(lx+1)-ms(lx))+qs(lx+1)	!approx. lin.
	  write(6,*)'interp. lin. q_int=',q_int
	 endif
	elseif(m_ou_r .eq. 'rayon')then
	 ri=r**2
	 call slinf(r,rs,ns,lx)		!entre lx et lx+1
	 q_int=(r-rs(lx+1))/(rs(lx+1)-rs(lx))+qs(lx+1)	!approx. lin.
	 corr=1.d30
	 ntour=0
	 do while(ntour .lt. 30 .and. abs(corr) .gt. 1.d-12)
	  call sbsp1dn(nes,f,qs,qst,ns,4,knots,.true.,q_int,lx,fx,dfdx)
	  corr=(fx(3)-ri)/dfdx(3)		!mq=mu(qn) et corr=f(Xn)/f'(Xn)
	  q_int=q_int-corr			!NR: Xn+1=Xn-f(Xn)/f'(Xn)
	  q_int=max(1.d0,min(q_int,dfloat(ns)))
	  ntour=ntour+1
	 enddo
c	 write(6,*)ntour,corr
	 if(ntour .gt. 30)then
	  call slinf(r,rs,ns,lx)		!entre lx et lx+1
	  write(6,*)'intertot_3 en r: pas conv. nwt; r,lx,q_int,ri,fx(3),corr',
     1	 r,lx,q_int, ri, fx(3),corr
	  q_int=(r-rs(lx+1))/(rs(lx+1)-rs(lx))+qs(lx+1)	!approx. lin.
	  write(6,*)'interp. lin. q_int=',q_int
	 endif
	endif
 
c	on retrouve la solution en q_int
 
c	print*,q_int,qs(1),qs(ns),ns
c	pause'q_int'
 
	call sbsp1dn(nes,f,qs,qst,ns,4,knots,.true.,q_int,lx,fx,dfdx)
	
c	print*,knots,qst(1),qst(knots)	
	
	p=exp(fx(1))
	t=exp(fx(2))
	r=sqrt(abs(fx(3)))
	l=sqrt(abs(fx(4)))**3
	m=sqrt(abs(fx(5)))**3
	grad_mj=dfdx(2)/dfdx(1)		!d ln P  /  d ln T
c	write(6,*)'ppm',p
c	write(6,2000)q_int,(fx(j),j=1,nes)
c	write(6,2000)q_int,p,t,r,l,m
c	write(6,2000)grad_mj,dfdx(2),dfdx(1)
 
c	la composition chimique interpolation en m**2/3
 
c	print*,n,ns,q_int,r
 
	if(fx(5) .le. mc(nc))then		!dans la structure interne
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,fx(5),lx,xchim,dxchim)
	 do i=1,nbelem
	  if(fx(5) .gt. 0.d0)then
	   dxchim(i)=dxchim(i)*2.d0/3.d0/sqrt(fx(5))
	  else
	   dxchim(i)=0.
	  endif
	 enddo
	else		!dans l'atmosphere
	 do i=1,nbelem	
	  xchim(i)=chim_atm(i)
	  dxchim(i)=0.d0
	 enddo
	endif
	
	if(iw .gt. 1 .and. fx(3) .ne. 0.)then		!rotation
	 w=xchim(iw)/fx(3)
	 dw=dxchim(iw)/fx(3)
	else
	 w=w_rot
	 dw=0.
	endif
 
c	print*,lim
c	write(6,2000)p,t,m,l,r,mstar
c	write(6,2000)(xchim(i),i=1,nbelem)
c	write(6,2000)(dxchim(i),i=1,nbelem)
 
	call thermo_3(p,t,xchim,m,l,r,.true.,dxchim,mstar,
     1	ro,drop,drot,drox,u,dup,dut,dux,lamb,r_zc,r_ov,lim,
     2	gradient,dgradp,dgradt,dgradl,dgradr,dgradm,dgradx,
     3	epsilon,depsp,depst,depsx,kap,dkapp,dkapt,dkapx,
     4	delta,deltap,deltat,deltax,cp,dcpp,dcpt,dcpx,
     5	hp,gradad,gradrad,gradconv,d_grad,w,dw,nh1,nhe1,nhe2,
     6	etat,opa,conv,nuc)
c	pause'apres thermo'
 
	alfa=p/ro*drop	!un peu de thermo
	beta=1.d0-aradias3*t**4/p
	
	do i=1,nchim
	 chimg(i)=xchim(i)
	 dchimg(i)=dxchim(i)
	enddo
	call chim_gram_3(chimg,dchimg,nuc_m)
	
	if(hp .le. 0.d0)then		!au centre
	 vaissala=0.
	
	else
	 	
	 vaissala=r*rsol/hp*delta*(gradad-gradient)
     1	-cte3*r**3*drox*dchimg(1)
	
	endif
 
c	write(6,*)'ro,drot,gradient,gradad,drox'
c	write(6,2000)ro,drot,gradient,gradad,drox
c	print*,r,q_int
c	pause'avant return'
 
	return
 
	end
