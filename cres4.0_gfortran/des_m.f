c*************************************************************************
 
	subroutine des_m(bp,q,n,qt,knot,chim,mc,nc,mct,knotc,teff,
     1	r2,m23,age,lim,lconv,m_zc,mstar)
	
c	dessin de variables en cours d'evolution
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	25 08 97 : mise en place des variables euleriennes
 
c entrees, significations evidentes
c	p,t,m,l,r,ro,grad_ad,grad_mj,alfa,delta,kap,cp,teff,age
c	epsilon(5,i): energie thermo et gravifique, voir cesam_3
c	n: nombre de couches
c	chim,mc,mct,nc,knotc: pour interpolation de la comp. chim.
c	lim: nombre de limites ZR/ZC
c	lconv: .true. si debut de ZC
c	m_zc: masse aux limites ZC/ZR
c	mstar: masse totale au temps du dessin
	
	implicit none
 
	include 'cesam_3.parametres'
	include 'modele_3.common'
	include 'evol_chim_3.common'
	include 'ctephy.common'
	
	integer pndes
	parameter(pndes=1000)
 
	integer i,ndes,nc,knotc,lx,j,lim,nxsub,nysub,long,n,knot,mm,pp

	external long

	real*8	chim(1),mc(1),mct(1),xdbl,teff,age,dc,dq,f(pne),dfdx(pne),
     1	xchim(pnelem),dxchim(pnelem),m_zc(1),qx(pndes),lib,
     2	bp(1),q(1),qt(1),nuc_m(pnelem),mstar,r2(1),m23(1)
	
	real*4	x(pndes),y(pnelem,pndes),yd(pndes),ymaxj(pnelem),
     1	xmin,xmax,ymin,ymax,inch,xleft,ybot,xmin1,xmax1,ymin1,ymax1,
     2	xtick,ytick,disp,coord,fjust,h,dh,ld,dl,pas,xloc,yloc(pnelem),
     3	xmax3,xmin3,ymax3,ymin3,xmin4,xmax4,lteff,llext,lteffp,llextp,
     4	ymin4,ymax4,x4(2),qd(pndes),pmax,tmax,rmax,lmax,mmax,yp(pndes),
     5	yt(pndes),yr(pndes),yl(pndes),ym(pndes),fac2,fac,fac1
 
	character*10 device,xopt,yopt,side
	character*15 htext(pnelem),char	
	character*23 tex(pnelem),ptex,ttex,rtex,ltex
	character*50 text
 
	logical init,lconv(1),efface, en_masse
 
	data init/.true./,efface/.false./
	data inch/2.54/,pas/-1./	
 
	save	
 
2000	format(1x,1p8d10.3)
 
	if(init)then
	 init=.false.
	
c	 fac facteur en r, si rtot > fac rtotp, on change d'echelle en R
 
	 fac=4./3.
	 fac1=.9*fac
	 fac2=.75
	
	 lib=2.5e6		!temperature de disparition du Li
	
	 device='/XSERVE'
c	 print*,'entrer le nom du device: /x, /wd, /SUN, /XINDOW, /PS ou ?'
c	 read(5,'(a)')device
c	 print*,device(:long(device))
c	 pause
	 if(device(:long(device)) .eq. '/PS')then
	  device=nom_fich2(:long(nom_fich2))//'.ps/PS'
	  print*,'nom du fichier du plot: ',nom_fich2(:long(nom_fich2))//'.ps'
	 endif
	 call pgbegin(0,device,1,1)
	 call pgscf(2)		!roman font
c	 call pgsch(1.2)	!hauteur des caracteres
	 call pgslw(3)		!epaisseur du trait
c	 call pgqcol(i,j)
c	 print*,'color',i,j
c	 call pgscr(0,.66,.66,.66)	!amenagement du fond (blanc)
c	 call pgscr(1,0.,0.,0.)	!amenagement du 1 (noir)	
	 call pgscr(4,.3,.5,1.)	!amenagement du bleu
		
	 ndes=pndes
c	 ndes=30
	 dc=mtot/ndes
 
	 h=7./inch	!hauteur du cadre
	 dh=1.5/inch	!saut entre cadres
	 ld=10./inch	!largeur du cadre
	 dl=2./inch	!espace en largeur
	 xleft=1.8/inch
	 ybot=1.8/inch
	
	 do i=1,pnelem
	  yloc(i)=1.-0.08*i
	 enddo
	
c	 premier cadre: composition chimique
 
	 xmin=-0.01
	 xmax=mtot*1.01
	 xloc=xmax*.65
	 ymax=1.05
	 ymin=.0	
	
	 call pgvsize(xleft,xleft+ld,ybot,ybot+h)
	 call pgwindow(xmin,xmax,ymin,ymax)
	 xopt='bcnst'
	 xtick=0.
	 nxsub=0
	 yopt='bcnst'
	 ytick=0.
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'	!abscisses
	 disp=2
	 coord=.5
	 fjust=.5
	 text='M en M\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)
	 side='l'	!ordonnees
	 disp=2
	 coord=.5
	 fjust=.5
	 text='Log10 X\\di'
	 call pgmtext(side,disp,coord,fjust,text)	
 
c	 second cadre: P, T, R, L, M en f(m)
 
	 xmin1=-.01
	 xmax1=xmax
	 ymin1=-.01
	 ymax1=1.1
	 call pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin1,xmax1,ymin1,ymax1)
	 xopt='bcnst'
	 xtick=0.
	 nxsub=0
	 yopt='bcnst'
	 ytick=0.
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'	!abscisses
	 disp=2
	 coord=.5
	 fjust=.5
	 text='M en M\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)
	 side='l'	!ordonnees
	 disp=2
	 coord=.5
	 fjust=.5
	 text='Log\\d10\\u P Log\\d10\\u T L/L\\d\\(2281)\\u R/R\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)	
c	 pause
	
	 side='t'	!sommet
	 disp=1
	 coord=.5
	 fjust=.5
	 call pgmtext(side,disp,coord,fjust,nom_fich2)
	
c	 troisieme cadre: diagramme HR
 
	 lteffp=log10(teff)
	 lx=n
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(n),lx,f,dfdx)
	 if(en_masse)then	
	  llextp=log10(f(4))*1.5
	 else
	  llextp=log10(f(4))	
	 endif
 
	 if(mtot .le. .95)then
	  xmin3=3.75
	  xmax3=3.5
	  ymax3=3.
	  ymin3=-1.
 
	 elseif(mtot .le. 1.5)then
	  xmin3=3.90
	  xmax3=3.5
	  ymax3=3.
	  ymin3=-0.5
	  	
	 elseif(mtot .le. 2.5)then
	  xmin3=4.05
	  xmax3=3.5
	  ymax3=3.5
	  ymin3=0.
	  	 	 	
	 elseif(mtot .le. 4.9)then
	  xmin3=4.3
	  xmax3=3.5
	  ymax3=4.5
	  ymin3=2.
 
	 else
	  xmin3=4.6
	  xmax3=3.5
	  ymax3=5.
	  ymin3=2.5
 
	 endif	
 
	 call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	 call pgwindow(xmin3,xmax3,ymin3,ymax3)
	 xopt='bcnst'
	 xtick=0.
	 nxsub=0
	 yopt='bcnst'
	 ytick=0.
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'	!abscisses
	 disp=2
	 coord=.5
	 fjust=.5
	 text='HR --- Log\\d10\\u T\\deff\\u'
	 call pgmtext(side,disp,coord,fjust,text)
	 side='l'	!ordonnees
	 disp=2
	 coord=.5
	 fjust=.5
	 text='HR --- Log\\d10\\u L/L\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)
 
c	 quatrieme cadre: ZC
 
	 xmin4=-3.
	 xmax4=600
	 ymin4=-.05
	 ymax4=mtot*1.05
	
	 call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin4,xmax4,ymin4,ymax4)
	 xopt='bcnst'
	 xtick=0.
	 nxsub=0
	 yopt='bcnst'
	 ytick=0.
	 nysub=0
	 call pgbox(xopt,xtick,nxsub,yopt,ytick,nysub)
	 side='b'	!abscisses
	 disp=2
	 coord=.5
	 fjust=.5
	 text='pas temporel'
	 call pgmtext(side,disp,coord,fjust,text)
	 side='l'	!ordonnees
	 disp=2
	 coord=.5
	 fjust=.5
	 text='ZC ---- M/M\\d\\(2281)\\u'
	 call pgmtext(side,disp,coord,fjust,text)	
c	 pause
	 side='t'	!sommet
	 disp=1
	 coord=.5
	 fjust=.5
	 call pgmtext(side,disp,coord,fjust,nom_fich2)
	 call pgslw(1)		!epaisseur du trait
	endif
 
c	les abondances
 
	if(efface)then
	 call pgsch(0.8)	!hauteur des caracteres
	
	 call pgvsize(xleft,xleft+ld,ybot,ybot+h)	
	 call pgwindow(xmin,xmax,ymin,ymax)
	
	 call pgsls(1)
	 call pgsci(0)
	 do j=1,nbelem
	  do i=1,ndes	
	   yd(i)=y(j,i)
	  enddo
	  call pgline(ndes,x,yd)
	  call pgtext(xloc,yloc(j),tex(j))
	 enddo
	 call pgsci(1)
	 call pgsch(1.)	!hauteur des caracteres
	endif		!fin d'effacage
	
	dc=mc(nc)/ndes
	do j=1,nbelem
	 ymaxj(j)=0
	enddo
	do i=1,ndes
	 xdbl=i*dc	
	 x(i)=sqrt(xdbl)**3
	 call sbsp1dn(nbelem,chim,mc,mct,nc,m_ch,
     1	knotc,.true.,xdbl,lx,xchim,dxchim)
	 call chim_gram_3(xchim,dxchim,nuc_m)
	 do j=1,nbelem
	  y(j,i)=max(xchim(j),1.d-30)	
	  ymaxj(j)=max(y(j,i),ymaxj(j))
	 enddo
	enddo
	do i=1,ndes
	 do j=1,nbelem
	  if(j .ne. 1 .and. j .ne. ihe4)y(j,i)=y(j,i)/ymaxj(j)		
	 enddo
	enddo
	call pgsls(1)
	call pgvsize(xleft,xleft+ld,ybot,ybot+h)	
	call pgwindow(xmin,xmax,ymin,ymax)
	
	call pgsch(0.8)	!hauteur des caracteres
	do j=1,nbelem
	 pp=log10(abs(ymaxj(j)))-2.	!pour avoir 2 chiffres
	 mm=abs(ymaxj(j))/10.**pp
	 call pgnumb(mm,pp,2,htext(j),15)
	 if(j .ne. 1 .and. j .ne. ihe4)then
	  tex(j)=nom_elem(j)//'='//htext(j)
	 else
	  tex(j)=nom_elem(j)//'= exact'
	 endif
c	 write(6,2000)ymaxj(j)
c	 print*,j,mm,pp,' ',htext(j),tex(j)
c	 pause	
	 do i=1,ndes
	  yd(i)=y(j,i)
	 enddo
	 call pgsci(j)
	 call pgline(ndes,x,yd)	
	 call pgtext(xloc,yloc(j),tex(j))
	enddo	
	call pgsci(1)		
	efface=.true.	
	call pgsch(1.)	!hauteur des caracteres
c	pause
 
c	P, T, R, L, R en fonction de M
 
	if(efface)then
	 call pgvsize(xleft,xleft+ld,ybot+h+dh,ybot+h+dh+h)
	 call pgwindow(xmin1,xmax1,ymin1,ymax1)
	 call pgsci(0)
	 call pgline(ndes,ym,yp)
	 call pgline(ndes,ym,yt)
	 call pgline(ndes,ym,yl)
	 call pgline(ndes,ym,yr)	
	 call pgsch(.8)	!hauteur des caracteres
	 call pgtext(xloc,yloc(1),ptex)
	 call pgtext(xloc,yloc(2),ttex)
	 call pgtext(xloc,yloc(3),ltex)
	 call pgtext(xloc,yloc(4),rtex)
	 call pgsch(1.)	!hauteur des caracteres	
	endif
	call pgsls(1)
	
	dq=(n-1.d0)/dble(ndes-1)	!pour les variables	
	do i=1,ndes
	 qx(i)=dq*(i-1)+1.d0
	 qd(i)=qx(i)
	enddo
	
	call pgsls(1)
	call pgsci(1)
	do i=1,ndes
	 call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,min(qx(i),q(n)),lx,f,dfdx)
	 yp(i)=f(1)
	 yt(i)=f(2)
	 if(en_masse)then	
	  yr(i)=sqrt(abs(f(3)))
	  yl(i)=sqrt(abs(f(4)))**3
	  ym(i)=sqrt(abs(f(5)))**3
	 else
	  yr(i)=abs(f(3))
	  yl(i)=f(4)
	  ym(i)=abs(f(5))
	 endif	 		
	enddo
	pmax=0
	tmax=0
	rmax=0
	lmax=0
	mmax=0
	do i=1,ndes
	 pmax=max(pmax,yp(i))
	 tmax=max(tmax,yt(i))	
	 rmax=max(rmax,yr(i))	
	 lmax=max(lmax,yl(i))
	 mmax=max(mmax,ym(i))
	enddo	
	do i=1,ndes
	 yp(i)=yp(i)/pmax	
	 yt(i)=yt(i)/tmax	
	 yr(i)=yr(i)/rmax	
	 yl(i)=yl(i)/lmax
c	 ym(i)=yr(i)/mmax
	enddo
	
	call pgsch(0.8)	!hauteur des caracteres
	
	pmax=exp(pmax)	
	pp=log10(pmax)-2.	!pour avoir 2 chiffres
	mm=pmax/10.**pp	
	call pgnumb(mm,pp,2,char,15)
	ptex='P\\dmax\\u='//char
	call pgsci(2)
	call pgline(ndes,ym,yp)
	call pgtext(xloc,yloc(1),ptex)
		
	tmax=exp(tmax)	
	pp=log10(tmax)-2.	!pour avoir 3 chiffres
	mm=tmax/10.**pp	
	call pgnumb(mm,pp,2,char,15)
	ttex='T\\dmax\\u='//char
	call pgsci(3)
	call pgline(ndes,ym,yt)
	call pgtext(xloc,yloc(2),ttex)
	
	pp=log10(lmax)-2.	!pour avoir 2 chiffres
	mm=lmax/10.**pp	
	call pgnumb(mm,pp,2,char,15)
	ltex='L\\dmax\\u='//char
	call pgsci(6)
	call pgline(ndes,ym,yl)
	call pgtext(xloc,yloc(3),ltex)
	
	pp=log10(rmax)-2.	!pour avoir 2 chiffres
	mm=rmax/10.**pp	
	call pgnumb(mm,pp,2,char,15)
	rtex='R\\dmax\\u='//char
	call pgsci(7)
	call pgline(ndes,ym,yr)
	call pgtext(xloc,yloc(4),rtex)
	
	call pgsch(1.)	!hauteur des caracteres
c	pause	
	call pgsci(1)	
	
c	diagramme HR
 
	lteff=log10(teff)
	lx=n
	call sbsp1dn(ne,bp,q,qt,n,mpr,knot,.true.,q(n),lx,f,dfdx)
	if(en_masse)then	
	 llext=log10(f(4))*1.5
	else
	 llext=log10(f(4))
	endif
	call pgvsize(xleft+ld+dl,xleft+2*ld+dl,ybot,ybot+h)
	call pgwindow(xmin3,xmax3,ymin3,ymax3)
	call pgmove(lteffp,llextp)
	call pgdraw(lteff,llext)
	lteffp=lteff
	llextp=llext
	
c	zones radiatives, convectives
 
c	print*,lim,(lconv(i),i=1,lim)
c	print*,(m_zc(i),i=1,lim),age
 
	pas=pas+1.
	x4(1)=pas
	x4(2)=pas
 
	if(lim .gt. 0)then	!lim=0: modele totalement radiatif pas de dessin	
	 call pgvsize(xleft+ld+dl,xleft+2.*ld+dl,ybot+h+dh,ybot+2.*h+dh)
	
c	 print*,lim,(lconv(i),i=1,lim)
c	 write(6,2000)(m_zc(i),i=1,lim),mstar,mtot
 
	 call pgwindow(xmin4,xmax4,ymin4,ymax4)
	 if(lim .eq. 1 .and. .not.lconv(1) .and.
     1	m_zc(1)/mstar .ge. .99)then!modele totalement convectif
	  yd(1)=0.
	  yd(2)=mstar
	  call pgline(2,x4,yd)
	 else
	  i=1
	  do while(i .le. lim)
	   if(lconv(i))then	!debut de ZC
	    yd(1)=m_zc(i)
	    if(i .eq. lim)then		!ZC externe
	     yd(2)=mstar
	    else
	     yd(2)=m_zc(i+1)
	     i=i+1
	    endif	
	    call pgline(2,x4,yd)
	   else
	    if(i .eq. 1)then			!fin de ZC centrale
	     yd(1)=0
	     yd(2)=m_zc(1)	
	     call pgline(2,x4,yd)
	    else
	     print*,'probleme'
	     call pgend
	     stop
	    endif
	   endif
	   i=i+1
	  enddo	!while
	 endif
	endif
	
	return
	
	end
 
 
