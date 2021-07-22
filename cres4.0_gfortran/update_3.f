 
c*******************************************************************
 
	subroutine update_3(next,bp,q,n,qt,knot,chim,mc,nc,mct,knotc,
     1	r2,m23,dts,m_zc,r_zc,r_ov,lim,jlim,lconv,mstar,age,
     2	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,
     3	r2_t,m23_t,dt,m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,
     4	mstar_t,bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,
     5	r2_,m23_,dt_,bp__,q__,n__,qt__,knot__,chim__,mc__,
     6	nc__,mct__,knotc__,r2__,m23__,dt__,
     7	old_m23,new_m23,nm,new_m23t,knotm,
     8	old_m23_,new_m23_,nm_,new_m23t_,knotm_,
     9	old_m23__,new_m23__,nm__,new_m23t__,knotm__)	
	
c	passage au pas temporel suivant ou reinitialisation
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	version 3
 
c entree
c   next=.true.	 ***_---> ***__, ***_t--->***_, ***--->***_t,
c		 dt_---> dt__, dt-->dt_, dts-->dt : on decale d'un dt
 
c	 .false. ***_t--->*** , dt-->dt/1.2, reinitialisation
 
c solution, age=t+dt :
c	bp,q,n,qt,knot,chim,mc,nc,mct,knotc,r2,m23,dts,
c	m_zc,r_zc,r_ov,lim,jlim,lconv,mstar
c solution, age=t :
c	bp_t,q_t,n_t,qt_t,knot_t,chim_t,mc_t,nc_t,mct_t,knotc_t,r2_t,m23_t,dt,
c	m_zc_t,r_zc_t,r_ov_t,lim_t,jlim_t,lconv_t,mstar_t
c	entre t+dt et t
c	old_m23,new_m23,nm,new_m23t,knotm : interpolation m(t+dt)-->m(t)
 
c solution, age=t-dt_ :
c	bp_,q_,n_,qt_,knot_,chim_,mc_,nc_,mct_,knotc_,r2_,m23_,dt_,
c	entre t et t-dt_
c	old_m23_,new_m23_,nm_,new_m23t_,knotm_ : interpolation m(t+dt)-->m(t)
 
c solution, age=t-dt_-dt__:
c	bp__,q__,n__,qt__,knot__,chim__,mc__,nc__,mct__,knotc__,r2__,m23__,dt__
c	entre t et t-dt_-dt__
c	old_m23__,new_m23__,nm__,new_m23t__,knotm__
 
 
c	dts: pas temporel estime pour integration entre age et age+1
c	dt : pas temporel utilise entre age-1 et age
c	dt_ : pas temporel utilise entre age-2 et age-1
c	dt__ : pas temporel utilise entre age-2 et age-3
c	lconv=.true. : debut de ZC
 
c sortie
c	bp: solution
c	dt: nouveau pas temporel
c	lnp$,lnt$,r2$,l23$,m23$,q_t,qt$_t,knot$_t,nt,mu,dt: solution en t
c	m_t,mt_t,n_t,knot_t,xchims_t: idem en comp. chim.
c	dt_,nt_,m_,mt_,n_,knot_,xchims_: idem en t-1
c	m_zct : masses des limites ZR/ZC au temps t
c	lim_t : nombre de limites ZR/ZC au temps t
c	jlim_t : limites ZR/ZC au temps t
c	lconv_t=.true. : debut de ZC temps t
 
	implicit none
 
	include 'cesam_3.parametres'
	include 'ctephy.common'
	include 'modele_3.common'
	include 'evol_chim_3.common'
 
	integer n,knot,nc,knotc,lim,jlim(1),
     1	n_t,knot_t,nc_t,knotc_t,lim_t,jlim_t(1),
     2	n_,knot_,nc_,knotc_,n__,knot__,nc__,knotc__,i,
     3	nm,nm_,nm__,knotm,knotm_,knotm__
		
	real*8	bp(1),q(1),qt(1),chim(1),mc(1),mct(1),
     1	r2(1),m23(1),dts,m_zc(1),r_zc(1),r_ov(1),mstar,age,
     2	bp_t(1),q_t(1),qt_t(1),chim_t(1),mc_t(1),mct_t(1),
     3	r2_t(1),m23_t(1),dt,m_zc_t(1),r_zc_t(1),r_ov_t(1),mstar_t,
     4	bp_(1),q_(1),qt_(1),chim_(1),mc_(1),mct_(1),
     5	r2_(1),m23_(1),dt_,bp__(1),q__(1),qt__(1),chim__(1),
     6	mc__(1),mct__(1),r2__(1),m23__(1),dt__,
     7	old_m23(1),new_m23(1),new_m23t(1),
     8	old_m23_(1),new_m23_(1),new_m23t_(1),	
     9	old_m23__(1),new_m23__(1),new_m23t__(1)	
			
	logical next,lconv(1),lconv_t(1)
 
2000	format((1x,1p8d10.3))
 
c	print*,'update', dt,dt_,dt__
 
	if(next)then
	 if(age .gt. 0.d0)then
	  if(dt_ .gt. 0.d0)then !***_--> ***__
	   do i=1,pbp
	    bp__(i)=bp_(i)
	   enddo		!i
	   n__=n_
	   do i=1,n_
	    q__(i)=q_(i)
	   enddo	
	   knot__=knot_	
	   do i=1,knot_
	    qt__(i)=qt_(i)
	   enddo
	   do i=1,n_
	    r2__(i)=r2_(i)
	    m23__(i)=m23_(i)
	   enddo
	   nc__=nc_			!pour la comp. chim.
	   do i=1,nc_
	    mc__(i)=mc_(i)
	   enddo	!i
	   knotc__=knotc_
	   do i=1,knotc_
	    mct__(i)=mct_(i)
	   enddo
	   do i=1,pchim
	    chim__(i)=chim_(i)
	   enddo
	   dt__=dt_
	  endif	
	
	  nm__=nm_
	  do i=1,nm_
	   old_m23__(i)=old_m23_(i)
	   new_m23__(i)=new_m23_(i)
	  enddo
	  knotm__=knotm_
	  do i=1,knotm_
	   new_m23t__(i)=new_m23t_(i)
	  enddo
	
	  do i=1,pbp		!***_t--> ***_
	   bp_(i)=bp_t(i)
	  enddo		!i
	  n_=n_t
	  do i=1,n_t
	   q_(i)=q_t(i)
	  enddo	
	  knot_=knot_t	
	  do i=1,knot_t
	   qt_(i)=qt_t(i)
	  enddo
	  do i=1,n_t
	   r2_(i)=r2_t(i)
	   m23_(i)=m23_t(i)
	  enddo
	  nc_=nc_t			!pour la comp. chim.
	  do i=1,nc_t
	   mc_(i)=mc_t(i)
	  enddo	!i
	  knotc_=knotc_t
	  do i=1,knotc_t
	   mct_(i)=mct_t(i)
	  enddo
	  do i=1,pchim
	   chim_(i)=chim_t(i)
	  enddo
	  dt_=dt
	
	 endif
	
	 do i=1,pbp		!*** ---> ***_t
	  bp_t(i)=bp(i)
	 enddo
	 n_t=n
	 do i=1,n
	  q_t(i)=q(i)
	 enddo	
	 knot_t=knot	
	 do i=1,knot
	  qt_t(i)=qt(i)
	 enddo
	 do i=1,n
	  r2(i)=bp(ne*(i-1)*m_qs+3)		!extraction de r2
	  m23(i)=bp(ne*(i-1)*m_qs+5)		!extraction de m23	
	  r2_t(i)=r2(i)
	  m23_t(i)=m23(i)
	 enddo
	 nc_t=nc			!pour la comp. chim.
	 do i=1,nc
	  mc_t(i)=mc(i)
	 enddo	!i
	 knotc_t=knotc
	 do i=1,knotc
	  mct_t(i)=mct(i)
	 enddo
	 do i=1,pchim
	  chim_t(i)=chim(i)
	 enddo
	 lim_t=lim
	 do i=1,lim
	  m_zc_t(i)=m_zc(i)
	  r_zc_t(i)=r_zc(i)
	  r_ov_t(i)=r_ov(i)	
	  lconv_t(i)=lconv(i)
	  jlim_t(i)=jlim(i)
	 enddo	
	 mstar_t=mstar
	
	 nm_=nm
	 do i=1,nm
	  old_m23_(i)=old_m23(i)
	  new_m23_(i)=new_m23(i)
	 enddo
	 knotm_=knotm
	 do i=1,knotm
	  new_m23t_(i)=new_m23t(i)
	 enddo
	 	
	 if(age .eq. 0.d0)then
	  dt_=0
	  dts=dt
	 else
	  dt_=dt
	 endif
	
	 dt=max(dtmin,min(dtmax,dts,dt*1.2,agemax-age))	!le nouveau dt
	
c	 write(6,2000)dtmin,dtmax,dts,age,agemax,dt
c	 write(6,2000)(q(i),i=1,n)
c	 pause'sortie update'
 
	else
c	 write(6,*)'UPDATE(2) dtmin,dtmax,dts,age,agemax,dt'
c	 write(6,2000)dtmin,dtmax,dts,age,agemax,dt
 
	 dt=max(dtmin,min(dt,agemax-age))	!reinitialisation
	
	 do i=1,pbp
	  bp(i)=bp_t(i)
	 enddo		!i
	 n=n_t
	 do i=1,n_t
	  q(i)=q_t(i)
	 enddo	
	 knot=knot_t	
	 do i=1,knot_t
	  qt(i)=qt_t(i)
	 enddo
	 do i=1,n
	  r2(i)=r2_t(i)
	  m23(i)=m23_t(i)
	 enddo
	 nc=nc_t			!pour la comp. chim.
	 do i=1,nc_t
	  mc(i)=mc_t(i)
	 enddo	!i
	 knotc=knotc_t
	 do i=1,knotc_t
	  mct(i)=mct_t(i)
	 enddo
	 do i=1,pchim
	  chim(i)=chim_t(i)
	 enddo
	 lim=lim_t
	 do i=1,lim_t
	  m_zc(i)=m_zc_t(i)
	  r_zc(i)=r_zc_t(i)
	  r_ov(i)=r_ov_t(i)	
	  lconv(i)=lconv_t(i)
	 enddo	
	 mstar=mstar_t	
	endif
	
c	write(6,*)'sortie UPDATE_3 next,dt,dt-',next,dt,dt_,dt__
 
	return
 
	end
 
