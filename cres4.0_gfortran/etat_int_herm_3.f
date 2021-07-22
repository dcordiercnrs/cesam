 
c****************************************************************
 
	subroutine etat_int_herm_3(p,t,xchim,deriv,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
 
c	interpolation de l'equation d'etat par interpolation Hermite 2D
 
c	Auteur: P. Morel Observatoire de NICE
c	version 3
 
c entree :
c	p : pression
c	t : temperature
c	xchim : composition chimique
c	deriv=.true. : calcul des derivees secondes
 
c sortie :
c	ro : densite et derivees
c	u : energie interne et derivees
c	nh1, nhe1, nhe2 : taux d'ionisation
c	lamb: degenerescence
 
	implicit none
 
	integer pnt,pnp,pntp,pnx,pntot
	parameter (pnt=100, pnp=180, pntp=40, pnx=5)
	parameter (pntot=pnp*pntp*pnx)
 
	include 'modele_3.common'
	include 'ctephy.common'
 
	integer np,nt,nx,i,j,pt(pnp),ntp,k,lx,ly,ind(2)
 
	real*8 p,t,xchim(1),ro,drop,drot,drox,drott,drotp,drotx,z,
     1	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb
 
	real*8 p_t(pnp),t_t(pnt),x(pnx),dx(pnx),ld(pnx),li(pnx),lip(pnx),der,
     3	ro_t(pntot),drop_t(pntot),drot_t(pntot),hx,h2x,h3x,tmax,pmax,
     4	drotp_t(pntot),u_t(pntot),dup_t(pntot),hy,h2y,h3y,tmin,pmin,
     5	dut_t(pntot),dutp_t(pntot),roz(pnx),drozp(pnx),herx(0:2,2,2),
     6	d2rozpt(pnx),drozt(pnx),d2rozt2(pnx),uz(pnx),ln10,lt,lp,
     7	duzp(pnx),d2uzpt(pnx),duzt(pnx),d2uzt2(pnx),xchi,hery(0:2,2,2)
 
	character*30 equation
 
	logical init,deriv,hors
c	data init/.true./
 
	character*40 fich(8)
	common/fich_etat/fich
 
	data init/.true./
 
2000	format((1x,1p8d10.3))
 
	if(init)then
	 init=.false.
	 open(unit=14,form='unformatted',status='old',file=fich(1))
	 read(14)ro_t,drop_t,drot_t,drotp_t,u_t,dup_t,dut_t,dutp_t,p_t,t_t,x,
     1	z,pt,ntp,np,nt,nx,equation
	 write(6,10)equation
	 write(2,10)equation
10	 format(1x,/,1x,'equation d''etat par interp. Hermite 2D de l''',
     1	a30,//)
	 close(unit=14)
	 if(z .ne. z0)then
	  write(6,111)z,z0
	  write(2,111)z,z0	
111	  format(1x,'ATTENTION, ATTENTION, ATTENTION, ATTENTION, ATTENTION',/,
     1	1x,'Z equation d''etat=',1pd10.3,' # Z0',1pd10.3,/,
     2	1x,' on utilisera Z0')
	 endif
	
	 ln10=log(10.d0)
c	 write(6,*)'x,z,np,nt,nx,ntp,equation',x,z,np,nt,nx,ntp,equation
	 do i=1,nx		!denominateurs de Lagrange
	  ld(i)=1
	  do j=1,nx
	   if(j .ne. i)ld(i)=ld(i)*(x(i)-x(j))
	  enddo
	 enddo
	 hx=p_t(2)-p_t(1)
	 hx=1./hx
	 h2x=hx**2
	 h3x=hx**3
	 hy=t_t(2)-t_t(1)
	 hy=1./hy
	 h2y=hy**2
	 h3y=hy**3
	 lx=10
	 ly=10
	 tmax=exp(t_t(pt(np-1)+ntp-1))
	 tmin=exp(t_t(1))
	 write(6,11)tmin,tmax
	 write(2,11)tmin,tmax
11	 format(1x,'utilisation de ETAT_EFF si T hors de [ ',1pd10.3,
     1	' , ',1pd10.3,' ]')
	 pmax=exp(p_t(np))
	 pmin=exp(p_t(1))
	 write(6,12)pmin,pmax
	 write(2,12)pmin,pmax
12	 format(1x,'utilisation de ETAT_EFF si P hors de [ ',1pd10.3,
     1	' , ',1pd10.3,' ]')
c	 pause
	endif
 
	hors=t .ge. tmax .or. t .le. tmin .or. p .ge. pmax .or. p .le. pmin
 
	if(hors)then
	 call etat_eff_ps_3(p,t,xchim,deriv,
     1	ro,drop,drot,drox,drott,drotp,drotx,
     2	u,dup,dut,dux,dutt,dutp,dutx,nh1,nhe1,nhe2,lamb)
	else		!base de Lagrange en X
	 if(xchim(1) .lt. 0.d0)then
	  xchi=0.
	 else
	  xchi=xchim(1)
	 endif
	 do i=1,nx
	  dx(i)=xchi-x(i)
	 enddo
	 do i=1,nx		!base de Lagrange
	  li(i)=1.d0
	  do j=1,nx
	   if(j .ne. i)li(i)=li(i)*dx(j)
	  enddo
	  li(i)=li(i)/ld(i)
	 enddo
	 do i=1,nx		!derivee de la base de Lagrange
	  lip(i)=0.d0
	  do k=1,nx
	   if(k .ne. i)then
	    if(dx(k). ne. 0.d0)then
	     lip(i)=lip(i)+li(i)/dx(k)
	    else
	     der=1
	     do j=1,nx
	      if(j .ne. i .and. j. ne. k)der=der*dx(j)
	     enddo
	     lip(i)=lip(i)+der/ld(i)
	    endif
	   endif
	  enddo
	 enddo
 
c	 write(6,*)'nx,xchi/x,li,lip',nx,xchi
c	 write(6,2000)(li(i),i=1,nx)
c	 write(6,2000)(lip(i),i=1,nx)
c	 write(6,2000)(x(i),i=1,nx)
	
	 lp=log(p)
	 lt=log(t)
c	 write(6,*)'lt,lp'
c	 write(6,2000)lt,lp
 
c	 recherche de l'encadrement
 
	 call slinf(lp,p_t,np,lx)
	 call slinf(lt,t_t,nt,ly)
	 call base_her_cte(p_t,lp,lx,herx,hx,h2x,h3x)
	 call base_her_cte(t_t,lt,ly,hery,hy,h2y,h3y)
 
c	 write(6,*)'lx,ly/p,t',lx,ly
c	 write(6,2000)lp,lt
c	 write(6,*)'base herx'
c	 do i=0,2
c	  write(6,2000)((herx(i,j,k),j=1,2),k=1,2)
c	 enddo		!i
c	 write(6,*)'base hery'
c	 do i=0,2
c	  write(6,2000)((hery(i,j,k),j=1,2),k=1,2)
c	 enddo		!i
 
	 ind(1)=ly-pt(lx)+1
	 ind(2)=ly-pt(lx+1)+1
 
	 do i=1,2
	  if(ind(i) .lt. 1)then
	   write(6,2001)lt/ln10,t_t(pt(lx+i-1))/ln10,p_t(lx+i-1)/ln10,lp/ln10,
     1	lx,ly,ind(i),i,pt(lx+i-1)
2001	   format(1x,'HERM_ETAT: hors limite en T: Log(T)=',1pd10.3,
     1	' < Log(T1i)=',1pd10.3,/,1x,' pour Log(Pi)=',1pd10.3,
     2	' Log(P)=',1pd10.3,' resultats incertains',/,1x,
     3	'lx=',i3,' ly=',i3,' ind(i)=',i3,' i=',i3,' pt(lx+i-1)=',i3)
	   ind(i)=1
c	   pause
	  elseif(ind(i) .ge. ntp)then
	   write(6,2002)lt/ln10,t_t(pt(lx+i-1)+ntp-1)/ln10,p_t(lx+i-1)/ln10,
     1	lp/ln10,lx,ly,ind(i),i,pt(lx+i-1),ntp
2002	   format(1x,'HERM_ETAT: hors limite en T: Log(T)=',1pd10.3,
     1	' > Log(Tni)=',1pd10.3,/,1x,' pour Log(Pi)=',1pd10.3,
     2	' Log(P)=',1pd10.3,' resultats incertains',/,1x,
     3	'lx=',i3,' ly=',i3,' ind(i)=',i3,' i=',i3,' pt(lx+i-1)=',i3,
     4	' ntp=',i3)
	   ind(i)=ntp-1
	  endif
	 enddo
 
	 call herm_etat(ro_t,drop_t,lx,roz,drozp,d2rozpt,herx,hery,
     1	nx,np,ly,drot_t,drotp_t,drozt,d2rozt2,deriv,ind)
	 call herm_etat(u_t,dup_t,lx,uz,duzp,d2uzpt,herx,hery,
     1	nx,np,ly,dut_t,dutp_t,duzt,d2uzt2,deriv,ind)
 
c	 write(6,*)'l,lx,ly',l,lx,ly
 
	 nh1=-100
	 nhe1=-100
	 nhe2=-100
	 lamb=5		!degenerescence non interpolee
 
c	 write(6,*)'roz,uz'
c	 write(6,2000)(roz(i),i=1,nx)
c	 write(6,2000)(drozp(i),i=1,nx)
c	 write(6,2000)(drozt(i),i=1,nx)
c	 write(6,2000)(uz(i),i=1,nx)
c	 write(6,2000)(duzp(i),i=1,nx)
c	 write(6,2000)(duzt(i),i=1,nx)
 
	 ro=0.d0		!interpolation en X
	 drop=0.d0
	 drot=0.d0
	 drox=0.d0
 
	 u=0.d0
	 dup=0.d0
	 dut=0.d0
	 dux=0.d0
 
	 drott=0.d0
	 drotp=0.d0
	 drotx=0.d0
 
	 dutt=0.d0
	 dutp=0.d0
	 dutx=0.d0
 
	 do i=1,nx
	  ro=ro+roz(i)*li(i)
	  drop=drop+drozp(i)*li(i)
	  drot=drot+drozt(i)*li(i)
	  drox=drox+roz(i)*lip(i)
 
	  u=u+uz(i)*li(i)
	  dup=dup+duzp(i)*li(i)
	  dut=dut+duzt(i)*li(i)
	  dux=dux+uz(i)*lip(i)
 
	  if(deriv)then
	   drott=drott+d2rozt2(i)*li(i)
	   drotp=drotp+d2rozpt(i)*li(i)
	   drotx=drotx+drozt(i)*lip(i)
 
	   dutt=dutt+d2uzt2(i)*li(i)
	   dutp=dutp+d2uzpt(i)*li(i)
	   dutx=dutx+duzt(i)*lip(i)
	  endif
	 enddo
 
c	 write(6,*)'ro,u'
c	 write(6,2000)ro,drop,drot,u,dup,dut
 
c	 tranformation inverse (theta,pi)-->(T,P)
 
	 ro=exp(ro)
	 drop=ro/p*drop
	 drot=ro/t*drot
	 drox=ro*drox
	 u=exp(u)
	 dup=u/p*dup
	 dut=u/t*dut
	 dux=u*dux
	 if(deriv)then
	  drott=(ro/t*drott+drot*(t/ro*drot-1.d0))/t
	  drotp=ro/p/t*drotp+drot*drop/ro
	  drotx=ro/t*drotx+drox*drot/ro
 
	  dutt=(u/t*dutt+dut*(t/u*dut-1.d0))/t
	  dutp=u/p/t*dutp+dut*dup/u
	  dutx=u/t*dutx+dux*dut/u
	 endif
c	 write(6,2004)p,t,ro,u,xchim(1),lx,ly,ind(1),ind(2),pt(lx+i-1)
2004	 format(1x,1p5d10.3,5i5)
	endif
 
	return
 
	end
