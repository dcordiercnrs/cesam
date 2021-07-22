 
c***********************************************************************
 
	subroutine t6rinterp(slr,slt)
 
c	adaptation a CESAM3 de la routine t6rinterp du package OPAL_EOS
 
c	Auteur: P. Morel, Departement J.D. Cassini, O.C.A., Observatoire de Nice
 
c	The purpose of this subroutine is to interpolate in T6 and rho
 
	implicit none
 
	integer mx,mv,nr,nt
	parameter (mx=5,mv=10,nr=121,nt=167)
 
	integer	iu,is,kx,iw
	
	real*4	slr,slt,dix,dix2,esact2,esactq,esactq2,quad
 
	real*4  epl(mx,nt,nr),xx(mx)
	common/ee/epl,xx
	
	real*4	q(4),h(4),xxh
	common/aa/q,h,xxh
 
	integer	m,mf
	real*4	xz(mx,mv,nt,nr),t6list(nr,nt),rho(nr),t6a(nt),esk(nt,nr),
     1	esk2(nt,nr),dfsx(mx),dfs(nt),dfsr(nr),xa(mx)
	common/a/xz,t6list,rho,t6a,esk,esk2,dfsx,dfs,dfsr,m,mf,xa
 
	integer l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
	common/bb/l1,l2,l3,l4,k1,k2,k3,k4,ip,iq
 
	real*4	esact,eos(mv)	
	common/e/esact,eos
 
	save
	
      iu=0
      is=0
 
      do kx=k1,k1+ip
          iw=1
        iu=iu+1
        h(iu)=quad(is,iw,slr,esk(kx,l1),esk(kx,l2),esk(kx,l3),
     x  rho(l1),rho(l2),rho(l3))
          if(iq. eq. 3) then
            iw=2
            q(iu)=quad(is,iw,slr,esk(kx,l2),esk(kx,l3),esk(kx,l4),
     x      rho(l2),rho(l3),rho(l4))
          endif
        is=1
      enddo
c
      is=0
      iw=1
c..... eos(i) in lower-right 3x3(i=i1,i1+2 j=j1,j1+2)
      esact=quad(is,iw,slt,h(1),h(2),h(3),t6a(k1),t6a(k2),t6a(k3))
        if (iq. eq. 3) then
c.....    eos(i) upper-right 3x3(i=i1+1,i1+3 j=j1,j1+2)
          esactq=quad(is,iw,slt,q(1),q(2),q(3),t6a(k1),t6a(k2),t6a(k3))
        endif
        if(ip .eq. 3) then
c.....    eos(i) in lower-left 3x3.
          esact2=quad(is,iw,slt,h(2),h(3),h(4),t6a(k2),t6a(k3),t6a(k4))
c.....    eos(i) smoothed in left 3x4
          dix=(t6a(k3)-slt)*dfs(k3)
          esact=esact*dix+esact2*(1.-dix)
c       endif   ! moved to loc a
        if(iq .eq. 3) then
 
c.....     eos(i) in upper-right 3x3.
          esactq2=quad(is,iw,slt,q(2),q(3),q(4),t6a(k2),t6a(k3),t6a(k4))
          esactq=esactq*dix+esactq2*(1.-dix)
        endif
        endif  ! loc a
c
        if(iq .eq. 3) then
          dix2=(rho(l3)-slr)*dfsr(l3)
            if(ip .eq. 3) then
c.....        eos(i) smoothed in both log(T6) and log(R)
              esact=esact*dix2+esactq*(1-dix2)
            endif
        endif
        if (esact .gt. 1.e+15) then
          write(*,'("Interpolation indices out of range",
     x              ";please report conditions.")')
          stop
        endif
 
      return
 
      end
