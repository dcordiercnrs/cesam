c
c*******************************************************************************
c
	subroutine hmnion(tl, eah, ehm, xhm, secder)
c
c  given ratio between succesive states of ionization, and
c  derivatives, for H in eah(1-10)
c  sets fraction of H-, and derivatives, into xhm(1-10)
c
c  Note: assumes that fraction of H- is always very small.
c
      implicit double precision (a-h, o-z)
	implicit integer(i-n)
      logical secder
      dimension eah(10), ehm(10),xhm(10), eahs(10)
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/ln10c/ amm,amm2,amm3
	common/dmudec/ dmu,dmuf,dmut,dmux,dmuff,dmuft,dmufx,dmutt,dmutx,
     1	dmuxx,idmu
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
      if(idgeos.ge.3) write(6,*) 'Entering hmnion'
      t=10.d0**tl
      tkev=ck1*t
c
      call zeroc(xhm,10)
c
      ext=exhm/tkev
      eea=ext-dmu
c
c  test for no h-
c
      if(eea.le.-100) then
        if(idgeos.ge.3) write(6,*) 'No H-'
        return
      end if
      eea=0.5*exp(eea)
      dnm=1+eah(1)
      xhm(1)=eea/dnm
c
c  to avoid over- and underflow, divide derivatives
c  of ea by dnm (modification 3/1/84)
c
      if(secder) then
        iders=10
      else
       iders=4
      end if
c
      do 10 i=2,iders
   10 eahs(i)=eah(i)/dnm
c
      ehm(1)=eea
      ehm(2)=-dmuf*eea
      ehm(3)=-(amm*ext+dmut)*eea
      ehm(4)=-dmux*eea
      if(secder) then
        ehm(5)=-dmuff*eea-dmuf*ehm(2)
        ehm(6)=-dmuft*eea-dmuf*ehm(3)
        ehm(7)=-dmufx*eea-dmux*ehm(2)
        ehm(8)=(amm2*ext-dmutt)*eea-(amm*ext+dmut)*ehm(3)
        ehm(9)=-dmutx*eea-dmux*ehm(3)
        ehm(10)=-dmuxx*eea-dmux*ehm(4)
      end if
c
c  derivatives of xh-
c
      ii=4
      do 20 l=1,3
      l1=l+1
      ehml=ehm(l1)
      eal=eahs(l1)
      xhm(l1)=(ehml-eea*eal)/dnm
      if(secder) then
        do 15 m=l,3
        m1=m+1
        ii=ii+1
   15   xhm(ii)=(ehm(ii)-(ehml*eahs(m1)+ehm(m1)*eal+eea*eahs(ii))
     *    +2*eea*eal*eahs(m1))/dnm
      end if
   20 continue
      if(idgeos.ge.3) write(6,*) 'x(H-) =',xhm(1)
      return
      end
