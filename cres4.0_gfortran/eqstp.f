 
c******************************************************************
 
	subroutine eqstp(pl,tl,x,y,z,fl)
 
c	routine pour EFF, Auteur J. Christensen Dalsgaard, adaptation:
c	A. Baglin, M. Auvergne, P. Morel
 
c  this routine iterates s/r eqstf in order to find the appropriate
c  fl value. at the end, eqstf is called to compute td-quantities.
c
	implicit real*8 (a-h,o-z)
 
      logical nosdit,notdit
      common/eqstd/ xii1(4),ane(10),rho(20),ht(20),pt(20),cp(4),dad(4),
     *  dlt(4),gm1,tprh,trhp,rhxp
      common/consts/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     *  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      dimension flini(2),wt(2)
	logical ok
	common /marche/ok
 
	integer iter,icase,iextr,i,iti
 
c
c
      iter = 15
      eps = 1.d-9
      epsf = 1.d-3
      carsav=car/3.
      beta=1.-carsav*10.**(4.*tl-pl)
      beta=max(beta,1.0d-6)
      pgl0=pl + log10(beta)
c...... initial guess of fl
c
c================ icase=1 for gas pressure as argument
      icase=1
      call inteff(tl,pgl0,rlbid,flini,wt,icase,iextr)
      flin=wt(1)*flini(1)+wt(2)*flini(2)
      fl=flin
      nosdit = .true.
      notdit = .true.
      car    = 0.
c........note that for the iteration no second derivatives are needed.
c........car is switched off temporarily to account for no radiation.
      call eqstf(fl,tl,x,y,z,nosdit,notdit)
c
      plt=log10(pt(1))
      s=pgl0-plt
      dfl=1.d1*sign(epsf,s)
c
c  iterate to desired pressure
      fl=fl+dfl
      do 10 i=1,iter
      iti = i
      call eqstf(fl,tl,x,y,z,nosdit,notdit)
      pln=log10(pt(1))
      if(i.eq.1) dlpdlf=(pln-plt)/dfl
      if(abs(dfl).lt.epsf) go to 5
      dlpdlf=(pln-plt)/dfl
      plt=pln
    5 dfl=(pgl0-pln)/dlpdlf
      fl=fl+dfl
      rdif=abs(dfl/fl)
      if(rdif.lt.eps) go to 20
   10 continue
      print 110,pgl0,tl,rdif
	ok=.false.
 20   nosdit = .false.
      notdit = .true.
c........this call is here to provide the second derivatives.
c........for this car is again switched on to account for radiation.
      car = carsav*3.
      call eqstf(fl,tl,x,y,z,nosdit,notdit)
      return
110   format(' >>>>>>>>> warning from eqstp: number of iterations',
     1	' not sufficient. pgl0,tl,rdif = ',2f8.4,e12.3,
     2	' utilisation de ETAT_GONG2')
 
	end
