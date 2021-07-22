c
c*****************************************************************************
c
	subroutine clmnew(east,eahat,dmuc,ane,x,y,ea,npar,nitdmu)
c
c  finds Newton iteration correction to ionization fractions
c  Input: east: current value of ionization fractions
c         eahat: Coulomb-independent part of ionization fractions
c         dmuc: Coulomb corrections (dmuc(i,1) contains correction,
c               dmuc(i,2) derivative wrt log rho, dmuc(i, 2+j)
c               derivative wrt ionization fraction no. j.
c         ane: electron number density (per unit mass)
c         x and y: hydrogen and helium abundances
c  Returns updated ionization fractions in ea.
c
c  Original version: 10/5/90
c
      implicit double precision (a-h, o-z)
	implicit integer(i-n)
      parameter(nparw=3, nparw1=nparw+1)
      dimension east(30), eahat(30), ea(30), dmuc(npar,1),
     *  rhld(nparw),w(nparw,nparw1)
      common/constsc/ av,ah,ahe,az,avda,avd1,ck1,ck2,exh,exhe,exhep,
     .  ct,crho,cpe,che,ca03,caa,ckh,car,ergev,exhm
      common/ln10c/ amm,amm2,amm3
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
c
c  set derivatives of log rho wrt ionization fractions
c
      anxfct=x*av/(amm*ane*ah)
      anyfct=y*av/(amm*ane*ahe)
c
c..      write(6,*) 'east:',east(1),east(11),east(21)
c..      write(6,*) 'eahat:',eahat(1),eahat(11),eahat(21)
      denom=1+east(1)
      denom=denom*denom
      rhld(1)=-anxfct/denom
c
      denom=1+east(11)*(1+east(21))
      denom=denom*denom
      rhld(2)=-anyfct*(1+2*east(21))/denom
      rhld(3)=-anyfct*east(11)*(2+east(11))/denom
c..      write(6,*) 'rhld:', rhld
c..      write(6,*) 'dmuc(1,.):',(dmuc(1,j),j=1,5)
c..      write(6,*) 'dmuc(2,.):',(dmuc(2,j),j=1,5)
c..      write(6,*) 'dmuc(3,.):',(dmuc(3,j),j=1,5)
c
c  set up linear equations
c
      do 20 i=1,npar
      is=1+10*(i-1)
      eea=eahat(is)*exp(dmuc(i,1))
      w(i,nparw1)=eea-east(is)
c
      do 15 j=1,npar
   15 w(i,j)=-eea*(dmuc(i,j+2)+dmuc(i,2)*rhld(j))
c
   20 w(i,i)=1+w(i,i)
c
      if(idgeos.ge.4) then
        write(6,110)
        do 22 i=1,npar
   22   write(6,120) (w(i,j),j=1,npar),w(i,nparw1)
      end if
c
c  solve linear equations
c
      call leq(w,w(1,nparw1),npar,1,nparw,nparw,err)
c
c  as pure ad hoc fudge, halve correction if there are convergence
c  problems
c
      if(nitdmu.gt.20) then
        do 24 i=1,npar
   24   w(i,nparw1)=0.5*w(i,nparw1)
      end if
c
c  set updated ea
c
      call storec(east,ea,30)
      do 30 i=1,npar
      is=1+10*(i-1)
   30 ea(is)=ea(is)+w(i,nparw1)
c
      if(idgeos.ge.4) then
        write(6,*) 'Corr. to ea', (w(i,nparw1),i=1,npar)
      end if
c
      return
  110 format(' equations in clmnew:')
  120 format(1p5e13.5)
      end
