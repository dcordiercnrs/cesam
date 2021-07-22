c
c***************************************************************************
c
	subroutine eqstpc(pl, tl, xh, yh, zh, nosd, notd, fl, nit)
c
c	Modification de eqstrh par Annie Baglin Juillet 91
c
c  equation of state routine, with independent variables
c    pl = log(p)
c    tl  = log(T)
c
c  Iterates, using eqstfc, to determine log(f), and hence set equation
c  of state variables.
c  log(f) is returned in fl.
c  nit returns number of iterations.
c
c  Accuracy is set in eps below, and should be adjusted for
c  computer.
c
c  If fl .le. -1.e10 on input, a trial value is set in the routine.
c  Otherwise the value of fl provided is used as trial.
c
c  As for s/r eqstfc, s/r setcnsc must be called before calling
c  eqstrh.
c
c  Note: to obtain the trial value of fl, a fit to mue**-1, as a
c  fonction of log T, is used, in the region of partial ionization.
c  this is currently based on solar composition and
c  log T - log rho relation, and hence may give problems under
c  different circumstances.
c
c  Original version: 07/08/87
c
c  This version is set up for use with Coulomb term version
c  of eqstfc. When starting from trial value of fl, make initial
c  iteration for fl without Coulomb terms.
c
c  Note that present this is the case regardless of input fl
c
c  Version from 1/5/90.
c
c  Double precision version.
c  ++++++++++++++++++++++++
c
c  Dated 9/3/90
c
c
 
c	MODIF : appel a EFF en cas de non convergence P. MOREL oct 1991
c	ajout du common /marche/ avec ok=.true. si c'est bon
 
	implicit real*8 (a-h,o-z)
 
	implicit integer(i-n)
	
	logical nosd, notd
	
      dimension flini(2), wt(2)
      common/eqstdc/ xii1(4), ane(10), rho(20), ht(20), p(20), cp(4),
     *  dad(4), dlt(4), gm1, tprh, trhp, rhxp,gmm1(4)
c
c  common defining standard input and output
c
      common/cstdio/ istdin, istdou, istdpr
      common/diagns/ idgbcs,idgrhs,idgeos,idgopc,idgeng
      common/ccoulm/ epsdmu, icoulm, iclmit, iclmsd, nitdmu, epssdr
 
	logical ok
	common /marche/ok
c
      save
c
      data eps /1.d-10/
c
c  storec original value of Coulomb case flag
c
      icoulp=icoulm
c
c  test for setting trial
c
c
c  try to use wd setting of trial fl in all cases
c
      icase = 1
      icoulm=-1
c	write(6,*) tl,pl
      call inteffc(tl,pl,rhl,flini,wt,icase,iextr)
      fl = wt(1)*flini(1) + wt(2)*flini(2)
c      write(6,*) 'tl, rhl, pl =', tl, rhl,pl
c      write(6,*) 'trial fl set to', fl
      if(fl.le.-1.d10) then
       xt=max(0.d0,tl-3.76119)
        xt=xt*xt
        xmufl=3.8*xt/(xt+0.035)-3.83015
        fl=7.829+xmufl+rhl-1.5*tl
c
      end if
c
c  start loop for iteration
c
      nit=0
   10 call eqstfc(fl, tl, xh, yh, zh, nosd, notd)
c
c
	if(p(1) .le. 0.d0 .or. p(2) .eq. 0.d0) then
         print*,'eqstpc en 105 istdpr,istdou',istdpr,istdou
         write(istdpr, 105) fl, tl, pl, p(1)
105	 format(//' **** error in eqstpc. p .le. 0 or p(2) = 0 for'/
     1	' log f =',f10.5,'  log T =',f10.5,' log p =', f10.5,
     2	' pi =',1pe13.5,' appel a ETAT_EFF')	
	 ok=.false.		!modif appel a gong2
         call dmpeqs
 	 if(.true.)return			!modif appel a gong2
         stop
	end if
	pli=log10(p(1))
	dfl=(pl - pli)/p(2)
	nit=nit+1
c
c  limit size of dfl
c
      if(abs(dfl) .gt. 0.4d0) dfl=sign(0.4d0,dfl)
c      if(idgeos .ge. 1) write(istdpr,*) nit,pl, tl, fl, dfl
c
c  test for continuing iteration
c
      if(abs(dfl) .lt. eps) then
        go to 20
      else if(nit .le. 60) then
        fl = fl+dfl
        go to 10
      else
c
c  diagnostics for failed convergence
c
        write(istdou,110) pl, tl, fl, dfl	!probleme
110	format(//' ***** Iteration failed in s/r eqstrh.'/
     1	7x,'log(p) =',f10.5,'  log(T) =',f10.5,
     2	' Last log(f) =',f10.5,'  dfl =',1pe13.5,' appel a ETAT_EFF')
 
	ok=.false.
	print*,'pb en 110 eqstpc istdou',istdou
        if(istdpr .ne. istdou)then		!modif: appel a GONG2
	 write(istdpr,110) pl, tl, fl, dfl
	 ok=.false.
	endif
	if(.true.)return
        fl = -1.e11
        go to 20
c
      end if
c
c  this is the end
c
   20 continue
c
c  test for repeat with Coulomb term
c
      if(icoulm .ne. icoulp) then
        icoulm=icoulp
        nit=0
        go to 10
      endif
c..      write(6,*) 'converged fl =',fl
 
      return
 
      end
