c---------------------------------------------------------------------
c
c   Personal version of 'Z_interp_IEOS.f' the original OPAL code
c
c---------------------------------------------------------------------
c D. Cordier, Warsaw, April 2003.
c

c    This program interpolates IEOS* files to make a table for the
c    specific Z of interest.  The only input required is the value of
c    Z.  The output file is named EOSdata.  The EOSdata file is
c    read by the interpolation code IEOSxtrin.f. A new EOSdata file
c    is needed everytime Z is changed.

      parameter (nr=170,ntt=200)
      parameter (nruse=169)
      real moles,molesz
      integer long ! D.C.

      character*1 blank, blank1
      character front*3,front1*7,front2*7,front3*7,middle*24,eend(5)*2
      character*24 tabfile,tabfilem,tabfileu,tabfilez
      character*100 file_out ! D.C.
      external long ! D.C.

      dimension z(3),moles(3),tmass(3),frac(3,6),fracz(7)
      dimension density(nr),t6(ntt),alogr(nr,ntt),p(nr,ntt)
     x ,s(nr,ntt),dedrhoa(nr,ntt),cv(nr,ntt),chir(nr,ntt)
     x ,chit(nr,ntt),gam1(nr,ntt),gam2p(nr,ntt),gam3p(nr,ntt)
     x ,icycuse(nr),en(nr,ntt),toff(nr),doff(nr)
      dimension densitym(nr),t6m(ntt),alogrm(nr,ntt),pm(nr,ntt)
     x ,sm(nr,ntt),dedrhoam(nr,ntt),cvm(nr,ntt),chirm(nr,ntt)
     x ,chitm(nr,ntt),gam1m(nr,ntt),gam2pm(nr,ntt),gam3pm(nr,ntt)
     x ,icycusem(nr),enm(nr,ntt),toffm(nr),doffm(nr)
      dimension densityu(nr),t6u(ntt),alogru(nr,ntt),pu(nr,ntt)
     x ,su(nr,ntt),dedrhoau(nr,ntt),cvu(nr,ntt),chiru(nr,ntt)
     x ,chitu(nr,ntt),gam1u(nr,ntt),gam2pu(nr,ntt),gam3pu(nr,ntt)
     x ,icycuseu(nr),enu(nr,ntt),toffu(nr),doffu(nr)


      write(6,*) ' '
      write(6,*) ' ------------------------------------------'
      write(6,*) ' '
      write(6,*) '   Z Interpolation program for OPAL EOS'
      write(6,*) ' '
      write(6,*) '          Personal version D.C.'
      write(6,*) ' '
      write(6,*) ' ------------------------------------------'
      write(6,*) ' '

      write(6,*) ' WARNING: leave all data files gzipped!'
      write(6,*) ' '

      write (*,'("type Z:  ",$)')
      read (*,*) zz
      write (*,'(f15.9)') zz
      write(6,*) ' '
      write(6,*) 'Output file name:'
      read(5,'(A)') file_out
      write(6,*) ' '
      write(6,*) ' Ok! The output file will be called: ',
     & file_out(:long(file_out))

      write(6,*) ' '

      blank1=' '
      front='eos'
      eend(1)='0x'
      eend(2)='2x'
      eend(3)='4x'
      eend(4)='6x'
      eend(5)='8x'
      front1='IEOS00z'
      front2='IEOS02z'
      front3='IEOS04z'
      tabfilez='EOSdata'
c      open (15,file=tabfilez ,status='new' ,err=383) ! D.C.
      open (15,file=file_out(:long(file_out)),status='new' ,err=383)

      do kk=1,5
      tabfile=front1//eend(kk)
      tabfilem=front2//eend(kk)
      tabfileu=front3//eend(kk)
      j=0
          if (kk .eq. 1 ) then
            call system ('gunzip IEOS00z0x')
            call system ('gunzip IEOS02z0x')
            call system ('gunzip IEOS04z0x')
          endif
          if (kk .eq. 2 ) then
            call system ('gunzip IEOS00z2x')
            call system ('gunzip IEOS02z2x')
            call system ('gunzip IEOS04z2x')
          endif
          if (kk .eq. 3 ) then
            call system ('gunzip IEOS00z4x')
            call system ('gunzip IEOS02z4x')
            call system ('gunzip IEOS04z4x')
          endif
          if (kk .eq. 4 ) then
            call system ('gunzip IEOS00z6x')
            call system ('gunzip IEOS02z6x')
            call system ('gunzip IEOS04z6x')
          endif
          if (kk .eq. 5 ) then
            call system ('gunzip IEOS00z8x')
            call system ('gunzip IEOS02z8x')
            call system ('gunzip IEOS04z8x')
          endif
            open (8,file=tabfile ,status='old' ,err=380)
      numprev=0   ! table op

      read (8,'(3x,f6.4,3x,f6.4,11x,f10.7,17x,f10.7)')
     x  x,z(1),moles(1),tmass(1)
      tmassz=gmass(x,z(1),amoles,eground,amassz,fracz)
      read (8,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,
     x 4x,e11.4)') (frac(1,i),i=1,6)
      read (8,'(a)') blank
      do 398 jcs=1,50000
      read (8,'(2i5,2f12.7,17x,e15.7)') numtot,icycuse(jcs)
     x  ,toff(jcs),doff(jcs),density(jcs)
      if(numtot .eq. 0)then
       backspace 8
      go to 397
      endif
      numprev=numprev+1   ! table op
      read(8,'(a)') blank
      read(8,'(a)') blank
      do 399 i=1,icycuse(jcs)
      read (8,'(f10.5,f7.2,3e14.6,6f8.4)')t6(i),alogr(jcs,i)
     x ,p(jcs,i),en(jcs,i),s(jcs,i),cv(jcs,i),chir(jcs,i)
     x ,chit(jcs,i),gam1(jcs,i),gam2p(jcs,i),gam3p(jcs,i)
      if (im1 .lt. 200) then
        im1=im1+1
      endif
      p0=t6(i)*density(jcs)
      p(jcs,i)=p(jcs,i)/p0
      en(jcs,i)=(en(jcs,i)*tmass(1)-eground)/(t6(i)*tmass(1))
  399 continue
      read(8,'(a)') blank
      read(8,'(a)') blank
      read(8,'(a)') blank
  398 continue
  397 continue

            open (9,file=tabfilem ,status='old' ,err=381)
      numprevm=0   ! table op

      read (9,'(3x,f6.4,3x,f6.4,11x,f10.7,17x,f10.7)')
     x  x,z(2),moles(2),tmass(2)
      tmassz=gmass (x,z(2),amoles,eground,amassz,fracz)
      read (9,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,
     x 4x,e11.4)') (frac(2,i),i=1,6)
      read (9,'(a)') blank
      do 598 jcs=1,50000
      read (9,'(2i5,2f12.7,17x,e15.7)') numtotm,icycusem(jcs)
     x  ,toffm(jcs),doffm(jcs),densitym(jcs)
      if(numtotm .eq. 0)then
       backspace 9
      go to 597
      endif
      numprevm=numprevm+1   ! table op
      read(9,'(a)') blank
      read(9,'(a)') blank
      do 599 i=1,icycusem(jcs)
      read (9,'(f10.5,f7.2,3e14.6,6f8.4)')t6m(i),alogrm(jcs,i)
     x ,pm(jcs,i),enm(jcs,i),sm(jcs,i),cvm(jcs,i)
     x ,chirm(jcs,i)
     x ,chitm(jcs,i),gam1m(jcs,i),gam2pm(jcs,i),gam3pm(jcs,i)

      p0=t6(i)*density(jcs)
      p0=t6m(i)*densitym(jcs)
      pm(jcs,i)=pm(jcs,i)/p0
      enm(jcs,i)=(enm(jcs,i)*tmass(2)-eground)
     x /(t6m(i)*tmass(2))
  599 continue
      read(9,'(a)') blank
      read(9,'(a)') blank
      read(9,'(a)') blank
  598 continue
  597 continue


            open (10,file=tabfileu ,status='old' ,err=382)
      numprevu=0   ! table op

      read (10,'(3x,f6.4,3x,f6.4,11x,f10.7,17x,f10.7)')
     x  x,z(3),moles(3),tmass(3)
      tmassz=gmass (x,z(3),amoles,eground,amassz,fracz)
      read (10,'(21x,e14.7,4x,e14.7,3x,e11.4,3x,e11.4,3x,e11.4,
     x 4x,e11.4)') (frac(3,i),i=1,6)
      read (10,'(a)') blank
      do 798 jcs=1,50000
      read (10,'(2i5,2f12.7,17x,e15.7)') numtotu,icycuseu(jcs)
     x  ,toffu(jcs),doffu(jcs),densityu(jcs)
      if(numtotu .eq. 0)then
       backspace 10
      go to 797
      endif
      numprevu=numprevu+1   ! table op
      read(10,'(a)') blank
      read(10,'(a)') blank
      do 799 i=1,icycuseu(jcs)
      read (10,'(f10.5,f7.2,3e14.6,6f8.4)')t6u(i),alogru(jcs,i)
     x ,pu(jcs,i),enu(jcs,i),su(jcs,i),cvu(jcs,i)
     x ,chiru(jcs,i)
     x ,chitu(jcs,i),gam1u(jcs,i),gam2pu(jcs,i),gam3pu(jcs,i)
      p0=t6(i)*density(jcs)
      p0=t6u(i)*densityu(jcs)
      pu(jcs,i)=pu(jcs,i)/p0
      enu(jcs,i)=(enu(jcs,i)*tmass(3)-eground)
     x /(t6u(i)*tmass(3))
  799 continue
      read(10,'(a)') blank
      read(10,'(a)') blank
      read(10,'(a)') blank
  798 continue
  797 continue

      do jcs=1,numprev
      do i=1,icycuse(jcs)
      alogr(jcs,i)=quad(zz,alogr(jcs,i),alogrm(jcs,i),alogru(jcs,i)
     x ,z(1),z(2),z(3))
      p(jcs,i)=quad(zz,p(jcs,i),pm(jcs,i),pu(jcs,i)
     x ,z(1),z(2),z(3))
      en(jcs,i)=quad(zz,en(jcs,i),enm(jcs,i),enu(jcs,i)
     x ,z(1),z(2),z(3))
      s(jcs,i)=quad(zz,s(jcs,i),sm(jcs,i),su(jcs,i)
     x ,z(1),z(2),z(3))
      cv(jcs,i)=quad(zz,cv(jcs,i),cvm(jcs,i),cvu(jcs,i)
     x ,z(1),z(2),z(3))
      chir(jcs,i)=quad(zz,chir(jcs,i),chirm(jcs,i),chiru(jcs,i)
     x ,z(1),z(2),z(3))
      chit(jcs,i)=quad(zz,chit(jcs,i),chitm(jcs,i),chitu(jcs,i)
     x ,z(1),z(2),z(3))
      gam1(jcs,i)=quad(zz,gam1(jcs,i),gam1m(jcs,i),gam1u(jcs,i)
     x ,z(1),z(2),z(3))
      gam2p(jcs,i)=quad(zz,gam2p(jcs,i),gam2pm(jcs,i),gam2pu(jcs,i)
     x ,z(1),z(2),z(3))
      gam3p(jcs,i)=quad(zz,gam3p(jcs,i),gam3pm(jcs,i),gam3pu(jcs,i)
     x ,z(1),z(2),z(3))
      enddo
      enddo

      tmassz=gmass (x,zz,amoles,eground,amassz,fracz)

      write (15,'(" X=",f6.4," Z=",f12.9,5x,"moles=",f10.7,
     x "  mass/mole(ion)=",f10.7,3x,"eground=",f12.7)') 
     x  x,zz,amoles,tmassz,eground
      write (15,'(" Number fractions: ","H=",1pe14.7," He=",e14.7,
     x " C=",e11.4," N=",e11.4," O=",e11.4," Ne=",e11.4)')
     x (fracz(i),i=6,1,-1)
      write(15, '(a)') blank

      kcs=0
      numuse=numprev
      if (numuse .gt. nruse) numuse=nruse
      do 498 jcs=1,numuse
      icycuse(jcs)=min0( icycuse(jcs), icycusem(jcs), icycuseu(jcs) )
      kcs=kcs+1
      write(15,'(2i5,2f12.7,"  density(g/cc)= ",1pe15.7,
     x" U=10e+12ergs/gm Us=U/(E+6K) P0=T6*density W=sqrt(P0)")')
     x  kcs,icycuse(jcs),toff(jcs),doff(jcs), density(jcs)
      write(15,
     x        '("  T6      log R",
     x "   P/P0         E(U)/T6      S(Us)         Cv   "
     x    ,"   chir  ","  chit  "
     x       ,"  g1    "," g2/(g2-1)", " g3-1",/)')
c     write (15,'(1h0)')


      do 499 i=1,icycuse(jcs)
c     read(9,'(1h0)')
      write (15,'(f9.5,1x,f6.2,1p3e13.5,0p6f8.4)')t6(i),
     x alogr(jcs,i),p(jcs,i),en(jcs,i),s(jcs,i),cv(jcs,i)
     x,chir(jcs,i),chit(jcs,i),gam1(jcs,i),gam2p(jcs,i),gam3p(jcs,i)
  499 continue
      write(15,'(a)') blank
      write(15,'(a)') blank
      write(15,'(a)') blank

  498 continue
      izero=0
      write (15,'(2i5)') izero,izero
          if (kk .eq. 1 ) then
            call system ('gzip IEOS00z0x')
            call system ('gzip IEOS02z0x')
            call system ('gzip IEOS04z0x')
          endif
          if (kk .eq. 2 ) then
            call system ('gzip IEOS00z2x')
            call system ('gzip IEOS02z2x')
            call system ('gzip IEOS04z2x')
          endif
          if (kk .eq. 3 ) then
            call system ('gzip IEOS00z4x')
            call system ('gzip IEOS02z4x')
            call system ('gzip IEOS04z4x')
          endif
          if (kk .eq. 4 ) then
            call system ('gzip IEOS00z6x')
            call system ('gzip IEOS02z6x')
            call system ('gzip IEOS04z6x')
          endif
          if (kk .eq. 5 ) then
            call system ('gzip IEOS00z8x')
            call system ('gzip IEOS02z8x')
            call system ('gzip IEOS04z8x')
          endif
      enddo

      close (15)
c      call system ('gzip EOSdata') ! D.C.
      stop
  380 write (*,'("File missing is ",a)') tabfile
      stop
  381 write (*,'("File missing is ",a)') tabfilem
      stop
  382 write (*,'("File missing is ",a)') tabfileu
      stop
  383 write (*,'("File already exists:  ",a)') tabfilez
      stop
      end

c************************************************************************
      function quad(x,y1,y2,y3,x1,x2,x3)
c..... this function performs a quadratic interpolation.
      save
      dimension  xx(3),yy(3)
      xx(1)=x1
      xx(2)=x2
      xx(3)=x3
      yy(1)=y1
      yy(2)=y2
      yy(3)=y3
        if(ic .eq. 0) then
          ic=1
          xx12=1./(xx(1)-xx(2))
          xx13=1./(xx(1)-xx(3))
          xx23=1./(xx(2)-xx(3))
          xx1sq=xx(1)*xx(1)
          xx1pxx2=xx(1)+xx(2)
        endif
      c3=(yy(1)-yy(2))*xx12
      c3=c3-(yy(2)-yy(3))*xx23
      c3=c3*xx13
      c2=(yy(1)-yy(2))*xx12-xx1pxx2*c3
      c1=yy(1)-xx(1)*c2-xx1sq*c3
      dkap=c2+(x+x)*c3
      quad=c1+x*(c2+x*c3)
      return
      end
c******************************************************************
      function gmass(x,z,amoles,eground,fracz,frac)
      dimension anum(6),frac(7), amas(7),eion(7)
      data (eion(i),i=1,6)/-3394.873554,-1974.86545,-1433.92718,
     x  -993.326315,-76.1959403,-15.28698437/
      data (anum(i),i=1,6)/10.,8.,7.,6.,2.,1./
      xc=0.247137766
      xn=.0620782
      xo=.52837118
      xne=.1624188
      amas(7)=1.0079
      amas(6)=4.0026
      amas(5)=12.011
      amas(4)=14.0067
      amas(3)=15.9994
      amas(2)=20.179
      amas(1)=0.00054858
      fracz=z/(xc*amas(5)+xn*amas(4)+xo*amas(3)+xne*amas(2))
      xc2=fracz*xc
      xn2=fracz*xn
      xo2=fracz*xo
      xne2=fracz*xne 
      xh=x/amas(7)
      xhe=(1.-x -z)/amas(6)
      xtot=xh+xhe+xc2+xn2+xo2+xne2
      frac(6)=xh/xtot
      frac(5)=xhe/xtot
      frac(4)=xc2/xtot
      frac(3)=xn2/xtot
      frac(2)=xo2/xtot
      frac(1)=xne2/xtot
      eground=0.0
      amoles=0.0
      do i=1,6
      eground=eground+eion(i)*frac(i)  
      amoles=amoles+(1.+anum(i))*frac(i)
      enddo
      anume=amoles-1.
      gmass=0.
      do i=2,7
      gmass=gmass+amas(i)*frac(i-1)
      enddo

      return
      end
