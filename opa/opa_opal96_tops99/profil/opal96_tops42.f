c--------------------------------------------------------------------------
  
      subroutine opal96_tops42(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)


c Routine d'opacite avec opacite TOPS version 4.0, uniquement pour
c etoiles avec Z0 = 0.02, 0.008 et 0.004

c Daniel, Mai 1998

c
c      calcul de l'opacite repris du code de Geneve
c      polynomes de Lagrange pour log T6 et log R
c      quadratique pour X et Z
c       Yveline (aout 1991-code de Geneve--> 10/95 CESAM)   
c

c Remarques : * l'opacite conductive est incluse!!!
c             * les routines supplementaires ont ete mises a la
c               fin du fichier.

c-------------------------------------------------------------------
c      entree :
c      xh(1)=X : comp. chim. en H
c      t : temperature K
c      ro : densite cgs
  
c      sortie :
c      kappa : opacite gr / cm2)
c      dkapdt : kappa / d t
c      dkapdr : kappa / d densite              
c      dkapdx : kappa / d xchim(1)
  
c      Z est obtenu par 1-X-Y
c-------------------------------------------------------------------
  
      implicit none

      include 'cesam_3.parametres'
      include 'modele_3.common'
      include 'evol_chim_3.common'

      logical tops

      real*8 t6, lr, xh(1), t, ro, kappa, dkapdt, dkapdr,
     +       dkapdx, xi, zi, y, dc12, do16, xc12_0, xo16_0,
     +       log_do16sdc12, ymax, delta, bid1, bid2, bid3,
     +       kap_opal, kap_tops

      character*50 file_opa_tops

c     ------ calcul des valeurs pour lesquelles on veut l'opacite ---
  
      t6=t*1.d-6			!T6
      lr=log10(ro/t6**3)		!log R
      xi=xh(1)		!determination de X et de Z

      y=xh(ihe4)+xh(ihe4-1)

      zi=1.d0-xi-xh(ihe4)-xh(ihe4-1)

      if( z0 .eq. 0.02d0 ) then
        xc12_0=0.34664719d-02
        xo16_0=0.96438571d-02
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa2
      end if

      if( z0 .eq. 0.008d0 ) then
        xc12_0=0.13865888d-02
        xo16_0=0.38575428d-02
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa3
      end if

      if( z0 .eq. 0.004d0 ) then
        xc12_0=0.69329439d-03
        xo16_0=0.19287714d-02
        dc12=xh(4)-xc12_0
        do16=xh(9)-xo16_0
        file_opa_tops=opa4
      end if

      if( dc12 .gt. 0.d0 .AND. do16 .gt. 0.d0 ) then
          log_do16sdc12=log10(do16/dc12)
          tops=.true.
      else
          tops=.false.
      end if


      if( tops ) then
         xi=0.d0
         call kappa_opal_3(lr,t6,ro,xi,z0,kap_opal,bid1,bid2,
     +                     bid3,opa1)
         ymax=1.d0-z0
         call kappa_tops_4(ymax,log_do16sdc12,t,ro,file_opa_tops,
     +           kap_tops,dkapdt,dkapdr)
         delta=kap_tops-kap_opal
         call kappa_tops_4(y,log_do16sdc12,t,ro,file_opa_tops,
     +           kappa,dkapdt,dkapdr)
         kappa=kappa-delta
         if( kappa .le. 0.d0 ) then
             print*, 'Pb. avec kappa .le. 0. !'
             stop
         end if
         dkapdx=0.d0
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      else
         call kappa_opal_3(lr,t6,ro,xi,zi,kappa,dkapdr,dkapdt,
     +                     dkapdx,opa1)
         call kappa_cond_3(xh,t,ro,kappa,dkapdt,dkapdr,dkapdx)
         return
      end if

      end
  
cc------------------------------------------------------------------------
c
      subroutine kappa_opal_3(lr,t6,ro,x,z,cap,capr,capt6,capx,opa1)
c
c----------------------------------------------------------------------
c
c     sous-programme de calcul de l'opacite radiative par interpolation
c     dans les tables de Livermore. Ces tables donnent le log10 de
c     l'opacite pour differentes valeurs de z, x, T6(T/1e6) et
c     log R (R=ro/T6^3). Les valeurs tabulees de l'opacite correspondent
c     a l'opacite radiative on entre en parametre nz, nx, nt, nr qui
c     sont respectivement le nombre de valeurs de z, de x, de T et de R.
c
c     parametres d'entree de kappa_opal:
c     t6 et log R     
c     x et z: H et metaux dans la couche consideree
c     opa1: nom de la table d'opacite
c     sortie: kappa, dlogkappa/dt6,dlogkappa/dlogR, dlogkappa/dx
c     
c     
c              Yveline (aout 1991-code de Geneve--> 10/95 CESAM)
c
c              Daniel (octobre 1996-version pour CESAM 2.x --> CESAM 3.1)
c
c----------------------------------------------------------------------
c
      implicit none

      integer nz,nx,nt,nr,pnz,pnx,pnt,pnr
      parameter(pnz=16,pnx=8,pnt=85,pnr=23)
c
      real*8 vz,vx,vlt6,vlr,vlk,lr,T6,lt6,ro
      common/val_opal/vz(pnz),vx(pnx),vlt6(pnt),vlr(pnr),
     &                vlk(pnz,pnx,pnt,pnr),nt,nr,nz,nx
c
      integer iz,iza,izb,ix,ixa,ixb,it,ita,ir,ira,ilag,ilin,i,j
      real*8 opa,opar,opat6,opax,cap,capr,capx,capt6,x,z,
     &       kapq,dkapq,qx,qk,qkk,qkl,opaxx,oparx,opat6x,opatx
c
      dimension ilag(3,3),ilin(3,3),opa(3,3),opar(3,3),opat6(3,3),
     &          opax(3),oparx(3),opatx(3),opaxx(3),qx(3),qk(3),
     &          qkk(3),qkl(3),opat6x(3)

      
      character*50 opa1

      logical i_lect,i_exit
      data i_lect/.true./
      data i_exit/.false./
  
c
c
c     ----- lecture des tables d'opacite ------------------------------
      if(i_lect) then
        i_lect=.false.
        call lect_opal(vz,vx,vlt6,vlr,vlk,opa1,nt,nr,nx,nz)
      endif
c
c
c     ----- localisation dans les tables du point pour lequel on cherche
c     l'opacite. Sortie eventuelle des tables. -------------------------
c
c     ---------- en composition chimique -------------------------------
  
      
      if(z.lt.vz(1)) then
        z=vz(1)
        write(6,10) x,z
      endif
  
      if(z.gt.vz(nz)) then
        z=vz(nz)
        write(6,10) x,z
      endif
  
      if(x.lt.vx(1)) then
        x=vx(1)
        write(6,10) x,z
      endif
  
      if(x.gt.vx(nx)) then
        x=vx(nx)
        write(6,10) x,z
      endif
      
10      format(/1X,' attempt to extrapolate beyond the domain of
     &  validity of the opacity tables , X=',F7.4,' Z=',F7.4)
  
      call pos_table_op(z,vz,nz,iz)
      call pos_table_op(x,vx,nx,ix)
c
c     ---------- en densite et temperature ----------------------------
  
      lt6=log10(t6)
  
      if(lr.lt.vlr(1)) then
        write(6,20) lr,T6
        lr=vlr(1)
      endif
  
      if(lr.gt.vlr(nr)) then
        write(6,20) lr,T6
        lr=vlr(nr)
      endif
  
      if(lT6.lt.vlt6(1)) then
        write(6,20) lr,T6
        lT6=vlt6(1)
      endif
  
      if(lT6.gt.vlt6(nt)) then
        write(6,20) lr,T6
        lT6=vlt6(nt)
      endif
  
20    format(/1X,' attempt to extrapolate beyond the domain of
     &  validity of the opacity tables , R_log=',F7.4,' T6=',F8.4)
  
      call pos_table_op(lr,vlr,nr,ir)
      call pos_table_op(lT6,vlt6,nt,it)
c
c
c     -----------------------------------------------------------------
c     In this version of the opacity tables, a 4-point double argument
c     interpolation is incorporated. This requires 16 given opacity
c     entries, and the input pair (rh,T) should be preferently in the
c     center square. The ideal configuration required is the following
c     (horizontal axis:log R, vertical axis: T6, *: entry given in the
c     table, X: position of the input (log R,T6)
c
c        *   *   *   *
c
c        *   *   *   *
c              x
c        *   *   *   *
c
c        *   *   *   *
c
c      the 16-point configuration is only used if all corresponding
c      entries are provided. If a 16-point configuration is impossible
c      a linear interpolation is made if the 4 surrounding points have
c      their opacity entries, otherwise error exit occurs. In the case
c      of very low densities, however, extrapolation is done by
c      constant continuation. In this way spurious increases of the
c      opacity due to 3-rd degree extrapolation are suppressed.
c     -----------------------------------------------------------------
c
c     ----- Peut-on faire une interpolation Lagrangienne ou seulement
c     une interpolation lineaire? -------------------------------------
c
c     ---------- test pour l'interpolation lagrangienne ---------------
c
      if(ir.ge.3.and.ir.le.(nr-1).and.
     &   it.ge.3.and.it.le.(nt-1))     then
        izb=0
        do 1 iza=iz-1,iz+1
          izb=izb+1
          ixb=0
          do 2 ixa=ix-1,ix+1
            ixb=ixb+1
            ilag(izb,ixb)=1
            do 3 ita=it-2,it+1
              do 4 ira=ir-2,ir+1
c                if(vlk(iza,ixa,ita,ira).eq.-100.) then
                if(vlk(iza,ixa,ita,ira).gt.98.) then
                  ilag(izb,ixb)=0
                  goto 100
                endif
4             continue
3           continue
100          continue
2         continue
1       continue
      else
        do 5 i=1,3
          do 6 j=1,3
            ilag(i,j)=0
6         continue
5       continue
      endif
c     ---------- sinon test pour l'interpolation lineaire -------------
c
      izb=0
      do 7 iza=iz-1,iz+1
        izb=izb+1
        ixb=0
        do 8 ixa=ix-1,ix+1
          ixb=ixb+1
          if(ilag(izb,ixb).eq.0) then
102         continue
            do 9 ita=it-1,it
101           continue
              do 11 ira=ir-1,ir
c                if(vlk(iza,ixa,ita,ira).eq.-100.) then
                if(vlk(iza,ixa,ita,ira).gt.98.) then
                  i_exit=.true.
                  if(ira.eq.ir.and.abs(lr-vlr(ir-1)).lt.1e-15.
c     &               and.vlk(iza,ixa,ita,ir-1).ne.-100.) then
     &               and.vlk(iza,ixa,ita,ir-1).le.98.) then
                    ir=ir-1
                    i_exit=.false.
                    goto 101
                  endif
                endif
11             continue
              if(i_exit) then
                if(it.ne.2.and.ita.eq.it.and.
     &             abs(lT6-vlt6(it-1)).lt.1e-15.and.
c     &             vlk(iza,ixa,it-1,ira).ne.-100.) then
     &             vlk(iza,ixa,it-1,ira).le.98.) then
                  it=it-1
                  i_exit=.false.
                  goto 102
                else
                  write(6,30) lr,T6
30                format(/1X,'attempt to extrapolate beyond the domain
     &                   of validity of opacity tables, vlk=99.,
     &                   R_log=',F7.4,' T6_log=',F8.4)
                   stop
c                  call exit
                endif
              endif
9           continue
            ilin(izb,ixb)=1
          endif
8       continue
7     continue
c
c
c     ----- calcul des opacites pour les points voisins ---------------
c
      izb=0
      do 12 iza=iz-1,iz+1
        izb=izb+1
        ixb=0
        do 13 ixa=ix-1,ix+1
          ixb=ixb+1
          if(ilag(izb,ixb).eq.1) then
            call intl_opal(iza,ixa,ir,it,lr,lT6,opa(izb,ixb),
     &                     opar(izb,ixb),opat6(izb,ixb))
          else
            if(ilin(izb,ixb).eq.1) then
              call intlin_opal(vlr(ir-1),vlr(ir),
     &                         vlt6(it-1),vlt6(it),
     &                         vlk(iza,ixa,it-1,ir-1),
     &                         vlk(iza,ixa,it-1,ir),
     &                         vlk(iza,ixa,it,ir-1),vlk(iza,ixa,it,ir),
     &                         lr,lT6,opa(izb,ixb),
     &                         opar(izb,ixb),opat6(izb,ixb))
            else
              print*,'pb opacites: ilin=ilag=0'
              stop
c              call exit
            endif
          endif
          opa(izb,ixb)=10.**opa(izb,ixb)
          opat6(izb,ixb)=opa(izb,ixb)/t6/1d6*(-3.*opar(izb,ixb)+
     &                    opat6(izb,ixb))
          opar(izb,ixb)=opa(izb,ixb)/ro*opar(izb,ixb)
13      continue
12    continue
c
c     -----------------------------------------------------------------
c     interpolation quadratique (K, X)
c     -----------------------------------------------------------------
  
      do 16 i=1,3
        qx(i)=vx(ix+i-2)
16    continue
  
      do 17 i=1,3
        do 70 j=1,3
          qk(j)=opa(i,j)
          qkk(j)=opar(i,j)
          qkl(j)=opat6(i,j)
70      continue
        call quad(x,qx,qk,kapq,dkapq)
        opax(i)=kapq
        opaxx(i)=dkapq
        call quad(x,qx,qkk,kapq,dkapq)
        oparx(i)=kapq
        call quad(x,qx,qkl,kapq,dkapq)
        opat6x(i)=kapq
17    continue
  
c
c     -----------------------------------------------------------------
c     interpolation quadratique (K, Z)
c     -----------------------------------------------------------------
  
      do 19 i=1,3
        qx(i)=vz(iz+i-2)
19    continue
  
      call quad(z,qx,opax,cap,dkapq)
      call quad(z,qx,opaxx,capx,dkapq)
      call quad(z,qx,oparx,capr,dkapq)
      call quad(z,qx,opat6x,capt6,dkapq)
  
      return
      end
  
c------------------------------------------------------------------------
cc------------------------------------------------------------------------
      subroutine intlin_opal(R1,R2,T1,T2,vk11,vk12,vk21,vk22,lr,
     &                       T6,opac,opacr,opact6)
c
      implicit none
c
      real*8 R1,R2,T1,T2,vk11,vk12,vk21,vk22,lr,T6,opac,opacr,
     &     opact6,dr,dt,x,y,dr1,dr2,dt1,ct1,ct2
c
      dr  = R2 - R1
      dt = T2 - T1
      x = lr - R1
      y = T6 - T1
      dr1 = (vk12 - vk11)/dr
      dr2 = (vk22 - vk21)/dr
      ct1 = vk11 + x*dr1
      ct2 = vk21 + x*dr2
      dt1 = (ct2 - ct1)/dt
      opac = ct1 + y*dt1
      opacr = dr1+ (dr2-dr1)*(y/dt)
      opact6 = dt1
      return
      end
c------------------------------------------------------------------------
cc------------------------------------------------------------------------
c
      subroutine intl_opal(iza,ixa,ir,it,lr,T6,opac,opacr,opact6)
c
c---- lagrangian interpolation of opac in (log R,T6) within a table
c
c
      implicit none
c
      integer pnz,pnx,pnt,pnr
      parameter(pnz=16,pnx=8,pnt=85,pnr=23)
c
      real*8 vz,vx,vlt6,vlr,vlk
      integer nz,nx,nt,nr
  
      common/val_opal/vz(pnz),vx(pnx),vlt6(pnt),vlr(pnr),
     &                vlk(pnz,pnx,pnt,pnr),nt,nr,nz,nx
  
      integer iza,ixa,i,ir,it,k
      real*8 xr,xt,lr,T6,pr,ppr,pt,ppt,s1,s2,
     &     opac,opacr,opact6
      dimension xr(4),xt(4),s1(4),s2(4),pr(4),ppr(4),pt(4),ppt(4)
c
      do 1 i=1,4
        xr(i)=vlr(ir-3+i)
        xt(i)=vlt6(it-3+i)
1     continue
c
      call lpol_op(xr,4,lr,pr,ppr)
      call lpol_op(xt,4,T6,pt,ppt)
c
      do 2 k=1,4
        s1(k)=0.
        s2(k)=0.
c
        do 3 i=1,4
          s1(k)= s1(k)+pr(i)*vlk(iza,ixa,it-3+k,ir-3+i)
          s2(k)= s2(k)+ppr(i)*vlk(iza,ixa,it-3+k,ir-3+i)
3       continue
2     continue
c
      opac=0.
      opacr=0.
      opact6=0.
      do 4 k=1,4
        opac=opac+pt(k)*s1(k)
        opacr=opacr+pt(k)*s2(k)
        opact6=opact6+ppt(k)*s1(k)
4     continue
c
c
      return
      end
c------------------------------------------------------------------------
cc------------------------------------------------------------------------
      subroutine quad(x,qx,qk,k,dk)
  
c     interpolation quadratique
  
      implicit none
c
      real*8 x,qx,qk,k,dk,x1,x2,x3,y1,y2,y3,denom,aa,bb,cc
  
      dimension qx(3),qk(3)
c
      x1=qx(1)                                         
      x2=qx(2)                                           
      x3=qx(3)
      y1=qk(1)
      y2=qk(2)
      y3=qk(3)                                         
      denom=(x2-x1)*(x3-x1)*(x2-x3)                       
      aa=(y2-y1)*(x3-x1)-(y3-y1)*(x2-x1)                
      aa=aa/denom                                       
      bb=(y3-y1)/(x3-x1)-aa*(x3+x1)
      cc=y1-aa*x1*x1-bb*x1                              
      k=aa*x*x+bb*x+cc                         
      dk=2.*aa*x+bb                            
  
      return
      end
cc------------------------------------------------------------------------
cc------------------------------------------------------------------------
c
      subroutine lect_opal(vz,vx,vlt6,vlr,vlk,opa1,nt,nr,nx,nz)
c
c----------------------------------------------------------------------
c     sous-programme de lecture des tables d'opacite de Livermore
c     les sorties sont vz=z, vx=x, vlt6=log T6, vlrro=log R et
c     vlk=log kappa (log decimaux). on affecte la valeur 99.
c     aux points pour lesquels l'opacite n'est pas donnee.
c----------------------------------------------------------------------
c
      implicit none
c
      integer pnz,pnx,pnt,pnr
      parameter(pnz=16,pnx=8,pnt=85,pnr=23)
c
      integer iunit_lect,jz,jx,jt,jr,nz,nx,nt,nr
      real*8 vz(pnz),vx(pnx),vlt6(pnt),vlr(pnr),
     &       vlk(pnz,pnx,pnt,pnr)
  
      character*50 opa1
c
c
c     ----- ouverture des fichiers d'opacite --------------------------
      iunit_lect=11      
      open(unit=iunit_lect,form='unformatted',file=opa1)      
c
c     ----- lecture des tables ----------------------------------------
      read(iunit_lect) nz,nx,nt,nr
  
      write(6,*) 'tables d''opacite:'
      write(6,*) nz,'valeurs de Z'
      write(6,*) nx,'valeurs de X'
      write(6,*) nt,'valeurs de T6'
      write(6,*) nr,'valeurs de log R'
  
      read(iunit_lect) (vlt6(jt),jt=1,nt)
      read(iunit_lect) (vlr(jr),jr=1,nr)
      do 1 jz=1,nz
        read(iunit_lect) vz(jz)
        do 2 jx=1,nx
          read(iunit_lect) vx(jx)
          read(iunit_lect) ((vlk(jz,jx,jt,jr),jr=1,nr),jt=1,nt)
2       continue
1     continue
  
      do 3 jt=1,nt
        vlt6(jt)=log10(vlt6(jt))
3     continue
        
c
c     ----- fermeture des fichiers d'opacite --------------------------
      close(iunit_lect)
c
      return
      end
  
  
c------------------------------------------------------------------------
cc------------------------------------------------------------------------
c
      subroutine lpol_op(x,n,x0,p,pp)
c
c     ------------------------------------------------------------------
c     polynomes de Lagrange (4 points) et leurs derivees (n=4)
c     x=coordonnees tabulees
c     x0=coordonnee pour laquelle p et pp sont demandes
c     p=valeur des polynomes en x0
c     pp=valeur des derivees de ces polynomes en x0
c     ------------------------------------------------------------------
c
      implicit none
c
      real*8 x,x0,p,pp,s
      integer i,j,k,l,n
c
      dimension x(4),p(4),pp(4),s(4)
c
      do 1 j=1,n
        p(j)=1.
        pp(j)=0.
        s(j)=1.
        do 2 i=1,n
          if(i.ne.j) then
            p(j)=p(j)*(x0-x(i))/(x(j)-x(i))
            s(j)=s(j)*(x(j)-x(i))
          endif
2       continue
        do 3 k=1,n
        if(k.ne.j) then
          do 4 l=1,n
          if (l.ne.j) then
            if (l.gt.k) then
            pp(j)=pp(j)+(x0-x(k))*(x0-x(l))
            endif
          endif
4         continue
        endif
3       continue
      pp(j)=pp(j)/s(j)
1     continue
      return
      end
c------------------------------------------------------------------------
cc------------------------------------------------------------------------
      subroutine pos_table_op(val,tab_val,nval,ival)
c
c     reperage dans une table
  
      implicit none
c
      integer i,ival,nval
      real*8 val,tab_val,t1,t2
      dimension tab_val(100)
c
      ival=0
      if((tab_val(2)).ge.val) then
        ival=2
        return
      endif
      do 1 i=3,nval
        t1=(tab_val(i-2)+tab_val(i-1))/2.
        t2=(tab_val(i-1)+tab_val(i))/2.
        if(t1.le.val.and.t2.gt.val) then
          ival=i-1
          return
        endif
1     continue
  
      ival=nval-1
      return
      end

c************************************************************************

      subroutine kappa_tops_4( y, log_do16sdc12, t, ro, file_opa,
     +                         kap, dkapt, dkapro)

c	Sous-programme de calcul de l'opacite radiative par interpolation
c	dans les tables fabriquees a l'aide du serveur :

c	http://t4.lanl.gov/opacity/tops.html

c	avec :		Z0 = 0.004, 0.008 et 0.02

c	constantes :	X = 0., log10( DNe22 / DC12 ) = -1.4

c	Melange Grevesse, Noels 1993 ( pour Z0 seulement !!)

c	ATTENTION : interpolation parabolique !!!

c	Par hypothese on a ici : Y+Z=1.


c Entrees : y               : fraction massique en He4
c           log_do16sdc12   : log10( D016 / DC12 ) 
c           t               : la temperature (en K)
c           ro              : la masse volumique (g.cm^-3)
c           file_opa        : le nom du fichier binaire ou sont les tables

c Sorties : kap	    : l'opacite (cm^2.g^-1)
c           dkapdt  : derivee par rapport a t
c           dkapdro : derivee par rapport a ro

c	Version :	4.0
c	Auteur	:	Daniel Cordier
c			22 Mai 1998


      implicit none

      integer n_y, n_o16, nt, nro
      integer jt, jr, j_y, j_o16
      integer pn_y, pn_o16, pnt, pnro, max 
     
      parameter( pn_y=4, pn_o16=3, pnt=26, pnro=30)

      parameter ( max = 30 )
c Remarque : max = max( pn_y, pn_o16, pnt, pnro ) !!!

      logical i_lect

      character*50 file_opa

      real*8 y, log_do16sdc12, z0

      real*8 vlt6(pnt), vlro(pnro), ro, lro, t6, lt6,
     + vlk(pn_y,pn_o16,pnt,pnro), lkappa, dkappadt, dkappadro,t,
     + kap, dkapt, dkapro, val_y(pn_y), val_o16(pn_o16)

      real*8 fro(pnro), kro(pn_y,pn_o16,pnt), dkro(pn_y,pn_o16,pnt),
     +       ft(max), krot(pn_y,pn_o16), dkrot(pn_y,pn_o16),
     +       dkrob(pn_y,pn_o16), bidon

      real*8 kroto16(pn_y),dkrobo16(pn_y),dkroto16(pn_y)

      save i_lect, vlt6, vlro, vlk, val_y, val_o16

      data i_lect /.true./

      lro=log10(ro)
      t6=t/10.**6
      lt6=log10(t6)

c Lecture des tables dans le fichier 'file_opa'

      if(i_lect) then
        call lect_tops_4(file_opa,n_y,n_o16,nt,nro,z0,val_y,val_o16,
     +                   vlt6,vlro,vlk)
c      print*, 'vz= ', vz
c      print*, 'vlt6= ', vlt6
c      print*, 'vlro= ', vlro
        i_lect=.false.
      end if

c Tests d'appartenance de ro, t aux differents domaines de validite.
c densite

      if(lro .lt. vlro(1)) then
        print*, 'Pb. avec opacite TOPS, ro= ', ro,' g.cm^-3'
        stop
      end if

      if(lro .gt. vlro(nro)) then
        print*, 'Pb. avec opacite TOPS, ro= ', ro,' g.cm^-3'
        stop
      end if

c temperature

      if(lt6 .lt. vlt6(1)) then
        print*, 'Pb. avec opacite TOPS, t= ', t, ' K'
        stop
      end if

      if(lt6 .gt. vlt6(nt)) then
        print*, 'Pb. avec opacite TOPS, t= ', t, ' K'
        stop
      end if

c      print*, 'OK 1'

c Valeurs de log10( DO16 / DC12 )

      if( log_do16sdc12 .lt. val_o16(1) ) then
        print*, 'Pb. avec opacite TOPS, log do16/dc12= ',
     +           log_do16sdc12
      end if

      if( log_do16sdc12 .gt. val_o16(n_o16) ) then
        print*, 'Pb. avec opacite TOPS, log do16/dc12= ',
     +           log_do16sdc12
      end if

c Interpolation par rapport a log10 ro

      do j_y= 1, n_y
         do j_o16= 1, n_o16
            do jt= 1, nt, 1

      do jr= 1, nro
         fro(jr)=vlk(j_y,j_o16,jt,jr)
      end do

c        10        20        30        40        50        60        70

      call int_ext_parab(fro,vlro,nro,lro,kro(j_y,j_o16,jt),
     +                   dkro(j_y,j_o16,jt))
 
            end do
          end do
      end do

c      write(6,*) 'OK 2'

c Interpolation par rapport a log10 t6


      do j_y= 1, n_y
         do j_o16= 1, n_o16

         do jt= 1, nt
         ft(jt)=kro(j_y,j_o16,jt)
         end do

        call int_ext_parab(ft,vlt6,nt,lt6,krot(j_y,j_o16),
     +                     dkrot(j_y,j_o16))

        end do
      end do

c      print*, 'OK 3'

c Construction du tableau 'dkrob' (valeurs de la derivee /ro)
c C'est la aussi une interpolation par rapport a t6


      do j_y= 1, n_y
         do j_o16= 1, n_o16

         do jt= 1, nt
         ft(jt)=dkro(j_y,j_o16,jt)
         end do

         call int_ext_parab(ft,vlt6,nt,lt6,dkrob(j_y,j_o16),bidon)

         end do
      end do

c      print*, 'OK 4'

c Interpolation de log10 kappa par rapport a log10( DO16 /D12 )

      do j_y= 1, n_y

         do j_o16= 1, n_o16
            ft(j_o16)=krot(j_y,j_o16)
         end do

         call int_ext_parab(ft,val_o16,n_o16,log_do16sdc12,kroto16(j_y),
     +                      bidon)

      end do

c Interpolation de dlog kap /dlog ro par rapport a log( DO16 /DC12 )

      do j_y= 1, n_y

         do j_o16= 1, n_o16
            ft(j_o16)=dkrob(j_y,j_o16)
         end do

      call int_ext_parab(ft,val_o16,n_o16,log_do16sdc12,dkrobo16(j_y),
     +                   bidon)

      end do

c Interpolation de dlog kap / dlog t6 par rapport a log( DO16 /DC12 )

      do j_y= 1, n_y

         do j_o16= 1, n_o16
            ft(j_o16)=dkrot(j_y,j_o16)
         end do

      call int_ext_parab(ft,val_o16,n_o16,log_do16sdc12,dkroto16(j_y),
     +                   bidon)

      end do

c Interpolation de log kappa par rapport a Y

      call int_ext_parab(kroto16,val_y,n_y,y,lkappa,bidon)

c Interpolation de dlog kap / dlog t6 par rapport a Y

      call int_ext_parab(dkroto16,val_y,n_y,y,dkappadt,bidon)

c Interpolation de dlog kap /dlog ro par rapport a Y

      call int_ext_parab(dkrobo16,val_y,n_y,y,dkappadro,bidon)

c Diverses conversions :

      kap=10.**lkappa

      dkapt=dkappadt*kap/t

      dkapro=dkappadro*kap/ro

      return

      end

c-----------------------------------------------------------------------

      subroutine lect_tops_4(file_opa, n_y, n_o16, nt, nr, z0,
     +                       val_y, val_o16, vlt6, vlr, vlk)

c Routine de lecture du fichier binaire d'opacite 'file_opa'
c pour les coeurs de Cepheides

c Entree : file_opa : nom du fichier binaire de data

c Sorties : n_y     : nombre de valeurs de Y
c           n_o16   : nombre de valeurs de log10(DO16/DC12)
c           nt      : nombre de valeurs de la temperature
c           nr      : nombre de valeus de ro
c           z0      : metallicite initiale du modele
c           val_y   : les valeurs de Y
c           val_o16 : les valeurs de log10(DO16/C12)
c           vlt6    : valeurs de log10(t6=t/10.**6)
c           vlr     : valeurs de log10(ro) ATTENTION ici ce n'est pas R !!!
c           vlk     : valeurs de log10(kappa)           

c Daniel Cordier, Mai 1998

      implicit none

      integer iunit_lect
      integer pnf, pnt, pnr, pn_y, pn_o16
      integer nf, nt, nr, n_y, n_o16
      integer jf, jt, jr, j_y, j_o16

      character*50 file_opa

      parameter( pnf=12, pnt=26, pnr=30, pn_y=4, pn_o16=3 )

      real*8 z0, val_y(pn_y), val_o16(pn_o16), vlt6(pnt), vlr(pnr), 
     +       k(pnf,pnt,pnr), kap(pn_y,pn_o16,pnt,pnr),
     +       vlk(pn_y,pn_o16,pnt,pnr)

c Ouverture du fichier d'opacite

      iunit_lect=71
      open(iunit_lect,form='unformatted',file=file_opa)

c Lecture des tables

      read(iunit_lect) nf, n_y, n_o16, nt, nr
      read(iunit_lect) z0
      read(iunit_lect) (val_y(j_y), j_y= 1, n_y)
      read(iunit_lect) (val_o16(j_o16), j_o16= 1, n_o16)

      if( nf .gt. pnf) then
        print*, 'Pb. avec ''lect_tops'' : nf .gt. pnf !'
        stop
      end if
      if( n_y .gt. pn_y) then
        print*, 'Pb. avec ''lect_tops'' : n_y .gt. pn_y !'
        stop
      end if
      if( n_o16 .gt. pn_o16 ) then
        print*, 'Pb. avec ''lect_tops'' : n_o16 .gt. pn_o16 !'
        stop
      end if
      if( nf .ne. n_y*n_o16 ) then
        print*, 'Pb. avec ''lect_tops'' : n_y*n_o16 .ne. nf !'
        stop
      end if

      write(6,*) 'Tables d''opacite << TOPS >> pour coeur de Cepheide'
      write(6,*) 'Version 4.0'
      write(6,*) 'X=0., Z0= ', z0, ' , log10(DNe22/DC12)= -1.4'
      write(6,*) n_y, ' valeurs de Y=X_He4 :'
      write(6,*) (val_y(j_y), j_y= 1, n_y)
      write(6,*) n_o16, ' valeurs de log10(DO16/DC12) :'
      write(6,*) (val_o16(j_o16), j_o16= 1, n_o16)
      write(6,*) nt, ' valeurs de t6'
      write(6,*) nr, ' valeurs de ro'

      read(iunit_lect) (vlt6(jt), jt= 1, nt)
      read(iunit_lect) (vlr(jr), jr= 1, nr)

      read(iunit_lect) (((k(jf,jt,jr),jr=1,nr),jt=1,nt),
     +                    jf=1,nf)

      close(iunit_lect)

c Rearangement dans un tableau ou la dimension "numero du fichier TOPS"
c est change en deux dimensions : j_y "valeur de Y" et j_o16
c "valeur de log D016 /DC12"

      do j_y= 1, n_y
         do j_o16= 1, n_o16
            do jt= 1, nt
               do jr= 1, nr
      kap(j_y,j_o16,jt,jr)=k((j_y-1)*n_o16+j_o16,jt,jr)
               end do
            end do
         end do
      end do

c Conversion kappa ---> log10 kappa

      do j_y= 1, n_y
         do j_o16= 1, n_o16
            do jt= 1, nt
               do jr= 1, nr
      vlk(j_y,j_o16,jt,jr)=log10(kap(j_y,j_o16,jt,jr))
               end do
            end do
          end do
      end do

      end


c------------------------------------------------------------------------
