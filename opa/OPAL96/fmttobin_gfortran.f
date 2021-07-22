      programme fmttobin
  
c     traduit une table d'opacite du binaire en ASCII
c     ou l'inverse (Yveline- octobre 1996)
c---------------------------------------------------------
  
      implicit none
c 
      integer pnz,pnx,pnt,pnr
      parameter(pnz=16,pnx=8,pnt=85,pnr=23)
  
      integer nt,nr,nx,nz,ix,iz,it,ir
  
      real*8 T6(pnt),logR(pnr),X(pnx),Z(pnz),kap(pnz,pnx,pnt,pnr)
  
      character answer*1, nomasc*40, nombin*40
  
c-------------------------------------------------------------------
      write(6,*) 'traduction ASCII ---> binaire (y/n?)     '
      read(5,'(a)') answer
  
      if(answer.eq.'y') then
        write(6,*) 'ASCII ---> binaire'
        write(6,*)'name of the existing ASCII opacity table : '
        read(5,'(a)') nomasc
        write(6,*) 'name of the new binary table: '
        read(5,'(a)') nombin

c Instruction d'origine d'Yveline :
c        open (unit=8,form='formatted',status='old',readonly,
c     &        name=nomasc)
c        open (unit=9,form='unformatted',status='new',
c     &        name=nombin)



        open (unit=8,form='formatted',status='old',
     &        file=nomasc)
        open (unit=9,form='unformatted',status='new',
     &        file=nombin)


        read(8,*) nz,nx,nt,nr
        read(8,*) (T6(it),it=1,nt)
        read(8,*) (logR(ir),ir=1,nr)
        write(9)  nz,nx,nt,nr
        write(9)  (T6(it),it=1,nt)
        write(9)  (logR(ir),ir=1,nr)
        do 1 iz=1,nz
          read(8,*) Z(iz)
          write(9)  Z(iz)
          do 2 ix=1,nx
            read(8,*) X(ix)
            read(8,*) ((kap(iz,ix,it,ir),ir=1,nr),it=1,nt)
            write(9)  X(ix)
            write(9)  ((kap(iz,ix,it,ir),ir=1,nr),it=1,nt)
2         continue
1       continue
      else
        write(6,*) 'binaire ---> ASCII'      
        write(6,*) 'name of the existing binary table: '
        read(5,'(a)') nombin
        write(6,*)'name of the new ASCII opacity table : '
        read(5,'(a)') nomasc

c Instructions d'origine d'Yveline :
c        open (unit=8,form='formatted',status='new',
c     &        name=nomasc)
c        open (unit=9,form='unformatted',status='old',readonly,
c     &        name=nombin)

        open (unit=8,form='formatted',status='new',
     &        file=nomasc)
        open (unit=9,form='unformatted',status='old',
     &        file=nombin)


        read(9)  nz,nx,nt,nr
        read(9)  (T6(it),it=1,nt)
        read(9)  (logR(ir),ir=1,nr)
        write(8,*) nz,nx,nt,nr
        write(8,*) (T6(it),it=1,nt)
        write(8,*) (logR(ir),ir=1,nr)
        do 3 iz=1,nz
          read(9)  Z(iz)
          write(8,*) Z(iz)
          do 4 ix=1,nx
            read(9) X(ix)
            read(9) ((kap(iz,ix,it,ir),ir=1,nr),it=1,nt)
            write(8,*) X(ix)
            write(8,*) ((kap(iz,ix,it,ir),ir=1,nr),it=1,nt)
4         continue
3       continue
      endif
  
      close(8)
      close(9)
  
      stop
      end      
  
