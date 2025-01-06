C----------------------------------------------------------------------
      subroutine adjusta(a,npar,adj)
C     adjust the elements in a for consistence using adj
C     a is a simple vector
      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR)
      integer npar,adj
C     local..
      real*8  f(DIMTOP,DIMTOP),pi
      real*8  tmprho,tmpgam,tmpbet
      real*8  f0(DIMTOP),lx(DIMTOP),ly(DIMTOP),lz(DIMTOP)
      real*8  rho(DIMTOP),beta(DIMTOP),agam(DIMTOP)
      integer itop,ift
      integer myand
      external myand
      pi=dacos(-1.0d0)

c      write(*,*)'adjusta start',adj,a(P1_F)
      do itop=1, ctlint(C_NTOP)
        ift=(itop-1)*DIMPIR

C     adjust 16: obtain beta and gamma from angx and angz 
        if (myand(adj,16).gt.0) then
          call reclg(a(P_BJ),a(P_BK),a(P_BD)
     $         ,tmpbet,tmpgam
     $         ,a(P1_ANGX+ift),a(P1_ANGZ+ift))
          if ((myand(ctlint(C_PRI),AP_PC).ne.0).and.(xde.ge.1))
     $         write(*,'(2(A,F12.6))')
     $         ' adjusting beta  from',a(P1_BETA+ift)
     $         ,' to',tmpbet ,
     $         ' adjusting gamma from',a(P1_GAMA+ift)
     $         ,' to',tmpgam 
          a(P1_BETA+ift)=tmpbet
          a(P1_GAMA+ift)=tmpgam
        end if
        a(P1_BETA+ift)=dmod(a(P1_BETA+ift),2.0d0*PI)        
        a(P1_GAMA+ift)=dmod(a(P1_GAMA+ift),2.0d0*PI) 
c        write(*,*) ift,a(P1_ANGX+ift),a(P1_ANGZ+ift)
c     $       ,a(P1_BETA+ift),a(P1_GAMA+ift)

C     adjust 8 : obtain rho from F0 (= 1 / I_alpha)
        if (myand(adj,8).gt.0) then
          call recrho(a(P_BJ),a(P_BK),a(P_BD)
     $         ,a(P1_BETA+ift),a(P1_GAMA+ift)
     $         ,a(P1_F0+ift),tmprho)
          if ((myand(ctlint(C_PRI),AP_PC).ne.0).and.(xde.ge.1))
     $         write(*,'(2(A,F12.6))')
     $         ' adjusting rho from',a(P1_RHO+ift)
     $         ,' to',tmprho 
          a(P1_RHO+ift)=tmprho
        end if
      end do

C     adjust 1,2, and 4  
      do itop=1,ctlint(C_NTOP)
        ift=(itop-1)*DIMPIR
        rho(itop)  =a(P1_RHO+ift)
        beta(itop) =a(P1_BETA+ift)
        agam(itop) =a(P1_GAMA+ift)
      end do
      call recaff(a(P_BJ),a(P_BK),a(P_BD)
     $       ,rho,beta,agam
     $       ,f,f0,lx,ly,lz,ctlint(C_NTOP),DIMTOP
     $       ,myand(adj,4))
c      write(*,*)'adjusta     a',adj,a(P1_F)

      do itop=1,ctlint(C_NTOP)
        ift=(itop-1)*DIMPIR
        if (myand(adj,1).gt.0) a(P1_F+ift)=f(itop,itop)
        if (myand(adj,8).eq.0) a(P1_F0+ift)=f0(itop)
        if (myand(adj,16).eq.0) then
          a(P1_ANGZ+ift)=acos(lz(itop))
          if ((lx(itop)**2+ly(itop)**2).gt.(0.0d0)) then 
            a(P1_ANGX+ift)=sign(1.0d0,ly(itop))*acos(lx(itop)
     $           /sqrt(lx(itop)**2+ly(itop)**2))
          else
            a(P1_ANGX+ift)=0.0d0
          end if
        end if
        a(P1_ANGX+ift)=dmod(a(P1_ANGX+ift),2.0d0*PI)        
        a(P1_ANGZ+ift)=dmod(a(P1_ANGZ+ift),2.0d0*PI) 
      end do
      if (myand(adj,2).gt.0) a(P_FF)=f(1,2)
      return
      end

C----------------------------------------------------------------------
      subroutine recalf(bj,bk,bd,rho,beta,gamma,f0,lx,ly,lz)
C     calc. f0 and angles (lx ly lz) from IAM-Parameters rho beta gamma
      implicit none
      real*8 bj,bk,bd,f0,rho,beta,gamma,lx,ly,lz
      real*8 bx,by,bz,rhox,rhoy,rhoz
      bx=bj+bd
      by=bj-bd
      bz=bj+bk
      if (rho.eq.0.0) stop 'ERROR: rho (or F0) can not be zero'
      rhoz=cos(beta)*rho
      rhox=sin(beta)*cos(gamma)*rho
      rhoy=sin(beta)*sin(gamma)*rho
      f0=1.0/sqrt((rhox/bx)**2+(rhoy/by)**2+(rhoz/bz)**2)
      lx=rhox*f0/bx
      ly=rhoy*f0/by
      lz=rhoz*f0/bz
C      r=1-bz/f0*lz**2-bx/f0*lx**2-by/f0*ly**2
C      f=f0/r
      return
      end

C----------------------------------------------------------------------
      subroutine recrho(bj,bk,bd,beta,gamma,f0,rho)
C     calc rho from F0, beta and gamma
      implicit none
      real*8 bj,bk,bd,f0,rho,beta,gamma
      real*8 bx,by,bz,rz,ry,rx
      if (f0.eq.0.0) stop 'ERROR in recrho: F0 = 0' 
      bx=bj+bd
      by=bj-bd
      bz=bj+bk
      rz=cos(beta)*f0/bz
      rx=sin(beta)*cos(gamma)*f0/bx
      ry=sin(beta)*sin(gamma)*f0/by
      rho=1.0/sqrt(rx**2+ry**2+rz**2)
      return
      end
      
C----------------------------------------------------------------------
      subroutine reclg(bj,bk,bd,beta,gamma,ax,az)
C     calc beta, and gamma from ax,az
      implicit none
      real*8 bj,bk,bd,beta,gamma,lx,ly,lz,ax,az
      real*8 bx,by,bz,rz,ry,rx,ac
      lz=cos(az)
      lx=sin(az)*cos(ax)
      ly=sin(az)*sin(ax)
      bx=bj+bd
      by=bj-bd
      bz=bj+bk
      rx=lx*bx
      ry=ly*by
      rz=lz*bz
      beta=acos(rz/sqrt(rx**2+ry**2+rz**2))
      ac=rx/sqrt(rx**2+ry**2)
      ac=min(ac,1.0d0)
      ac=max(ac,-1.0d0)
      gamma=sign(1.0d0,ry)*acos(ac)
      return
      end
      
C----------------------------------------------------------------------
      subroutine recaff(bj,bk,bd,rho,beta,gamma,
     $     f,f0,lx,ly,lz,ntop,DIMTOP,singl)
      implicit none
      integer ntop,DIMTOP
      real*8 bj,bk,bd
      real*8 rho(DIMTOP),beta(DIMTOP),gamma(DIMTOP)
      real*8 f(DIMTOP,DIMTOP)
      real*8 f0(DIMTOP)
      real*8 lx(DIMTOP),ly(DIMTOP),lz(DIMTOP)
C     local ..
      real*8 bx,by,bz,fmat(4,4),fimat(4,4)!Herbers2024
      real*8 d(4),e(4)
      integer itp1,itp2,singl

      if (DIMTOP.gt.4) stop ' Dimension Error in recaff !'
      bx=bj+bd!Herbers2024 .gt.4
      by=bj-bd
      bz=bj+bk

      do itp1=1,ntop
        call recalf(bj,bk,bd,rho(itp1),beta(itp1),gamma(itp1)
     $       ,f0(itp1),lx(itp1),ly(itp1),lz(itp1))
      end do

      do itp1=1,ntop      
        do itp2=1,ntop
          fmat(itp1,itp2)=0.0
        end do
      end do
      do itp1=1,ntop
        fmat(itp1,itp1)=(f0(itp1)
     $       -bx*lx(itp1)**2
     $       -by*ly(itp1)**2
     $       -bz*lz(itp1)**2)/(f0(itp1)**2)
      end do
      if (singl.eq.0) then
        do itp1=1,ntop
          do itp2=itp1+1,ntop
            fmat(itp1,itp2)=
     $           -(lx(itp1)*lx(itp2)*bx
     $           +ly(itp1)*ly(itp2)*by
     $           +lz(itp1)*lz(itp2)*bz)/(f0(itp1)*f0(itp2))
            fmat(itp2,itp1)=fmat(itp1,itp2)
          end do
        end do
      end if
C      write(*,*)'pre-Inverse'
C      do itp1=1, ntop
C        write(*,*) (fmat(itp1,itp2),itp2=1,ntop)
C      end do
      call svdsydc(fmat,d,e,ntop,4)!Herbers2024
      call svdsyinv(fmat,d,ntop,4,fimat)!Herbers2024
      do itp1=1,ntop
        do itp2=1,ntop
          f(itp1,itp2)=fimat(itp1,itp2)
        end do
      end do
C      write(*,*)'Inverse'
C      do itp1=1, ntop
C        write(*,*) (f(itp1,itp2),itp2=1,ntop)
C      end do
      
      return
      end
      
