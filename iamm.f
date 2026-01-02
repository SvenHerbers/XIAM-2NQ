C----------------------------------------------------------------------
      subroutine calcm(sigma,h,evalv,ovv,mvec
     $     ,am,qmv,ifit,k,maxm,minv,sizev,ai)
C     calculation of the eigenvalues of one matrix with specified sigma
C     the evalues are put in the field of evalv(1..sizev)
C     the matrix_elements are in ovv(1..sizev,1..sizev,x)
      implicit none
      include 'iam.fi'
      integer sigma,k,maxm,minv,sizev
      real*8  h(DIMTOT,DIMTOT),evalv(DIMV),ovv(DIMV,DIMV,DIMOVV)
      real*8  am(DIMPM), ai(DIMPIR)!Herbers2026
      real*8  mvec(DIMM,DIMV)
      integer qmv(DIMV),ifit(DIMOVV)

C     work
      real*8  e(DIMTOT),d(DIMTOT)
      real*8  ao(DIMOVV)
      integer qm(DIMTOT)
      integer sizem,ir,ic,i,ierr

      sizem=2*maxm+1
      if (sizem.lt.sizev+minv)  stop 'ERROR: size m < size v'
      if (sizev.gt.DIMV)   stop 'ERROR: sizev > DIMV in calcm'
      if (sizem.gt.DIMM)   stop 'ERROR: sizem > DIMM in calcm'
      if (sizem.gt.DIMTOT) stop 'ERROR: sizem > DIMTOT in calcm'
      do i=1, DIMPM
        ao(i)=am(i)
      end do

c      write(*,*) (am(ic),ic=1,DIMPM)
c      write(*,*) (qmv(ic),ic=1,DIMV)
      call buildm(sigma,h,am,ao,qm,sizem,k,ai)
c      write(*,'(2I3,40F12.4)') k,sigma,(h(ic,ic),ic=1,sizem)
      call hdiag(DIMTOT,sizem,h,d,e,ierr)
      if (ierr.ne.0) then
        write(0,'(A,3I3)') 'Error: HDIAG in CALCM: IERR =',ierr,k,sigma
      end if
      call eigsrt(d,h,sizem,DIMTOT)
      call phasem(sigma,h,d,qm,sizem,sizev,minv,qmv,k)
      call assgnm(sigma,h,d,qm,sizem,sizev,minv,qmv,k)
      call calovv(ifit,am,qm,sizem,sizev,minv,h,ovv,k,ai)

      do ir=1, sizev
        evalv(ir)=d(ir+minv-1)
C       ovv(ir,ir,PM_E)=d(ir+minv-1)
      end do
      do ir=1, sizem
        do ic=1, sizev
          mvec(ir,ic)=h(ir,ic+minv-1)
        end do
      end do
      return
      end      

C----------------------------------------------------------------------
      subroutine buildm(sigma,h,am,ao,qm,sizem,k,ai)
      implicit none
      include 'iam.fi'
      integer sigma, qm(DIMTOT), sizem, k
      real*8  h(DIMTOT,DIMTOT)
      real*8  am(DIMPM),ao(DIMOVV), ai(DIMPIR)!Herbers2026
 
C     work
      integer im,iv,mm
      real*8  v(DIMTOT), vo(DIMTOT)

      if (sizem.gt.DIMTOT) stop 'Dimension Error in BUILDM'

C     initialize the quantum no  qm
      mm=sizem/2
      do im=1, sizem
        qm(im)=ctlint(C_NFOLD)*(im-mm-1)+sigma
      end do

      do im=DIMPM+1, DIMOVV
        ao(im)=0.0
      end do
      ao(PM_RHO)=0.0
      do im=1, sizem
        do iv=1, sizem
          v(iv)=0.0
          vo(iv)=0.0
        end do
        v(im)=1.0
        call multm(v,sizem,am,ao,qm,vo,k,0,ai)
        do iv=1, sizem
          h(im,iv)=vo(iv)
        end do
      end do
      return
      end

C----------------------------------------------------------------------
      subroutine multm(v,sizem,am,ao,qm,vo,k,ip,ai)
      implicit none
      include 'iam.fi'
      integer sizem,k,ip
      integer qm(DIMTOT)
      real*8  v(DIMTOT), vo(DIMTOT)
      real*8  am(DIMPM), ao(DIMOVV), ai(DIMPIR)!Herbers2026
C     work
      real*8  dm,t,dk,dm1
      integer im,off

      if (sizem.gt.DIMTOT) stop 'Dimension Error in MULTM'
      
C     diagonal m/m 
      dk=dble(k)
      do im=1, sizem
        dm=dble(qm(im))
        t =
     $       + (ao(PM_F  ) 
     $       + ai(PI_FMK)*am(PM_RHO)*dk*dm
     $       + ai(PI_DFM2)*dm**2) ! ai() Parameters are discarded in calcovv
C     $       + ai(PI_FMK)*am(PM_RHO)**2*dk**2)
     $       * (dm - am(PM_RHO)*dk)**2  
     $       + ao(PM_VN1) * 0.5
     $       + ao(PM_VN2) * 0.5
     $       + ao(PM_RHO) * 2.0*am(PM_F)*dk*(am(PM_RHO)*dk-dm)
c     $                      +4.0*am(PM_DPI4)*dk*(am(PM_RHO)*dk-dm)**3)
     $       + ao(PM_PI ) * (dm - am(PM_RHO)*dk)
     $       + ao(PM_DPI4)* (dm - am(PM_RHO)*dk)**4
     $    + ao(PM_MK3)* am(PM_RHO)**3*dk**3*dm
     $    + ao(PM_M3K)* am(PM_RHO)*dk*dm**3
C     $       + ao(PM_M3K)* (dm - am(PM_RHO)*dk)**2*dm**2
C     $    + ao(PM_MK3)* (dm - am(PM_RHO)*dk)**2*am(PM_RHO)**2*dk**2
c     $       + ao(PM_M)   *  dm
c     $       - ao(PM_RK)  *  am(PM_RHO)*dk
        vo(im)=vo(im)+v(im)*t
      end do
      
C     off diagonal m/m+1 
      if ((ao(PM_VN1).ne.0.0)) then
        off=1
        do im=1, sizem-off
          t=   - ao(PM_VN1)*0.25
          vo(im)    =vo(im    )+v(im+off)*t
          vo(im+off)=vo(im+off)+v(im    )*t
        end do
      end if
C     off diagonal m/m+2 
      if ((ao(PM_VN2).ne.0.0)) then
        off=2
        do im=1, sizem-off
          t=   - ao(PM_VN2)*0.25
          vo(im)    =vo(im    )+v(im+off)*t
          vo(im+off)=vo(im+off)+v(im    )*t
        end do
      end if
C     off diagonal m/m+1   i * (sin n alpha)
      if ((ao(PM_SIN).ne.0.0)) then
        off=1
        do im=1, sizem-off
          t=   ao(PM_SIN)*0.5
          vo(im)    =vo(im    )-v(im+off)*t
          vo(im+off)=vo(im+off)+v(im    )*t
        end do
      end if
C     off diagonal m/m+1   (cos n alpha)
      if ((ao(PM_COS).ne.0.0)) then
        off=1
        do im=1, sizem-off
          t=   ao(PM_COS)*0.5
          vo(im)    =vo(im    )+v(im+off)*t
          vo(im+off)=vo(im+off)+v(im    )*t
        end do
      end if
c      if ((ao(PM_DPIC).ne.0.0).or.(ao(PM_RHO).ne.0.0)) then
c        off=1
c        dm=dble(qm(im))
c        dm1=dble(qm(im+off))
c        do im=1, sizem-off
c          t= ao(PM_DPIC)*((dm-am(PM_RHO)*dk)**2+(dm1-am(PM_RHO)*dk)**2)
c     $        -ao(PM_RHO)*2.0*am(PM_DPIC)*dk*
c     $         ((dm-am(PM_RHO)*dk)+(dm1-am(PM_RHO)*dk))
c          vo(im)    =vo(im    )+v(im+off)*t
c          vo(im+off)=vo(im+off)+v(im    )*t
c        end do
c      end if

      return        
      end

C---------------------------------------------------------------------- 
      subroutine assgnm(sigma,h,eval,qm,sizem,sizev,minv,qmv,k)
      implicit none
      include 'iam.fi'
      real*8 h(DIMTOT,DIMTOT),eval(DIMTOT)
      integer sigma,qm(DIMTOT),qmv(DIMV)
      integer sizem, sizev,k,minv
C     work
      integer icc
      integer bestm(DIMV), scndm(DIMV)
      real*8  besth(DIMV), scndh(DIMV)
      real*8  difh,degn,sym

      if (sizem.gt.DIMTOT) stop ' ERROR: DIM 2M1 exceeded'
      degn=2.0d-4

      call maxof(h,DIMTOT,DIMTOT,sizem,1,minv+sizev-1,bestm,scndm,
     $       besth,scndh)
      do icc=1, sizev
        difh=1.0
        if (sizem.gt.1) then
          difh=(dabs(besth(icc))-dabs(scndh(icc)))
        end if
        if (difh.lt.degn) then
C         if ((besth(icc).lt.0.0d0).and.(scndh(icc).lt.0.0d0)) then
C           do im=1, sizem
C             h(im,icc)=-h(im,icc)
C           end do
C         end if
          sym=dsign(10.0d0,besth(icc)*scndh(icc))
          qmv(icc)=isign(qm(bestm(icc)),int(sym))
        else
C         if (besth(icc).lt.0.0d0) then
C           do im=1, sizem
C             h(im,icc)=-h(im,icc)
C           end do
C         end if
          qmv(icc)=qm(bestm(icc))
        end if
  
      end do  
      return
      end

C---------------------------------------------------------------------- 
      subroutine phasem(sigma,h,eval,qm,sizem,sizev,minv,qmv,k)
      implicit none
      include 'iam.fi'
      real*8 h(DIMTOT,DIMTOT),eval(DIMTOT)
      integer sigma,qm(DIMTOT),qmv(DIMV)
      integer sizem, sizev,k, minv
C     work
      integer iv,im,ilim
      real*8  sum,maxh,sigsn

      if (sizem.gt.DIMTOT) stop ' ERROR: DIM 2M1 exceeded'
      if ((k.ne.0).or.(sigma.ne.0)) then
        ilim=sizem
      else
        ilim=sizem/2+1
      end if
      sigsn=1.0
      if (sigma.lt.0) sigsn=-1.0
C      do iv=minv, minv+sizev-1
      iv=1
      maxh=0.0
      do im=1,ilim
        if (abs(maxh).lt.abs(h(im,iv))) then
          maxh=h(im,iv)
        end if
      end do
      if (maxh.lt.0.0) then
        do im=1,sizem
          h(im,iv)=-h(im,iv)
        end do
      end if
C      end do
C     calculate < v | p | (v+1) >
C     sigma=0,1 <p> .gt. 0
C     sigma=-1  <p> .lt. 0
C      do iv=1,sizev-1
C     do iv=minv,minv+sizev-2      
      do iv=1,minv+sizev-2
        sum=0.0
        do im=1,sizem
          sum=sum+h(im,iv)*h(im,iv+1)*qm(im)
        end do
        sum=sum*sigsn
        if (sum.lt.0.0) then
          do im=1,sizem
            h(im,iv+1)=-h(im,iv+1)
          end do
        end if
      end do
  
      return
      end

C----------------------------------------------------------------------
      subroutine calovv(ifit,am,qm,sizem,sizev,minv,h,ovv,k,ai)
C     This Subroutine calculates the operator representations in the RAS.
C     Ai should not be used here.
      implicit none
      include 'iam.fi'
      integer sizem,sizev,k,minv
      integer qm(DIMTOT),ifit(DIMOVV)
      real*8  am(DIMPM),h(DIMTOT,DIMTOT),ovv(DIMV,DIMV,DIMOVV)
      real*8  ai(DIMPIR)!Herbers2026
C     work
      real*8  vo(DIMTOT),da(DIMOVV),t(DIMTOT,DIMV)
      real*8  noai(DIMPIR)!Herbers2026
      integer i,it,im,ir,ic
      do i=1, DIMPIR !Herbers2026
        noai(i) = 0.0 !Herbers2026
      end do          !Herbers2026
      if (sizem.gt.DIMTOT) stop 'Dimension Error in CALOVV'
C     clear ovv
      do ir=1, sizev
        do ic=1, sizev
          do i=1, DIMOVV
            ovv(ir,ic,i)=0.0
          end do
        end do
      end do
      do it=1, DIMOVV
        if (it.eq.PM_PI2) goto 99
        do i=1,DIMOVV
          da(i)=0.0
        end do
        da(it)=1.0
        do i=1, sizev
          do im=1, sizem
            vo(im)=0.0
          end do
          call multm(h(1,i+minv-1),sizem,am,da,qm,vo,k,1,noai)
          do im=1, sizem
            t(im,i)=vo(im)
          end do
        end do
        do ir=1, sizev
          do ic=1, sizev
            do im=1, sizem
              ovv(ir,ic,it)=ovv(ir,ic,it)+h(im,ir+minv-1)*t(im,ic)
            end do
          end do
        end do
      end do
      return
c --------
 99   continue
      do ir=1, sizev
        do ic=1, sizev
          ovv(ir,ic,it)=ovv(ir,ic,PM_PI)**2
        end do
      end do
      return 
      end

