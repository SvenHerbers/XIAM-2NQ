C----------------------------------------------------------------------
      subroutine bld2vjk(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,complex)
C     this subroutine builds the rotation + internal rotation matrix 
C     h(DIMTOT,DIMTOT). 
C      The energie E_{k,v,sigma} (given in the array evalv(v,sigma,k,top))
C      for each top and the operators, which will be needed for the top-top
C      coupling terms, 
C      (given in the array ovv(v,v',op_no,sigma,k,top)) are transformed into
C      the principal axes system with the rotation matrix rotm(k,k',1,top)
C      and written in the array rott(k,k',v,v',top). 
C      this is provided by the subroutines are rotovv (for the operators) 
C      and roteval (for the eigenenergies).
C      the second step is to write the arrays rott (where each top in independent) 
C      into the matrix h using the subroutines addo1 and addovv.

      implicit none
      include 'iam.fi'
      integer j,gam,f
      logical complex
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C      real*8            erk(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      integer itop,ift
      logical rrir1,rrir2,rrir3,rrir4 !Herbers2024
      external myand
      integer myand

      if (ctlint(C_NTOP).ge.2) then
        if (a(P_FF).ne.0.0) then
          call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,PM_PI,PM_PI
     $         ,.false.)
          if (complex) then
            call haddovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,a(P_FF))
          else
            call addovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,a(P_FF))
          end if
        end if
        if (a(P_VSS).ne.0.0) then 
          if (myand(ctlint(C_WOODS),1024).eq.0) then
            call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,PM_SIN,PM_SIN
     $           ,.false.)
            if (complex) then
              call haddovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,-5.d-1*a(P_VSS))
            else
              call addovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,-5.d-1*a(P_VSS))
            end if
          else
            call addiop(j,gam,f,qvk,ruse,h,a,ovv,rotm,rott,tori
     $           ,-1.0*a(P_VSS),PM_SIN,.true.)
          end if
        end if
        if (a(P_VCC).ne.0.0) then 
          if (myand(ctlint(C_WOODS),1024).eq.0) then
            call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,PM_COS,PM_COS
     $           ,.false.)
            if (complex) then
              call haddovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,+5.d-1*a(P_VCC))
            else
              call addovv(j,gam,f
     $             ,qvk,ruse,h,a,ovv,rotm,rott,tori,+5.d-1*a(P_VCC))
            end if
          else
            call addiop(j,gam,f,qvk,ruse,h,a,ovv,rotm,rott,tori
     $           ,a(P_VCC),PM_COS,.true.)
          end if
        end if
      end if
      rrir1=.false.
      rrir2=.false.
      rrir3=.false.
	  rrir4=.false.!Herbers2024
      do itop=1, ctlint(C_NTOP)
        ift=(itop-1)*DIMPIR
        if (a(P1_DPIJ+ift).ne.0.0) rrir1=.true.
        if (a(P1_DPIK+ift).ne.0.0) rrir2=.true.
        if (a(P1_DPID+ift).ne.0.0) rrir3=.true.
		if (a(P1_DPK2+ift).ne.0.0) rrir4=.true.!Herbers2024
c        write(*,'(50F10.4)')
c     $       (ovv(1,1,PM_PI2,gamma(gam,itop),ift,itop),ift=-j,j)
      end do
      if (rrir1.or.rrir2.or.rrir3.or.rrir4)
     $     call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,PM_F,PM_F
     $     ,.true.)
      if (rrir1) then
        call haddjp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori
     $       ,PI_DPIJ)!  
      end if
      if (rrir2) then
        call haddkp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori
     $       ,PI_DPIK)
      end if
      if (rrir3) then
        call hadddp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori
     $       ,PI_DPID)
      end if
      if (rrir4) then!Hebers2024
        call haddkkp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori
     $       ,PI_DPK2)!Herbers2024
      end if!Herbers2024
      rrir1=.false.
      rrir2=.false.
      rrir3=.false.
	  rrir4=.false.!Herbers2024
      do itop=1, ctlint(C_NTOP)
        ift=(itop-1)*DIMPIR
        if (a(P1_DC3J+ift).ne.0.0) rrir1=.true.
        if (a(P1_DC3K+ift).ne.0.0) rrir2=.true.   !Herbers2018
        if (a(P1_DC3D+ift).ne.0.0) rrir3=.true.   !Herbers2018
		if (a(P1_D3K2+ift).ne.0.0) rrir4=.true.   !Herbers2024
C       ...
C       ...
      end do
      if (rrir1.or.rrir2.or.rrir3.or.rrir4)                                              !Herbers
     $     call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,PM_COS,PM_COS
     $     ,.true.)
      if (rrir1) then
        call haddjp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori
     $       ,PI_DC3J)! 
      end if
      if (rrir2) then                                                       !Herbers
        call haddkp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori           !Herbers
     $       ,PI_DC3K)                                                      !Herbers
      end if                                                                !Herbers
      if (rrir3) then                                                       !Herbers
        call hadddp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori           !Herbers
     $       ,PI_DC3D)                                                      !Herbers
      end if                                                                !Herbers
      if (rrir4) then                                                       !Herbers
        call haddkkp(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori   !Herbers2024
     $       ,PI_D3K2)                                               !Herbers2024
      end if 
      call rotevl(j,gam,f,qvk,a,evalv,rotm,rott,tori)
      if (complex) then
        call haddo1(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,1.0d0)
      else
        call addo1(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,1.0d0)
      end if
      return
      end 

C ---------------------------------------------------------------------
      subroutine rotevl(j,gam,f,qvk,a,evalv,rotm,rott,tori)
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      real*8  a(DIMPAR)
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      integer itop,ikr,ikc,ivr,ivc,ik,iv,itest,off
      external myand
      integer myand
      real*8  tt

      do itop=1, ctlint(C_NTOP)
        do ikr=1, size(S_K)
          do ikc=1, size(S_K)
            do ivr=1, size(S_V+itop)
              do ivc=1, size(S_V+itop)
                rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)=0.0
              end do
            end do
          end do
        end do
      end do
      itest=0

C     rotate evalv into rho system for each top with torsional integrals
C      if ((myand(ctlint(C_WOODS),1).ne.0).and.      ! torsional integrals
C     $     (myand(ctlint(C_WOODS),4).ne.0).and.     ! use tors int. in rot. matrix
C     $     (myand(ctlint(C_WOODS),64).eq.0)) then   ! don't use Demaisons method
      if (myand(ctlint(C_WOODS),4).ne.0) then ! torsional integrals
        itest=itest+1
        do itop=1, ctlint(C_NTOP)
          off=size(S_MINV+itop)-size(S_FIRV+itop)
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ik=1, size(S_K)
                tt=rotm(qvk(ik,Q_K),qvk(ikr,Q_K),1,itop)
     $            *rotm(qvk(ik,Q_K),qvk(ikc,Q_K),1,itop)
                do ivr=1, size(S_V+itop)
                  do ivc=1, size(S_V+itop)
                    do iv=1, size(S_MAXV+itop)
                      rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $              = rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $              + tt
     $              * tori(qvk(ik,Q_K),qvk(ikr,Q_K),iv,ivr+off
     $                     ,gamma(gam,itop),itop)
     $              * evalv(iv,gamma(gam,itop)
     $                     ,qvk(ik,Q_K),itop)
     $              * tori(qvk(ik,Q_K),qvk(ikc,Q_K),iv,ivc+off
     $                     ,gamma(gam,itop),itop)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

C     rotate evalv into rho system for each top without torsional integrals
C      if ((myand(ctlint(C_WOODS),1).eq.0).or.
C     $     (myand(ctlint(C_WOODS),4).eq.0).or.    ! don't use tors int. in rot. matrix
C     $     (myand(ctlint(C_WOODS),64).ne.0)) then ! use Demaisons method
      if (myand(ctlint(C_WOODS),4).eq.0) then
        itest=itest+1
        do itop=1, ctlint(C_NTOP)
          off=size(S_MINV+itop)-size(S_FIRV+itop)
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ik=1, size(S_K)
                tt=rotm(qvk(ik,Q_K),qvk(ikr,Q_K),1,itop)
     $            *rotm(qvk(ik,Q_K),qvk(ikc,Q_K),1,itop)
                do iv=1, size(S_V+itop)
                  rott(qvk(ikr,Q_K),qvk(ikc,Q_K),iv,iv,itop)
     $                 = rott(qvk(ikr,Q_K),qvk(ikc,Q_K),iv,iv,itop)
     $                 + tt
     $                 * evalv(iv+off,gamma(gam,itop),qvk(ik,Q_K),itop)
                end do
              end do
            end do
          end do
        end do
      end if
      if (itest.ne.1) stop 'Error: rotating in rotevl' 
      itest=0

C     multiply torsional integrals of one top like Demaison
      if (myand(ctlint(C_WOODS),64).ne.0) then 
        itest=itest+1
        do itop=1, ctlint(C_NTOP)
          off=size(S_MINV+itop)-size(S_FIRV+itop)
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $                 = rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $                 * tori(qvk(ikr,Q_K),qvk(ikc,Q_K)
     $                 ,ivr+off,ivc+off
     $                 ,gamma(gam,itop),itop)
                end do
              end do
            end do
          end do
        end do
      end if
      if (itest.gt.1) stop ' ERROR: tori two times multiplied '

      return
      end

C----------------------------------------------------------------------
      subroutine addo1(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ap)
C     add the rotated eigenvalues (in rott) to the main matrix h
C     real version
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR),ap
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     
      integer qv(DIMTOT) 
C     work
      real*8 rt,tt
      integer iv,itop,it,ir,ic,voff
      external myand
      integer myand
      
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in addo1'
      do iv=1,size(S_H)
        qv(iv)=int((iv-1)/size(S_K))+1
      end do

      do itop=1, ctlint(C_NTOP)
        voff=size(S_MINV+itop)-1
        if (myand(ctlint(C_PRI),AP_MH).ne.0)
     $       write(*,'(A,I2)') ' Add one-top-operator Top',itop
        do ir=1,size(S_H)
          do ic=1,ir
            if (myand(ctlint(C_WOODS),8).eq.0) then 
              tt=dble(ruse(qv(ir),qv(ic),itop)) ! don't mult. with tor. int. of the other tops  
            else                                ! use kronecker 
              tt=1.0
              do it=1,ctlint(C_NTOP)      ! supply tor.int. of the other tops
                if (it.ne.itop) tt=tt*tori(qvk(ir,Q_K)
     $               ,qvk(ic,Q_K)
     $               ,qvk(ir,Q_V+it)-size(S_MINV+it)+1
     $               ,qvk(ic,Q_V+it)-size(S_MINV+it)+1
     $               ,gamma(gam,it),it)
              end do
            end if
            rt=rott(qvk(ir,Q_K)
     $           ,qvk(ic,Q_K)
     $           ,qvk(ir,Q_V+itop)-voff
     $           ,qvk(ic,Q_V+itop)-voff
     $           ,itop)
     $           *tt
            h(ir,ic)=h(ir,ic)+rt*ap 
            if (myand(ctlint(C_PRI),AP_MH).ne.0) then
              if (abs(rt).lt.1000.0) then
                write(*,'(F10.5,$)') rt*ap
              else
                write(*,'(F10.2,$)') rt*ap
              end if
            end if
          end do
          if (myand(ctlint(C_PRI),AP_MH).ne.0) write(*,*)
        end do
      end do

c      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
c        write(*,*) ' H_ir '
c        do ir=1, size(S_H)
c          do ic=1, size(S_H)
c            if (abs(h(ir,ic)).lt.1000.0) then
c              write(*,'(F10.5,$)') h(ir,ic)
c            else
c              write(*,'(F10.2,$)') h(ir,ic)
c            end if
c          end do
c          write(*,*)
c        end do
c       write(*,*)
c     end if

      return
      end
C----------------------------------------------------------------------
      subroutine haddjp(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ipm)
      implicit none
      include 'iam.fi'
      integer j,gam,f,ipm
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     
      integer qv(DIMTOT) 
C     work
      real*8 rt,tt
      real*8 dcgam, dsgam
      integer iv,itop,ir,ic,voff
      external myand
      integer myand
            
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in ADDJP'
      do iv=1,size(S_H)
        qv(iv)=int((iv-1)/size(S_K))+1
      end do
      do itop=1, ctlint(C_NTOP)
        voff=size(S_MINV+itop)-1
C        if (ctlint(C_PRI).gt.11) write(*,'(/,A,I2)') 'H_Djp',itop
        do ir=1,size(S_H)
          do ic=1,ir
            tt=dble(ruse(qv(ir),qv(ic),itop)) ! don't mult. with tor. int. of the other tops  
            rt=rott(qvk(ir,Q_K)
     $           ,qvk(ic,Q_K)
     $           ,qvk(ir,Q_V+itop)!-voff
     $           ,qvk(ic,Q_V+itop)!-voff
     $           ,itop)
     $           *tt
     $           *2.0*(a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $                *dble(j*(j+1)))
            dcgam=cos(a(P1_GAMA+(itop-1)*DIMPIR)
     $           *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
            h(ir,ic)=h(ir,ic)+dcgam
            if (qvk(ir,Q_K).ne.qvk(ic,Q_K)) then
              dsgam=sin(a(P1_GAMA+(itop-1)*DIMPIR)
     $             *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
              h(ic,ir)=h(ic,ir)+dsgam
            end if
          end do
        end do
      end do
      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
        write(*,*) ' H_DJPM '
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            if (abs(h(ir,ic)).lt.1000.0) then
              write(*,'(F10.5,$)') h(ir,ic)
            else
              write(*,'(F10.2,$)') h(ir,ic)
            end if
          end do
          write(*,*)
        end do
        write(*,*)
      end if

      return
      end

C----------------------------------------------------------------------
      subroutine haddkp(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ipm)
      implicit none
      include 'iam.fi'
      integer j,gam,f,ipm
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     
      integer qv(DIMTOT) 
C     work
      real*8 rt,tt
      real*8 dcgam, dsgam
      integer iv,itop,ir,ic,voff
      external myand
      integer myand
      
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in HaddKP'
      do iv=1,size(S_H)
        qv(iv)=int((iv-1)/size(S_K))+1
      end do
      do itop=1, ctlint(C_NTOP)
        voff=size(S_MINV+itop)-1
C        if (ctlint(C_PRI).gt.11) write(*,'(/,A,I2)') 'H_DKp',itop
        do ir=1,size(S_H)
          do ic=1,ir
            tt=dble(ruse(qv(ir),qv(ic),itop)) ! don't mult. with tor. int. of the other tops  
            rt=rott(qvk(ir,Q_K)
     $           ,qvk(ic,Q_K)
     $           ,qvk(ir,Q_V+itop)!-voff
     $           ,qvk(ic,Q_V+itop)!-voff
     $           ,itop)
     $           *tt
     $           *(a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $           *(dble(qvk(ir,Q_K))**2+dble(qvk(ic,Q_K))**2))
C            h(ir,ic)=h(ir,ic)+rt
            dcgam=cos(a(P1_GAMA+(itop-1)*DIMPIR)
     $           *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
            h(ir,ic)=h(ir,ic)+dcgam
            if (qvk(ir,Q_K).ne.qvk(ic,Q_K)) then
              dsgam=sin(a(P1_GAMA+(itop-1)*DIMPIR)
     $             *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
              h(ic,ir)=h(ic,ir)+dsgam
            end if
          end do
        end do
      end do
      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
        write(*,*) ' H_DKP '
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            if (abs(h(ir,ic)).lt.1000.0) then
              write(*,'(F10.5,$)') h(ir,ic)
            else
              write(*,'(F10.2,$)') h(ir,ic)
            end if
          end do
          write(*,*)
        end do
        write(*,*)
      end if

      return
      end
C----------------------------------------------------------------------
      subroutine haddkkp(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ipm)
      implicit none!test to implement Dc3KK Dpi2KK
      include 'iam.fi'
      integer j,gam,f,ipm
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     
      integer qv(DIMTOT) 
C     work
      real*8 rt,tt
      real*8 dcgam, dsgam
      integer iv,itop,ir,ic,voff
      external myand
      integer myand
      
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in HaddKP'
      do iv=1,size(S_H)
        qv(iv)=int((iv-1)/size(S_K))+1
      end do
      do itop=1, ctlint(C_NTOP)
        voff=size(S_MINV+itop)-1
C        if (ctlint(C_PRI).gt.11) write(*,'(/,A,I2)') 'H_DKp',itop
        do ir=1,size(S_H)
          do ic=1,ir
            tt=dble(ruse(qv(ir),qv(ic),itop)) ! don't mult. with tor. int. of the other tops 
            rt=rott(qvk(ir,Q_K)
     $           ,qvk(ic,Q_K)
     $           ,qvk(ir,Q_V+itop)!-voff
     $           ,qvk(ic,Q_V+itop)!-voff
     $           ,itop)
     $           *tt
     $           *(a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $           *(dble(qvk(ir,Q_K))**4+dble(qvk(ic,Q_K))**4))
C            h(ir,ic)=h(ir,ic)+rt
            dcgam=cos(a(P1_GAMA+(itop-1)*DIMPIR)
     $           *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
            h(ir,ic)=h(ir,ic)+dcgam
            if (qvk(ir,Q_K).ne.qvk(ic,Q_K)) then
              dsgam=sin(a(P1_GAMA+(itop-1)*DIMPIR)
     $             *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))*rt
              h(ic,ir)=h(ic,ir)+dsgam
            end if
          end do
        end do
      end do
      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
        write(*,*) ' H_DKKP'
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            if (abs(h(ir,ic)).lt.1000.0) then
              write(*,'(F10.5,$)') h(ir,ic)
            else
              write(*,'(F10.2,$)') h(ir,ic)
            end if
          end do
          write(*,*)
        end do
        write(*,*)
      end if

      return
      end


C----------------------------------------------------------------------
      subroutine hadddp(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ipm)
      implicit none
      include 'iam.fi'
      integer j,gam,f,ipm
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     
      integer qv(DIMTOT) 
C     work
      real*8 tmp(DIM2J1,DIM2J1),djj1,dk,dff
      integer itop,ir,ic,voff,ivr,ivc,ikr,ikc
      external myand
      integer myand

      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in HAddDP'
      do ivr=1,size(S_H)
        qv(ivr)=int((ivr-1)/size(S_K))+1
      end do
      djj1=dble(j*(j+1))
      do itop=1, ctlint(C_NTOP)
        voff=size(S_MINV+itop)-1 
        do ivr=1, size(S_VV)
          do ivc=1, ivr
            if (ruse(ivr,ivc,itop).ne.0) then 
              do ikr=1, size(S_K)
                do ikc=1, size(S_K)
                  tmp(ikr,ikc)=0.0
                end do
              end do
              do ikc=1, size(S_K)-2
                ic=ikc+(ivc-1)*size(S_K)
                dk=dble(qvk(ic,Q_K)) 
                dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))
     $               *(djj1-(dk+1.0)*(dk+2.0)))
                do ikr=1, size(S_K)
                  ir=ikr+(ivr-1)*size(S_K)
                  tmp(ikr,ikc)=tmp(ikr,ikc)+
     $                 rott(qvk(ir,Q_K)     ,qvk(ic+2,Q_K)
     $                 ,qvk(ir,Q_V+itop),qvk(ic+2,Q_V+itop)
     $                 ,itop)
     $                 *dff
                  tmp(ikr,ikc+2)=tmp(ikr,ikc+2)+
     $                 rott(qvk(ir,Q_K)     ,qvk(ic,Q_K)
     $                 ,qvk(ir,Q_V+itop),qvk(ic,Q_V+itop)
     $                 ,itop)
     $                 *dff
                end do
              end do
              do ikr=1, size(S_K)-2
                ir=ikr+(ivr-1)*size(S_K)
                dk=dble(qvk(ir,Q_K)) 
                dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))
     $               *(djj1-(dk+1.0)*(dk+2.0)))
                do ikc=1, size(S_K)
                  ic=ikc+(ivc-1)*size(S_K)
                  tmp(ikr,ikc)=tmp(ikr,ikc)+
     $                 rott(qvk(ir+2,Q_K)     ,qvk(ic,Q_K)
     $                 ,qvk(ir+2,Q_V+itop),qvk(ic,Q_V+itop)
     $                 ,itop)
     $                 *dff
                  tmp(ikr+2,ikc)=tmp(ikr+2,ikc)+
     $                 rott(qvk(ir,Q_K)     ,qvk(ic,Q_K)
     $                 ,qvk(ir,Q_V+itop),qvk(ic,Q_V+itop)
     $                 ,itop)
     $                 *dff
                end do
              end do
              if (ivr.eq.ivc) then
                do ikr=1, size(S_K)
                  ir=ikr+(ivr-1)*size(S_K)
                  do ikc=1, ikr
                    ic=ikc+(ivc-1)*size(S_K)
                    h(ir,ic)=h(ir,ic)
     $                 +a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $                 *cos(a(P1_GAMA+(itop-1)*DIMPIR)
     $                 *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))
     $                 *tmp(ikr,ikc)
                    h(ic,ir)=h(ic,ir)
     $                 +a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $                 *sin(a(P1_GAMA+(itop-1)*DIMPIR)
     $                 *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))
     $                 *tmp(ikr,ikc)
                  end do
                enddo
              else
                do ikr=1, size(S_K)
                  ir=ikr+(ivr-1)*size(S_K)
                  do ikc=1, size(S_K)
                    ic=ikc+(ivc-1)*size(S_K)
                    h(ir,ic)=h(ir,ic)
     $                 +a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $                 *cos(a(P1_GAMA+(itop-1)*DIMPIR)
     $                 *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))
     $                 *tmp(ikr,ikc)
                    h(ic,ir)=h(ic,ir)
     $                 +a(DIMPRR+(itop-1)*DIMPIR+ipm)
     $                 *sin(a(P1_GAMA+(itop-1)*DIMPIR)
     $                 *dble(qvk(ir,Q_K)-qvk(ic,Q_K)))
     $                 *tmp(ikr,ikc)
                  end do
                end do
              end if
            end if
          end do
        end do
        if (myand(ctlint(C_PRI),AP_MH).ne.0) then
          write(*,'(A,I3)') ' H_DKD   top',itop
          do ir=1, size(S_H)
            do ic=1, size(S_H)
              if (abs(h(ir,ic)).lt.1000.0) then
                write(*,'(F10.5,$)') h(ir,ic)
              else
                write(*,'(F10.2,$)') h(ir,ic)
              end if
            end do
            write(*,*)
          end do
          write(*,*)
        end if
      end do

      return
      end

C----------------------------------------------------------------------
      subroutine haddo1(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,ap)
C     add the rotated matrixelements (in rott) to the main matrix h
C     hermitian version, uses gamma
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR),ap
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C
      integer qv(DIMTOT) 
C     work
      real*8 rt, dcgam, dsgam, tt, gam1
      integer iv,itop,it,ir,ic,qkr,qkc,vr1,vc1
      external myand
      integer myand
      
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in haddo1'
      do iv=1,size(S_H)
        qv(iv)=int((iv-1)/size(S_K))+1
      end do
      do itop=1, ctlint(C_NTOP)
        gam1=a(P1_GAMA+(itop-1)*DIMPIR)
        do ir=1,size(S_H)
          do ic=1,ir
            qkr=qvk(ir,Q_K)
            qkc=qvk(ic,Q_K)
            vr1=qvk(ir,Q_V+itop)-size(S_MINV+itop)+1
            vc1=qvk(ic,Q_V+itop)-size(S_MINV+itop)+1
            
            if (myand(ctlint(C_WOODS),8).eq.0) then 
              tt=dble(ruse(qv(ir),qv(ic),itop)) ! don't mult. with tor. int. of the other tops  
            else                                ! use kronecker 
              tt=1.0
              do it=1,ctlint(C_NTOP)      ! supply tor.int. of the other tops
                if (it.ne.itop) tt=tt*tori(qkr,qkc
     $               ,qvk(ir,Q_V+it)-size(S_MINV+it)+1
     $               ,qvk(ic,Q_V+it)-size(S_MINV+it)+1
     $               ,gamma(gam,it),it)
              end do
            end if
            rt=rott(qkr,qkc,vr1,vc1,itop)*tt*ap
            dcgam=cos(gam1*dble(qkr-qkc))
            dsgam=sin(gam1*dble(qkr-qkc))
            h(ir,ic)=h(ir,ic)+dcgam*rt ! real part of H
            h(ic,ir)=h(ic,ir)+dsgam*rt ! imag. part
            if (myand(ctlint(C_PRI),AP_MH).ne.0) then
              if (abs(dcgam*ap).lt.1000.0) then
                write(*,'(F10.5,$)') dcgam*ap
              else
                write(*,'(F10.2,$)') dcgam*ap
              end if
              if (abs(dsgam*ap).lt.1000.0) then
                write(*,'(F10.5,A,$)') dsgam*ap,'i'
              else
                write(*,'(F10.2,A,$)') dsgam*ap,'i'
              end if
            end if
          end do
          if (myand(ctlint(C_PRI),AP_MH).ne.0) write(*,*)
        end do
      end do
      
      return
      end

C ---------------------------------------------------------------------
      subroutine addiop(j,gam,f,qvk,ruse,h,a,ovv,rotm,rott,tori,ap,iop,
     $     offv)
      implicit none
      include 'iam.fi'
      integer j,gam,f,iop
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      logical offv
      real*8  h(DIMTOT,DIMTOT),ap
      real*8  a(DIMPAR)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      integer it1,it2,ivr,ivc,ik,ir,ic

      do it1=1, ctlint(C_NTOP)-1
        do it2=it1+1, ctlint(C_NTOP)
          do ivr=1, size(S_VV)
            do ivc=1, ivr
             if ((ivr.ne.ivc).or.offv) then
              do ik=1, size(S_K)
                ir=(ivr-1)*size(S_K)+ik
                ic=(ivc-1)*size(S_K)+ik
                h(ir,ic)=h(ir,ic)
     $               +ap
     $               *ovv(qvk(ir,Q_V+it1)
     $                   ,qvk(ic,Q_V+it1)
     $                   ,iop,gamma(gam,it1)
     $                   ,qvk(ir,Q_K),it1) 
     $               *ovv(qvk(ir,Q_V+it2)
     $                   ,qvk(ic,Q_V+it2)
     $                   ,iop,gamma(gam,it2)
     $                   ,qvk(ir,Q_K),it2) 
              end do
             end if
            end do
          end do
        end do
      end do
     
      return
      end

C ---------------------------------------------------------------------
      subroutine addovv(j,gam,f
     $     ,qvk,ruse,h,a,ovv,rotm,rott,tori,ap)
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT)
      real*8  a(DIMPAR),ap
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      integer it1,it2,ii,ir,ic,qkr,qkc,qki,vr1,vc1,vr2,vc2
      real*8  rt1,rt2
      external myand
      integer myand
      
      if (myand(ctlint(C_PRI),AP_MH).ne.0)
     $     write(*,'(/,A)') ' Add Ovv'
      do it1=1, ctlint(C_NTOP)-1
        do it2=it1+1, ctlint(C_NTOP)
          do ir=1, size(S_H)
            do ic=1, ir
              rt1=0.0
              rt2=0.0
              qkr=qvk(ir,Q_K)
              qkc=qvk(ic,Q_K)
              vr1=qvk(ir,Q_V+it1)-size(S_MINV+it1)+1
              vc1=qvk(ic,Q_V+it1)-size(S_MINV+it1)+1
              vr2=qvk(ir,Q_V+it2)-size(S_MINV+it2)+1
              vc2=qvk(ic,Q_V+it2)-size(S_MINV+it2)+1
              do ii=1, size(S_K)
                qki=qvk(ii,Q_K)
                rt1=rt1  
     $               + rott(qkr,qki,vr1,vc1,it1)
     $               * rott(qki,qkc,vr2,vc2,it2)
                rt2=rt2
     $               + rott(qkr,qki,vr2,vc2,it2)
     $               * rott(qki,qkc,vr1,vc1,it1)
              end do
              h(ir,ic)=h(ir,ic)+ap*(rt1+rt2)
              if (myand(ctlint(C_PRI),AP_MH).ne.0) then
                if (abs(ap*(rt1+rt2)).lt.1000.0) then
                  write(*,'(F10.5,$)') ap*(rt1+rt2)
                else
                  write(*,'(F10.2,$)') ap*(rt1+rt2)
                end if
              end if
            end do
            if (myand(ctlint(C_PRI),AP_MH).ne.0) write(*,*)
          end do
        end do
      end do
      return
      end

C ---------------------------------------------------------------------
      subroutine haddovv(j,gam,f
     $     ,qvk,ruse,h,a,ovv,rotm,rott,tori,ap)
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT)
      real*8  a(DIMPAR),ap
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      integer it1,it2,ii,ir,ic,qkr,qkc,qki,vr1,vc1,vr2,vc2
      real*8  tr1,tr2,ti1,ti2,gam1,gam2,t1,t2
      external myand
      integer myand
      
      if (myand(ctlint(C_PRI),AP_MH).ne.0)
     $     write(*,'(A)') ' H AddOvv' 
      do it1=1, ctlint(C_NTOP)-1
        do it2=it1+1, ctlint(C_NTOP)
          gam1=a(P1_GAMA+(it1-1)*DIMPIR)
          gam2=a(P1_GAMA+(it2-1)*DIMPIR)
          do ir=1, size(S_H)
            do ic=1, ir
              tr1=0.0
              tr2=0.0
              ti1=0.0
              ti2=0.0
              qkr=qvk(ir,Q_K)
              qkc=qvk(ic,Q_K)
              vr1=qvk(ir,Q_V+it1)-size(S_MINV+it1)+1
              vc1=qvk(ic,Q_V+it1)-size(S_MINV+it1)+1
              vr2=qvk(ir,Q_V+it2)-size(S_MINV+it2)+1
              vc2=qvk(ic,Q_V+it2)-size(S_MINV+it2)+1
              do ii=1, size(S_K)
                qki=qvk(ii,Q_K)
                t1=
     $               rott(qkr,qki,vr1,vc1,it1)
     $               * rott(qki,qkc,vr2,vc2,it2)
                tr1=tr1+t1
     $               * cos(gam1*dble(qkr-qki)
     $                    +gam2*dble(qki-qkc))
                ti1=ti1+t1
     $               * sin(gam1*dble(qkr-qki)
     $                    +gam2*dble(qki-qkc))
                t2=
     $               rott(qkr,qki,vr2,vc2,it2)
     $               * rott(qki,qkc,vr1,vc1,it1)
                tr2=tr2+t2
     $               * cos(gam2*dble(qkr-qki)
     $                    +gam1*dble(qki-qkc))
                ti2=ti2+t2
     $               * sin(gam2*dble(qkr-qki)
     $                    +gam1*dble(qki-qkc))
              end do
              h(ir,ic)=h(ir,ic)+ap*(tr1+tr2)
              h(ic,ir)=h(ic,ir)+ap*(ti1+ti2)
              if (myand(ctlint(C_PRI),AP_MH).ne.0)
     $             write(*,'(2F11.5,A,$)')
     $             ap*(tr1+tr2),ap*(ti1+ti2),'i'
            end do
            if (myand(ctlint(C_PRI),AP_MH).ne.0) write(*,*)
          end do
        end do
      end do
      return
      end

C ---------------------------------------------------------------------
      subroutine rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,ifs1,ifs2,sc)
      implicit none
      include 'iam.fi'
      integer j,gam,f,ifs1,ifs2,ifs
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      logical sc
      real*8  a(DIMPAR)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      integer itop,ikr,ikc,ivr,ivc,ik,iv,off
      real*8  tt
      real*8  rtmp(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV)
      real*8  ovt(DIMV,DIMV,DIM2J1)
      external myand
      integer myand

      do itop=1, ctlint(C_NTOP)
        do ikr=1, size(S_K)
          do ikc=1, size(S_K)
            do ivr=1, size(S_V+itop)
              do ivc=1, size(S_V+itop)
                rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)=0.0
              end do
            end do
          end do
        end do
      end do

      if (myand(ctlint(C_WOODS),16).ne.0) then
        do itop=1, ctlint(C_NTOP)
          if (itop.eq.1) then
            ifs=ifs1
          else
            ifs=ifs2
          endif
          off=size(S_MINV+itop)-size(S_FIRV+itop)
          if (sc) then
            do ik=1, size(S_K)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  ovt(ivr,ivc,ik)=
     $                 ovv(ivr+off,ivc+off,ifs,gamma(gam,itop)
     $                 ,qvk(ik,Q_K),itop)
     $                 -ovv(ivr+off,ivc+off,ifs,gamma(1,itop)
     $                 ,qvk(ik,Q_K),itop)
                end do
              end do
            end do
          else
            do ik=1, size(S_K)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  ovt(ivr,ivc,ik)=
     $                 ovv(ivr+off,ivc+off,ifs,gamma(gam,itop)
     $                 ,qvk(ik,Q_K),itop)
                end do
              end do
            end do
          end if
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  rtmp(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc)=0.0
                end do
              end do
            end do
          end do
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  do iv =1, size(S_MAXV+itop)
                    rtmp(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc)=
     $                   rtmp(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc)
     $                   + ovt(ivr,iv,ikr)
     $                   * tori(qvk(ikr,Q_K),qvk(ikc,Q_K),iv,ivc+off
     $                          ,gamma(gam,itop),itop)
     $                   * rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,itop)
                  end do
                end do
              end do
            end do
          end do
          do ikr=1, size(S_K)
            do ikc=1, size(S_K)
              do ik=1, size(S_K)
                do ivr=1, size(S_V+itop)
                  do ivc=1, size(S_V+itop)
                    do iv =1, size(S_MAXV+itop)
                      rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)=
     $                     rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $                   + tori(qvk(ik ,Q_K),qvk(ikr,Q_K),iv ,ivr+off
     $                          ,gamma(gam,itop),itop)
     $                   * rotm(qvk(ik ,Q_K),qvk(ikr,Q_K),1,itop)
     $                   * rtmp(qvk(ik ,Q_K),qvk(ikc,Q_K),iv ,ivc)
                    end do
                  end do
                end do
              end do
            end do
          end do
        end do
      end if

C     rotate ifs into rho system for each top without torsional integrals
      if (myand(ctlint(C_WOODS),16).eq.0) then

      do itop=1, ctlint(C_NTOP)
        if (itop.eq.1) then
          ifs=ifs1
        else
          ifs=ifs2
        endif
        off=size(S_MINV+itop)-size(S_FIRV+itop)

        if (sc) then
          do ik=1, size(S_K)
            do ivr=1, size(S_V+itop)
              do ivc=1, size(S_V+itop)
                ovt(ivr,ivc,ik)=
     $               ovv(ivr+off,ivc+off,ifs,gamma(gam,itop)
     $               ,qvk(ik,Q_K),itop)
     $              -ovv(ivr+off,ivc+off,ifs,gamma(1,itop)
     $               ,qvk(ik,Q_K),itop)
              end do
            end do
          end do
        else
          do ik=1, size(S_K)
            do ivr=1, size(S_V+itop)
              do ivc=1, size(S_V+itop)
                ovt(ivr,ivc,ik)=
     $               ovv(ivr+off,ivc+off,ifs,gamma(gam,itop)
     $               ,qvk(ik,Q_K),itop)
              end do
            end do
          end do
        end if

        do ikr=1, size(S_K)
          do ikc=1, size(S_K)
            do ik=1, size(S_K)
              tt=rotm(qvk(ik,Q_K),qvk(ikr,Q_K),1,itop)
     $          *rotm(qvk(ik,Q_K),qvk(ikc,Q_K),1,itop)
              do ivr=1, size(S_V+itop)
                do ivc=1, size(S_V+itop)
                  rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $          = rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
     $                   + tt
     $                   * ovt(ivr,ivc,ik)
                end do
              end do
            end do
          end do
        end do
      end do
      end if

      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
      do itop=1, ctlint(C_NTOP)
        write(*,'(/,A)')
     $       'Matrixelements after rotation into principal axes system'
        write(*,'(2(A,I2))')
     $       '  K  V | Operator No.',ifs,'    Top No.',itop
        do ivr=1, size(S_V+itop)
          do ikr=1, size(S_K)
            write(*,'(2I3,$)') qvk(ikr,Q_K),ivr
            do ivc=1, size(S_V+itop)
              do ikc=1, size(S_K)
                write(*,'(F10.4,$)')
     $               rott(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc,itop)
              end do
            end do
            write(*,*)
          end do
        end do
      end do
      end if

      return
      end
C =========================================
C ---------------------------------------------------------------------
      subroutine prrott(j,gam,f,qvk,a,ovv,rotm,rott,tori,ifs1,ifs2,sc)
      implicit none
      include 'iam.fi'
      integer j,gam,f,ifs1,ifs2
      integer it,ir,ic,it1
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      logical sc
      real*8  a(DIMPAR)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C work
      real*8 tt,hh

      call rotovv(j,gam,f,qvk,a,ovv,rotm,rott,tori,ifs1,ifs2,.false.)
      write(*,*) ' rott '
      do it=1,ctlint(C_NTOP)
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            tt=1.0
            do it1=1,ctlint(C_NTOP)
              if (it.ne.it1)
     $             tt=tt
     $             * tori(qvk(ir,Q_K)
     $             ,qvk(ic,Q_K)
     $             ,qvk(ir,Q_V+it1)
     $             ,qvk(ic,Q_V+it1)
     $             ,gamma(gam,it),it1)
            end do
            hh=rott(qvk(ir,Q_K)
     $           ,qvk(ic,Q_K)
     $           ,qvk(ir,Q_V+it)
     $           ,qvk(ic,Q_V+it)
     $           ,it)
     $           * tt
            if (abs(hh).lt.1000.0) then
              write(*,'(F10.5,$)') hh
            else
              write(*,'(F10.2,$)') hh
            end if
          end do
          write(*,*)
        end do
        write(*,*)
      end do
      return
      end
