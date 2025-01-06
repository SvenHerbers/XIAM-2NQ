C                   ! I (Sven) made some changes to the intensity calculations for Xiam_2023
C                   ! The intensities from XIAM didnt seem to make a lot of sense to me?
C                   ! The column log10_total gives now populationdifference * reduced dipole^2
C                   ! While the column SPCAT should reproduce the logaritmic intensities from SPCAT
C                   ! My current understanding is that log10_total might be better for coherence experiments
C                   ! When efield amplitudes are recorded
C                   ! Intensities in spcat might be better for absorption experiments. But I am not sure.
C                   !
C                    ! Picketts calcat.c       uses  
C                    !          
C                    ! str = dgn * strr * fac * frq * (1. - exp(tmq * frq));
C                    ! strlg = log10(str + zero) + tmql * elow;
C                    ! 
C                    ! with strr being the transition dp**2
C                    ! dgn being - no idea but not the degeneracy of the lower level.
C                    ! strlg is the intensity given in cat
C                    ! elow is energy of lower state
C                    ! fac = (4.16231e-5) / qrot
C                    ! frq = frequency
C                    ! tmq = tmc / (Temperature * cmc)
C                    ! tmql = tmq * 0.43429448
C                    !
C                    !Params defined above:
C                    ! tmc = -1.43878
C                    ! cmc = 29979.2458
C                    !
C                    ! intensities following picketts definition will be written as 
C                    ! spcat (name in code "spcat"
C                    !
C                    ! intensities following Populationdifference*DP**2 will be denoted
C                    ! as total (name in code "bolt") 
C----------------------------------------------------------------------
      subroutine intens(ix,ints,mux,muy,muz,tori)
      implicit none
      integer ix
      include 'iam.fi'
      real*8  ints,muz, mux, muy
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer sizeup(DIMSIZ),sizelo(DIMSIZ)
      integer qvkup(DIMTOT,Q_K:Q_V+DIMTOP),qvklo(DIMTOT,Q_K:Q_V+DIMTOP)
      integer qup(DIMQLP),qlo(DIMQLP) 
      real*8  vrup(DIMTOT),vrlo(DIMTOT),viup(DIMTOT),vilo(DIMTOT)
C
      real*8 eup,elo
      real*8  eloM(DIMTOT),eupM(DIMTOT)! for matrix reading
      real*8 sj,sj2
      real*8 sjsum
      real*8 sj2sum
      integer F1slo(DIMQ2)
C      real*8 nF1slo(DIMQ2) ! here only DIMQ, because only vector is read not matrix 
      integer F1sup(DIMQ2)
C      real*8 nF1sup(DIMQ2)
      real*8 nF1slo(DIMTOT,DIMQ2) ! I swaped to full matrix reading
      real*8 nF1sup(DIMTOT,DIMQ2)
      real*8  hrup(DIMTOT,DIMTOT),hrlo(DIMTOT,DIMTOT) !needed in full matrix reading
      real*8  hiup(DIMTOT,DIMTOT),hilo(DIMTOT,DIMTOT) !needed in full matrix reading
      
      integer jup,jlo,fup,flo,tup,tlo,sup,slo,bup,blo
      integer f1lo,f1up !h2024
      integer i,j ! do loop inticies
      
      jup=qlin(ix,Q_J,Q_UP)
      tup=qlin(ix,Q_T,Q_UP)
      sup=qlin(ix,Q_S,Q_UP)
      fup=qlin(ix,Q_F,Q_UP)
      f1up=qlin(ix,Q_F1,Q_UP)!h2024
      bup=qlin(ix,Q_B,Q_UP)

      jlo=qlin(ix,Q_J,Q_LO)
      tlo=qlin(ix,Q_T,Q_LO)
      slo=qlin(ix,Q_S,Q_LO)
      flo=qlin(ix,Q_F,Q_LO)
      f1lo=qlin(ix,Q_F1,Q_LO)!h2024
      blo=qlin(ix,Q_B,Q_LO)

      if (sup.ne.slo) then
        ints=-99.0
        return
      end if
      if (abs(jup-jlo).gt.1) then
        ints=-98.0
        return
      end if
      
          
      if (jlo.gt.jup) then
        call rdvec(jlo,tlo,slo,flo,blo,f1lo,sizeup
     $            ,qvkup,vrup,viup,eup,qup)
        call rdvec(jup,tup,sup,fup,bup,f1up,sizelo
     $            ,qvklo,vrlo,vilo,elo,qlo)
        call calir(ints,tori,jlo,jup,slo,muz, mux, muy
     $     ,sizeup,qvkup,vrup,viup
     $     ,sizelo,qvklo,vrlo,vilo)
      end if                                   
      if (jup.gt.jlo) then                     
        call rdvec(jup,tup,sup,fup,bup,f1up,sizeup
     $            ,qvkup,vrup,viup,eup,qup)
        call rdvec(jlo,tlo,slo,flo,blo,f1lo,sizelo
     $            ,qvklo,vrlo,vilo,elo,qlo)
        call calir(ints,tori,jup,jlo,slo,muz, mux, muy
     $     ,sizeup,qvkup,vrup,viup
     $     ,sizelo,qvklo,vrlo,vilo)
      end if                                   
      if (jup.eq.jlo) then                     
        call rdvec(jup,tup,sup,fup,bup,f1up,sizeup
     $            ,qvkup,vrup,viup,eup,qup)
        call rdvec(jlo,tlo,slo,flo,blo,f1lo,sizelo
     $            ,qvklo,vrlo,vilo,elo,qlo)
        call caliq(ints,tori,jup,jlo,slo,muz, mux, muy
     $     ,sizeup,qvkup,vrup,viup
     $     ,sizelo,qvklo,vrlo,vilo)
      end if


      if ((ctlint(C_SPIN2).ne.0).and.(ctlint(C_SPIN).ne.0)
     $                                   .and.(flo.ne.-1)) then ! approximate to nuclei intensities, by weighting the wigner symbol prorducts by the contribution of the various F1s.
C     There seems to be problems with rdvecQ... I replaced it with rdmatQ for now. - Sven 04.01.2025
C        call rdvecQ(jlo,tlo,slo,flo,blo,f1lo,sizelo
C     $            ,qvklo,vrlo,vilo,elo,qlo,F1slo,nF1slo)
C        call rdvecQ(jup,tup,sup,fup,bup,f1up,sizeup
C     $            ,qvkup,vrup,viup,eup,qup,F1sup,nF1sup)
        call rdmatQ(jup,sup,fup,bup,f1up,sizeup
     $            ,qvkup,hrup,hiup,eupM,qup,F1sup,nF1sup)
        call rdmatQ(jlo,slo,flo,blo,f1lo
     $             ,sizelo,qvklo,hrlo,hilo,eloM,qlo,F1slo,nF1slo)
     
      sjsum=0.0
      sj2sum=0.0
      do i=1,DIMQ2 
        do j=1, DIMQ2 
          if ((F1sup(i).ne.(-1)).and.(F1slo(j).ne.(-1))) then
           call sixj(2*jup,F1sup(i),ctlint(C_SPIN),F1slo(j),2*jlo,2,sj)
           call sixj(F1sup(i),fup,ctlint(C_SPIN2),flo,F1slo(j),2,sj2)
C           sj=sj**2*(F1sup(i)+1)*(F1slo(j)+1)
C           sj2=sj2**2*(fup+1)*(flo+1)
C           sjsum=sjsum+sj*sj2*(nF1slo(j)*nF1sup(i))**2
           sj=sj*((F1sup(i)+1)*(F1slo(j)+1))**0.5
           sj2=sj2*((fup+1)*(flo+1))**0.5
           sjsum=sjsum+abs(sj*sj2*(nF1slo(tlo,j)*nF1sup(tup,i))) !added tlo tup indicies for matrix reading.
          end if
        end do
      end do
C      ints=ints*sjsum
       ints=ints*sjsum**2
      else
       if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
         if ((ctlint(C_SPIN2).eq.0).or.(flo.eq.-1)) then
          sj2=1.0
         else
          call sixj(f1up,fup,ctlint(C_SPIN2),flo,f1lo,2,sj2)
          sj2=sj2**2
          sj2=abs(sj2*(fup+1)*(flo+1))
         end if
         if (f1lo.eq.-1) then
          sj=1.0
         else
          call sixj(2*jup,f1up,ctlint(C_SPIN),f1lo,2*jlo,2,sj)
          sj=sj**2
          sj=abs(sj*(f1up+1)*(f1lo+1))!/(ctlint(C_SPIN)+1)) !Pickett does not use devision by 2I+1
         end if
         ints=ints*sj*sj2
       end if
      end if

      return
      end

C ----------------------------------------------------------------------
      subroutine caliq(ints,tori,jup,jlo,gam,muz, mux, muy
     $     ,sizeup,qvkup,vrup,viup
     $     ,sizelo,qvklo,vrlo,vilo)
      implicit none
      include 'iam.fi'
      integer jup,jlo,gam
      real*8 ints,muz, mux, muy
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer sizeup(DIMSIZ),sizelo(DIMSIZ)
      integer qvkup(DIMTOT,Q_K:Q_V+DIMTOP),qvklo(DIMTOT,Q_K:Q_V+DIMTOP)
      real*8  vrup(DIMTOT),vrlo(DIMTOT),viup(DIMTOT),vilo(DIMTOT)
C
      real*8  intr,inti,dq,dj,tt
      integer ik,ik1,ik2,qk,qk1,qk2,iv,iv1,iv2,i,i1,i2
      integer itop,ntop
     
      if (jlo.ne.jup) then  !Herbers2024 warning removal
      write(0,*) 'CALIQ: intens error: jlo <> jup'
      end if
      ntop=ctlint(C_NTOP)
      inti=0.0
      intr=0.0
      dj=dble(jlo)
C     muz
      do ik=1, sizelo(S_K)
        do iv=1, sizelo(S_VV)
          i=ik+(iv-1)*sizelo(S_K)
          qk=qvklo(i,Q_K)
          dq=dble(qk)
          intr=intr+2.0*muz*dq*(vrlo(i)*vrup(i)+vilo(i)*viup(i))
          inti=inti+2.0*muz*dq*(vrlo(i)*viup(i)-vilo(i)*vrup(i))
        end do
      end do
C     mux/y
      do ik2=1, sizelo(S_K)-1
        qk2=qvkup(ik2,Q_K)
        ik1=ik2+1
        qk1=qvklo(ik1,Q_K)
        dq=sqrt((dj-dble(qk2))*(dj+dble(qk2)+1.0))
        do iv1=1, sizelo(S_VV)
          i1=ik1+(iv1-1)*sizelo(S_K)
          do iv2=1, sizeup(S_VV)
            i2=ik2+(iv2-1)*sizeup(S_K)
            tt=1.0
            if (gam.ne.0) then
              do itop=1, ntop
                tt=tt*tori(qk1,qk2
     $               ,qvklo(i1,Q_V+itop)-size(S_MINV+itop)+1
     $               ,qvkup(i2,Q_V+itop)-size(S_MINV+itop)+1
     $               ,gamma(gam,itop),itop)
              end do
            else
              if (iv1.ne.iv2) tt=0.0
            end if
            intr=intr
     $           -tt*mux*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            inti=inti
     $           -tt*mux*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
            inti=inti
     $           -tt*muy*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            intr=intr
     $           +tt*muy*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
          end do
        end do
      end do
      do ik2=2, sizelo(S_K)
        qk2=qvkup(ik2,Q_K)
        ik1=ik2-1
        qk1=qvklo(ik1,Q_K)
        dq=sqrt((dj+dble(qk2))*(dj-dble(qk2)+1.0))
        do iv1=1, sizelo(S_VV)
          i1=ik1+(iv1-1)*sizelo(S_K)
          do iv2=1, sizeup(S_VV)
            i2=ik2+(iv2-1)*sizeup(S_K)
            tt=1.0
            if (gam.ne.0) then
              do itop=1, ntop
                tt=tt*tori(qk1,qk2
     $               ,qvklo(i1,Q_V+itop)-size(S_MINV+itop)+1
     $               ,qvkup(i2,Q_V+itop)-size(S_MINV+itop)+1
     $               ,gamma(gam,itop),itop)
              end do
            else
              if (iv1.ne.iv2) tt=0.0
            end if
            intr=intr
     $           -tt*mux*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            inti=inti
     $           -tt*mux*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
            inti=inti
     $           +tt*muy*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            intr=intr
     $           -tt*muy*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
          end do
        end do
      end do
      ints=(intr**2+inti**2)*(2.0*dj+1.0)/(2.0*dj*(2.0*dj+2.0))
      return
      end

C ----------------------------------------------------------------------
C     Delta J = 1 ----------------------------------------------------------------
      subroutine calir(ints,tori,jup,jlo,gam,muz, mux, muy
     $     ,sizeup,qvkup,vrup,viup
     $     ,sizelo,qvklo,vrlo,vilo)
      implicit none
      include 'iam.fi'
      integer jup,jlo,gam 
      real*8  ints,muz, mux, muy
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer sizeup(DIMSIZ),sizelo(DIMSIZ)
      integer qvkup(DIMTOT,Q_K:Q_V+DIMTOP),qvklo(DIMTOT,Q_K:Q_V+DIMTOP)
      real*8  vrup(DIMTOT),vrlo(DIMTOT),viup(DIMTOT),vilo(DIMTOT)
C
      real*8  intr,inti,dq,dj,tt
      integer ik,ik1,ik2,qk,qk1,qk2,iv,iv1,iv2,i,i1,i2
      integer itop,ntop
     
      if ((jlo+1).ne.jup) then !Herbers2024 warning removal
      write(0,*) 'CALIR: intens error: jlo+1 <> jup'
      end if
      ntop=ctlint(C_NTOP)
      inti=0.0
      intr=0.0
      dj=dble(jlo)
C     mu_z
      sizelo(S_K)=2*jlo+1
      do ik2=2, sizelo(S_K)+1 !z types.
        qk2=qvkup(ik2,Q_K)
        
        ik1=ik2-1 
        qk1=qvklo(ik1,Q_K)
        do iv=1, sizelo(S_VV)
          i1=ik1+(iv-1)*sizelo(S_K)
          i2=ik2+(iv-1)*sizeup(S_K)
c          qk1=qvklo(i1,Q_K)
c          qk2=qvkup(i2,Q_K)
c         if (qk1.eq.qk2) write(0,*) 'Matching K numbers        '
c    $    , qk1, qk2, jlo, jup  
          if (qk1.ne.qk2) write(0,*) 'Error in CALIR : K numbers'
     $    , qk1, qk2, jlo, jup 
          dq=2.0*sqrt((dj+dble(qk1)+1.0)*(dj-dble(qk1)+1.0))
          intr=intr-muz*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
          inti=inti-muz*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
        end do
      end do
C     mux/y 
      do ik2=1, sizelo(S_K) ! x/y types.
        qk2=qvkup(ik2,Q_K)
        ik1=ik2  
        qk1=qvklo(ik1,Q_K)
        
        dq=sqrt((dj-dble(qk2))*(dj-dble(qk2)+1.0))
        do iv1=1, sizelo(S_VV)
          i1=ik1+(iv1-1)*sizelo(S_K)
          do iv2=1, sizeup(S_VV)
            i2=ik2+(iv2-1)*sizeup(S_K)
            tt=1.0
            if (gam.ne.0) then
              do itop=1, ntop
                tt=tt*tori(qk1,qk2
     $               ,qvklo(i1,Q_V+itop)-size(S_MINV+itop)+1
     $               ,qvkup(i2,Q_V+itop)-size(S_MINV+itop)+1
     $               ,gamma(gam,itop),itop)
              end do
            else
              if (iv1.ne.iv2) tt=0.0
            end if
            intr=intr
     $           +tt*mux*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            inti=inti
     $           +tt*mux*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
            inti=inti
     $           +tt*muy*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            intr=intr
     $           -tt*muy*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
          end do
        end do
      end do
      do ik2=3, sizelo(S_K)+2
        qk2=qvkup(ik2,Q_K)
        ik1=ik2-2
        qk1=qvklo(ik1,Q_K)
        dq=sqrt((dj+dble(qk2))*(dj+dble(qk2)+1.0))
        do iv1=1, sizelo(S_VV)
          i1=ik1+(iv1-1)*sizelo(S_K)
          do iv2=1, sizeup(S_VV)
            i2=ik2+(iv2-1)*sizeup(S_K)
            tt=1.0
            if (gam.ne.0) then
              do itop=1, ntop
                tt=tt*tori(qk1,qk2
     $               ,qvklo(i1,Q_V+itop)-size(S_MINV+itop)+1
     $               ,qvkup(i2,Q_V+itop)-size(S_MINV+itop)+1
     $               ,gamma(gam,itop),itop)
              end do
            else
              if (iv1.ne.iv2) tt=0.0
            end if
            intr=intr
     $           -tt*mux*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            inti=inti
     $           -tt*mux*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
            inti=inti
     $           +tt*muy*dq*(vrlo(i1)*vrup(i2)+vilo(i1)*viup(i2))
            intr=intr
     $           -tt*muy*dq*(vrlo(i1)*viup(i2)-vilo(i1)*vrup(i2))
          end do
        end do
      end do
      ints=(intr**2+inti**2)/(2.0*(2.0*dj+2.0))        
      return
      end

C ----------------------------------------------------------------------
      subroutine binnam(ij,gam,if,ib,if1,binfname)
      implicit none
      include 'iam.fi'
      integer ij,if,gam,ib,if1 !h2024
      character*(*) binfname     
C      write(0,*) if1, int(if1/10), mod(if1,10)
      if (if.ge.0) then
       if (if1.ge.0) then 
        write(binfname,'(A,2I1,A,2I1,A,2I1,A,I1,A,2I1)')!h2024
     $       'j',int(ij/10),mod(ij,10),           !h2024
     $       'f',int(if/10),mod(if,10),           !h2024
     $       'f1',int(if1/10),mod(if1,10),        !h2024
     $       'b',ib,                              !h2024
     $       '.s',int(gam/10),mod(gam,10)         !h2024
       else
              write(binfname,'(A,2I1,A,2I1,A,I1,A,2I1)')
     $       'j',int(ij/10),mod(ij,10),
     $       'f',int(if/10),mod(if,10),
     $       'b',ib,
     $       '.s',int(gam/10),mod(gam,10)
       end if
      else
        write(binfname,'(A,2I1,A,I1,A,2I1)')
     $       'j',int(ij/10),mod(ij,10),
     $       'b',ib,
     $       '.s',int(gam/10),mod(gam,10)
      end if
      return
      end

C ----------------------------------------------------------------------
      subroutine wrvec(zr,zi,d,ij,gam,if,ib,if1,qvk,qmk)
      implicit none
      include 'iam.fi'
      integer ij,if,gam,ib,if1
      real*8 zr(DIMTOT,DIMTOT),zi(DIMTOT,DIMTOT)
      real*8 d(DIMTOT)
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP),qmk(DIMTOT,DIMQLP)
      integer ie,id,itop,bu,ios
      character*30 binfname     
      integer myand
      external myand

      call binnam(ij,gam,if,ib,if1,binfname)
      call getfu(bu)
      if (myand(ctlint(C_PRI),AP_ST).ne.0)
     $     write(*,*)' Writing ',binfname
      
      open(bu,file=binfname,status='unknown',err=99,iostat=ios
     $     ,form='unformatted')
      write(bu,err=99,iostat=ios) (size(ie),ie=1,DIMSIZ)
      write(bu,err=99,iostat=ios) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,size(S_H))
      do id=1, size(S_H)
c        write(bu,'(F18.9,$)',err=99,iostat=ios) d(id)
c        write(bu,err=99,iostat=ios) d(id)
        write(bu,err=99,iostat=ios)  d(id), (qmk(id,ie),ie=1,DIMQLP),
     $       (zr(ie,id),zi(ie,id),ie=1,size(S_H))
      end do
      close (bu)
      return
 99   continue
      write(*,*) 'wrvec Error iostat', ios,ij,gam,if,ib,if1,binfname
      return
      end

C----------------------------------------------------------------------
      subroutine getfu(myunit)
      implicit none
      integer myunit,funit
      data funit/10/
      save funit

      if (funit.gt.50) funit=10
      funit=funit+1
      myunit=funit
c      write(0,*)'funit ',myunit
      return
      end

C----------------------------------------------------------------------
      subroutine rdvec(ij,it,is,if,ib,if1,mysize,qvk,vr,vi,e,qmk)
      implicit none
      include 'iam.fi'
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer mysize(DIMSIZ)
      integer ij,it,is,if,ib,i,ie,itop,bu,if1
      integer qmk(DIMQLP) 
      real*8  vr(DIMTOT),vi(DIMTOT),e
      character*30 binfname     

      call getfu(bu)
c      write(0,*)'rdvec',bu
      call binnam(ij,is,if,ib,if1,binfname)
      open(bu,file=binfname,status='unknown',form='unformatted')
      read(bu) (mysize(ie),ie=1,DIMSIZ)
      read(bu) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,mysize(S_H))
      do i=1, mysize(S_H)
        read(bu) e, (qmk(ie),ie=1,DIMQLP),
     $       (vr(ie),vi(ie),ie=1,mysize(S_H))
        if (qmk(Q_T).eq.it) goto 10
      end do
      write(0,*) 'eigenvector not found' !wraning removal Herbers2024
 10   continue
      close (bu)
      return
      end

C ----------------------------------------------------------------------
      subroutine rdmat(ij,is,if,ib,if1,mysize,qvk,vr,vi,e,qmk)
      implicit none
      include 'iam.fi'
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer mysize(DIMSIZ)
      integer ij,is,if,ib,i,ie,itop,bu,if1
      integer qmk(DIMTOT,DIMQLP) 
      real*8  vr(DIMTOT,DIMTOT),vi(DIMTOT,DIMTOT),e(DIMTOT)
      character*30 binfname     

      call getfu(bu)
c      write(0,*)'rdmat',bu
      call binnam(ij,is,if,ib,if1,binfname)
      open(bu,file=binfname,status='unknown',form='unformatted')
      read(bu) (mysize(ie),ie=1,DIMSIZ)
      read(bu) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,mysize(S_H))
      do i=1, mysize(S_H)
        read(bu) e(i), (qmk(i,ie),ie=1,DIMQLP),
     $       (vr(ie,i),vi(ie,i),ie=1,mysize(S_H))
      end do
      close (bu)
      return
      end
C----------------------------------------------------------------------
C ----------------------------------------------------------------------
      subroutine wrvecQ(zr,zi,d,ij,gam,if,ib,if1,qvk,qmk,F1s,normis)
      implicit none
      include 'iam.fi'
      integer ij,if,gam,ib,if1
      real*8 zr(DIMTOT,DIMTOT),zi(DIMTOT,DIMTOT)
      real*8 d(DIMTOT)
      real*8 normis(DIMTOT,DIMQ2) ! stores the eigenvector**2 for each F1
      integer F1s(DIMQ2) ! stores the F1 quantum number contributions
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP),qmk(DIMTOT,DIMQLP)
      integer ie,id,itop,bu,ios
      character*30 binfname     
      integer myand
      external myand

      call binnam(ij,gam,if,ib,if1,binfname)
      call getfu(bu)
      if (myand(ctlint(C_PRI),AP_ST).ne.0)
     $     write(*,*)' Writing ',binfname
      
      open(bu,file=binfname,status='unknown',err=99,iostat=ios
     $     ,form='unformatted')
      write(bu,err=99,iostat=ios) (size(ie),ie=1,DIMSIZ)
      write(bu,err=99,iostat=ios) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,size(S_H))
      do id=1, size(S_H)
c        write(bu,'(F18.9,$)',err=99,iostat=ios) d(id)
c        write(bu,err=99,iostat=ios) d(id)
        write(bu,err=99,iostat=ios)  d(id), (qmk(id,ie),ie=1,DIMQLP),
     $       (zr(ie,id),zi(ie,id),ie=1,size(S_H))
      end do
      write(bu,err=99,iostat=ios) (F1s(ie),ie=1,DIMQ2)
      do id=1, size(S_H)
        write(bu,err=99,iostat=ios) 
     $       (normis(id,ie),ie=1,DIMQ2)
      end do
      
      
      close (bu)
      return
 99   continue
      write(*,*) 'wrvec Error iostat', ios,ij,gam,if,ib,if1,binfname
      return
      end

C----------------------------------------------------------------------
      subroutine rdvecQ(ij,it,is,if,ib,if1,mysize,qvk,vr,vi,e,qmk
     $ ,F1s,normis) ! somethings seems not to work with rdvecQ
      implicit none
      include 'iam.fi'
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer mysize(DIMSIZ)
      integer ij,it,is,if,ib,i,ie,itop,bu,if1
      real*8 normis(DIMQ2)! stores the eigenvector**2 for each F1
      integer F1s(DIMQ2) ! stores the F1 quantum number contributions
      integer qmk(DIMQLP) 
      real*8  vr(DIMTOT),vi(DIMTOT),e
      integer done, found
      character*30 binfname     
      done =0 
      call getfu(bu)
c      write(0,*)'rdvec',bu
      call binnam(ij,is,if,ib,if1,binfname)
      open(bu,file=binfname,status='unknown',form='unformatted')
      read(bu) (mysize(ie),ie=1,DIMSIZ)
      read(bu) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,mysize(S_H))
      do i=1, mysize(S_H)
        if (done.eq.0) then
         read(bu) e, (qmk(ie),ie=1,DIMQLP),
     $        (vr(ie),vi(ie),ie=1,mysize(S_H))
         found=i
         if (qmk(Q_T).eq.it) then
         done=1
         end if
        else
         read(bu) ! stop reading if done-1, if vr and vi are found.
        end if
      end do
      if (done.eq.0)then
       write(0,*) 'eigenvector not found' !wraning removal Herbers2024
      end if
      read(bu) (F1s(ie),ie=1,DIMQ2)
      do i=1,  mysize(S_H)
        read(bu) 
     $       (normis(ie),ie=1,DIMQ2)
      if (i.eq.found) goto 20
      end do
 20   close (bu)
      return
      end
C ----------------------------------------------------------------------
      subroutine rdmatQ(ij,is,if,ib,if1,mysize,qvk,vr,vi,e,qmk
     $ ,F1s,normis)
      implicit none
      include 'iam.fi'
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer mysize(DIMSIZ)
      integer ij,is,if,ib,i,ie,itop,bu,if1,id
      real*8 normis(DIMTOT,DIMQ2)! stores the eigenvector**2 for each F1
      integer F1s(DIMQ2) ! stores the F1 quantum number contributions
      integer qmk(DIMTOT,DIMQLP) 
      real*8  vr(DIMTOT,DIMTOT),vi(DIMTOT,DIMTOT),e(DIMTOT)
      character*30 binfname     

      call getfu(bu)
c      write(0,*)'rdmat',bu
      call binnam(ij,is,if,ib,if1,binfname)
      open(bu,file=binfname,status='unknown',form='unformatted')
      read(bu) (mysize(ie),ie=1,DIMSIZ)
      read(bu) (qvk(ie,Q_K),
     $     (qvk(ie,Q_V+itop),itop=1,DIMTOP)
     $     ,ie=1,mysize(S_H))
      do i=1, mysize(S_H)
        read(bu) e(i), (qmk(i,ie),ie=1,DIMQLP),
     $       (vr(ie,i),vi(ie,i),ie=1,mysize(S_H))
      end do
      read(bu) (F1s(ie),ie=1,DIMQ2)
      do id=1, mysize(S_H)
        read(bu) (normis(id,ie),ie=1,DIMQ2)
        
      end do
      
      close (bu)
      return
      end
C----------------------------------------------------------------------
      subroutine wrtori(tori,ib)
      implicit none
      include 'iam.fi'
      integer ib
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer ie,itop,ntop,isig,sigma,ivr,ivc,ikr,ikc,bu
      character*30 binfname     
      integer myand
      external myand

      call getfu(bu)
      ntop=ctlint(C_NTOP)
      write(binfname,'(A,1I1)')
     $     'tori.b',ib
      if (myand(ctlint(C_PRI),AP_ST).ne.0)
     $     write(*,*)' Writing ',binfname
      open(bu,file=binfname,status='unknown',form='unformatted')
      write(bu) (size(ie),ie=1,DIMSIZ)
      write(bu)
     $     ((gamma(isig,itop),itop=1,ntop),isig=1,size(S_G))

      do itop=1, ntop
        do isig=1, size(S_G)
          sigma=gamma(isig,itop)
          do ivr=1, size(S_V+itop)
            do ikr=-size(S_MAXK), size(S_MAXK)
              do ivc=1, size(S_V+itop)
                do ikc=-size(S_MAXK), size(S_MAXK)
                  write(bu) tori(ikr,ikc,ivr,ivc,sigma,itop)
     $                 ,ikr,ikc,ivr,ivc
                end do
              end do
            end do
          end do
        end do
      end do
      close (bu)
      
      return
      end

C----------------------------------------------------------------------
      subroutine rdtori(tori,ib)
      implicit none
      include 'iam.fi'
      integer ib
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer ie,itop,ntop,isig,sigma,ivr,ivc,ikr,ikc,bu
      character*30 binfname     
      integer myand
      external myand

      call getfu(bu)
      ntop=ctlint(C_NTOP)
      write(binfname,'(A,I1)')
     $     'tori.b',ib
      if (myand(ctlint(C_PRI),AP_ST).ne.0)
     $     write(*,*)' Reading Tori ',binfname
      open(bu,file=binfname,status='unknown',form='unformatted')
      read(bu) (size(ie),ie=1,DIMSIZ)
      read(bu)
     $     ((gamma(isig,itop),itop=1,ntop),isig=1,size(S_G))
      do itop=1, ntop
        do isig=1, size(S_G)
          sigma=gamma(isig,itop)
          do ivr=1, size(S_V+itop)
            do ikr=-size(S_MAXK), size(S_MAXK)
              do ivc=1, size(S_V+itop)
                do ikc=-size(S_MAXK), size(S_MAXK)
                  read(bu) tori(ikr,ikc,ivr,ivc,sigma,itop)
                end do
              end do
            end do
          end do
        end do
      end do
      close (bu)
      
      return
      end

C     -------------------------------------------------------------
      subroutine intal2(mux,muy,muz,temp)
      implicit none
      include 'iam.fi'
      real*8  muz, mux, muy
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer sizeup(DIMSIZ),sizelo(DIMSIZ)
      integer qvkup(DIMTOT,Q_K:Q_V+DIMTOP),qvklo(DIMTOT,Q_K:Q_V+DIMTOP)
C     real*8  vrup(DIMTOT),vrlo(DIMTOT),viup(DIMTOT),vilo(DIMTOT)
      real*8  tmql,tmq,strr,strlg,str,fac,frq,dgn,cmc,tmc,elow !Herbers2023
      integer qup(DIMTOT,DIMQLP),qlo(DIMTOT,DIMQLP)
      real*8  elo(DIMTOT),eup(DIMTOT),ints,freq
      real*8  hrup(DIMTOT,DIMTOT),hrlo(DIMTOT,DIMTOT)
      real*8  hiup(DIMTOT,DIMTOT),hilo(DIMTOT,DIMTOT)
      integer jup,jlo
      integer j1,j2,km1,km2,kp1,kp2,k1,k2,t1,t2,i,t01,t02,f1,f2
      integer fup,flo,tup,tlo,sup,slo,bup,blo,vup,vlo,isf,io
      integer f1up,f1lo !h2024
      integer iup(DIMVV),ilo(DIMVV)
      integer qsp(10,3*DIMTOT),ino,nno,iino,sp(0:DIMGAM)
      real*8  dsp(4,3*DIMTOT),sfrq
      real*8  kghz,bolt,temp,ener,statw,popl,hv,spcat !Herbers2023
      real*8  delo,deup,de1,de2 !Herbers2023/2024 
      real*8  enull(size(S_G)+1)
      real*8  sj,sj2       
      integer f11, f12
      parameter (kghz=20.8366191)!Herbers2023 updated from 20.8364 
      parameter (cmc=29979.2458)!Herbers2023
      parameter (tmc=-1.43878)!Herbers2023
      integer endidx
      integer stepflo,stepf1lo
      integer stepfup,stepf1up
      integer runf,runf1
      integer skip

      real*8 dump(DIMTOT), dump2(DIMTOT)
      real*8 sjsum
      real*8 sj2sum
      integer F1slo(DIMQ2)
      real*8 storesj(DIMQ2,DIMQ2) ! for storage of sixj values in the two nucleus approximation case to avoid redundant calculations.
      real*8 nF1slo(DIMTOT,DIMQ2) ! here only DIMQ, because only vector is read not matrix 
      integer F1sup(DIMQ2)
      real*8 nF1sup(DIMTOT,DIMQ2)
      integer jq,iq! runnin indicies for sj sj2 treatment
      
        
C      character*6 fno

      enull=0.0
      fup=-1
      flo=-1
      isf=1
c      if (ctlint(C_NTOP).eq.0) isf=0
      do blo=1,size(S_NB) ! run over all B
        write(*,*)
        write(*,*) 'Linestrengths**2 str**2 are reduced'
        write(*,*) 'transition dipole moments squared'
        write(*,*)
        write(*,*) 'd_pop is log base 10 of calculated population'
        write(*,*) 'difference, based on qsum'
        write(*,*)
        write(*,*) ' spcat is intensity as expected from spcat'
        write(*,'(A,I2,28X,2A12,A12,4A12)') '-- B',blo,'Freq','Split  '
     $       ,'log10_str**2','log10_total','stat.w.','d_pop.','spcat'
        call rdtori(tori,blo) !calculate tori for B
        bup=blo               !bup = blo by default in  intensity calculations 
        do vlo=1,size(S_VV)   !Run over all V
          vup=vlo
              if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
               runf1=ctlint(C_SPIN)
               endidx=size(S_MAXK)-ctlint(C_SPIN)-ctlint(C_SPIN2)
               if ((ctlint(C_SPIN2).ne.0)) then
                runf=ctlint(C_SPIN2)
               else
                runf=0
               end if
              else
              endidx=size(S_MAXK)
              runf1=0
              runf=0
              end if
            do jlo=0, endidx
            do jup=jlo,min(jlo+1,endidx)
            do stepf1lo=0,runf1
            do stepf1up=0,runf1
            do stepflo=0,runf
            do stepfup=0,runf
             nF1slo=0.0
             nF1sup=0.0
             F1slo=-1
             F1sup=-1
             if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
              f1lo=2*jlo-ctlint(C_SPIN)+2*(stepf1lo)
              f1up=2*jup-ctlint(C_SPIN)+2*(stepf1up)
             endif
             if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN2).ne.0)) then
             
              flo=f1lo-ctlint(C_SPIN2)+2*(stepflo)
              fup=f1up-ctlint(C_SPIN2)+2*(stepfup)
             else
              flo=f1lo
              fup=f1up
             endif
            
            if ((ctlint(C_DW).ne.3).or.((abs(fup-flo).le.2).and.
     $      (fup.ge.0).and.(flo.ge.0).and.(f1up.ge.0).and.(f1lo.ge.0))
     $      .or.((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).eq.0))) then !Will end if at very end
            skip=0
            
            if ((ctlint(C_SPIN2).ne.0).and.(ctlint(C_DW).eq.3)) then !IF-CHECK-EXIST-SPIN2
            
            if (flo.lt.ctlint(C_SPIN2))then ! checking for states that do not exist
             if ((abs(flo+f1lo)).lt.ctlint(C_SPIN2))then
              skip=1
             end if
            end if
            
            if (fup.lt.ctlint(C_SPIN2))then ! checking for states that do not exist
             if ((abs(fup+f1up)).lt.ctlint(C_SPIN2))then
              skip=1
             end if
            end if
            
            end if!END-CHECK-EXIST-SPIN2
            
            if ((ctlint(C_SPIN).ne.0).and.(ctlint(C_DW).eq.3)) then !IF-CHECK-EXIST-SPIN
            
            if (f1lo.lt.ctlint(C_SPIN))then ! checking for states that do not exist
             if ((abs(f1lo+2*jlo)).lt.ctlint(C_SPIN))then
              skip=1
             end if
            end if
            
            if (f1up.lt.ctlint(C_SPIN))then ! checking for states that do not exist
             if ((abs(f1up+2*jup)).lt.ctlint(C_SPIN))then
              skip=1
             end if
            end if
            
            end if !END-CHECK-EXIST-SPIN
            
            if (skip.ne.1) then!IF-SKIP

             
              ino=0
              do slo=0, size(S_G)              ! run over all species starting with slo=0 and sup is equal to slo
                sup=slo
                sp(slo)=ino                    !what is this? ino starts as zero

            if (((ctlint(C_SPIN).ne.0).and.(ctlint(C_SPIN2).ne.0))
     $          .and.((flo.ne.-1).and.(fup.ne.-1))) then
            call rdmatQ(jup,sup,fup,bup,f1up,sizeup
     $                ,qvkup,hrup,hiup,eup,qup,F1sup,nF1sup)
            call rdmatQ(jlo,slo,flo,blo,f1lo
     $                 ,sizelo,qvklo,hrlo,hilo,elo,qlo,F1slo,nF1slo)
            else 
            call rdmat(jup,sup,fup,bup,f1up
     $                 ,sizeup,qvkup,hrup,hiup,eup,qup)  !rdmat(ij,is,if,ib,mysize,qvk,vr,vi,e,qmk) does
            call rdmat(jlo,slo,flo,blo,f1lo
     $                 ,sizelo,qvklo,hrlo,hilo,elo,qlo)  !read(bu) e(i), (qmk(i,ie),ie=1,DIMQLP),    $       (vr(ie,i),vi(ie,i),ie=1,mysize(S_H)) which should aquire energies elo eup for certain quantum numbers
            end if
            if ((jlo.eq.0).and.(jup.eq.0)) then ! Herbers2023 ... referebce energy level for pop calculation
             if (enull(sup+1).ne.0.0) then
C             DoNothing!
             else
              enull(sup+1) = elo(1)!Herbers2023 sup+1 is required since fortran arrays start at 1... so 1 is S 0...
             end if
!             if (sup.eq.1) then
!             enull(sup) = elo(1) ! for some reason rigid E=0 takes an akward value of 0. I will put reference level to S 1, though this is not exactly correct, it should be okay. I think this could be fixed looking at the J=0 case in rdmat for sup=0 
!             end if              ! 27.08.2024: this got fixed by reloading the matrices in the second loop.
            end if
                do i=1, DIMVV
                  ilo(i)=0 !set all index_lo to zero
                end do
                
                storesj=-100.0 ! for two nucleus quadrupole coupling, the sj values can be stored for calculations within each J block to avoid redundant calculations.
                do tlo=1, sizelo(S_H) ! S_H is likely the size of J Block in the Hamiltonian (Size_Hamiltonian)
                  ilo(qlo(tlo,Q_V1))=ilo(qlo(tlo,Q_V1))+1 !set index_lo corresponding to related tlo to +1
                  do i=1, DIMVV
                    iup(i)=0 ! set all index_up to zero
                  end do
                  

                   
                  do tup=1, sizeup(S_H)
                  

                   
                    iup(qup(tup,Q_V1))=iup(qup(tup,Q_V1))+1 ! set index_up dependent on tup
                    freq=eup(tup)-elo(tlo) !frequency is difference in energy
                    if ((dabs(freq).lt.ctlpar(C_FRQLO))
     $                   .or.(dabs(freq).gt.ctlpar(C_FRQUP))) goto 10 ! if freq not in range, then next steps are skipped.
                    if ((tup.le.tlo).and.(jup.eq.jlo)) goto 10 ! this is some quantum number mismatch condition to avoid double counting in q branch for negative and positive freqs.
                    if ((qup(tup,Q_V1).ne.vup).or.
     $                   (qlo(tlo,Q_V1).ne.vlo)) goto 10 ! some condition regarding the torisional quantum number.

                    if (jup.gt.jlo) then                     
                      call calir(ints,tori,jup,jlo,sup,muz, mux, muy !calculate intensities r band transition delta J=1
     $                     ,sizeup,qvkup,hrup(1,tup),hiup(1,tup)
     $                     ,sizelo,qvklo,hrlo(1,tlo),hilo(1,tlo))
                    end if                                   
                    if (jup.eq.jlo) then                     
                      call caliq(ints,tori,jup,jlo,sup,muz, mux, muy !calculate intensities q brand transition delta J=0
     $                     ,sizeup,qvkup,hrup(1,tup),hiup(1,tup)
     $                     ,sizelo,qvklo,hrlo(1,tlo),hilo(1,tlo))
                    end if
C                 !EXTRA sj calculations for six j symbols in case of hyperfines   (Sorry for missing indentation.)                 
C                 !EXTRA sj calculations for six j symbols in case of hyperfines                    
C                 !EXTRA sj calculations for six j symbols in case of hyperfines                    
      if ((ctlint(C_SPIN2).ne.0).and.(ctlint(C_SPIN).ne.0)
     $                                   .and.(flo.ne.-1)) then ! approximate to nuclei intensities, by weighting the wigner symbol prorducts by the contribution of the various F1s.
      sjsum=0.0
      sj2sum=0.0
      
C            call rdvecQ(jlo,tlo,slo,flo,blo,f1lo,sizelo
C     $             ,qvklo,dump,dump2,elo(tlo),qlo(tlo,:)
C     $             ,F1slo,nF1slo(tlo,:))
C            call rdvecQ(jup,tup,sup,fup,bup,f1up,sizeup
C     $             ,qvkup,dump,dump2,eup(tup),qup(tup,:)
C     $             ,F1sup,nF1sup(tup,:))
      
      do iq=1,DIMQ2 
        do jq=1, DIMQ2 
          if ((F1sup(iq).ne.(-1)).and.(F1slo(jq).ne.(-1))) then
           if (storesj(iq,jq).lt.(-90.0)) then
           call sixj(2*jup,F1sup(iq),ctlint(C_SPIN),F1slo(jq)
     $                                              ,2*jlo,2,sj)
           call sixj(F1sup(iq),fup,ctlint(C_SPIN2),flo,F1slo(jq)
     $                                                   ,2,sj2)
           sj=sj*((F1sup(iq)+1)*(F1slo(jq)+1))**0.5
           sj2=sj2*((fup+1)*(flo+1))**0.5
           storesj(iq,jq)=abs(sj*sj2)
           end if
           sjsum=sjsum+abs(storesj(iq,jq)
     $                        *(nF1slo(tlo,jq)*nF1sup(tup,iq)))
          end if
        end do
      end do
C      ints=ints*sjsum
       ints=ints*sjsum
      
      else
       if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
         if ((ctlint(C_SPIN2).eq.0).or.(flo.eq.-1)) then
          sj2=1.0
         else
          call sixj(f1up,fup,ctlint(C_SPIN2),flo,f1lo,2,sj2)
          sj2=sj2**2
          sj2=abs(sj2*(fup+1)*(flo+1))
         end if
         if (f1lo.eq.-1) then
          sj=1.0
         else
          call sixj(2*jup,f1up,ctlint(C_SPIN),f1lo,2*jlo,2,sj)
          sj=sj**2
          sj=abs(sj*(f1up+1)*(f1lo+1))!/(ctlint(C_SPIN)+1)) !Pickett does not use devision by 2I+1
         end if
         ints=ints*sj*sj2
       end if
      end if
C                 !FINISHED EXTRA sj calculations for six j symbols in case of hyperfines                    
C                 !FINISHED EXTRA sj calculations for six j symbols in case of hyperfines                    
C                 !FINISHED EXTRA sj calculations for six j symbols in case of hyperfines   

                    if ((flo.ge.0).and.(ctlint(C_SPIN).ne.0)) then
                    statw=min(flo,fup)+1 !f is already 2 F
                    else
                    statw=2*jlo+1
                    endif
                    delo=elo(tlo)-enull(sup+1)!Herbers2023
                    deup=eup(tup)-enull(sup+1)!Herbers2023

                    ener=min(elo(tlo),eup(tup))
                    if (freq.gt.0.0) then
                        hv  =1.0d0-exp(-freq/(kghz*temp))!Herber2023 moved
                    else
                        hv  =1.0d0-exp(freq/(kghz*temp))!Herber2023 moved
                    end if
                    popl=(statw)*exp(-delo/(kghz*temp))
                    popl=popl-(statw)*exp(-deup/(kghz*temp))!delta m=0 , degeneracy should be entered the same here for jup and jlo.
                    popl=abs(popl)/ctlpar(C_QSUM)                         
                    bolt=popl*ints!Herber2023 moved
                    if (ints.eq.0) then
                    bolt=-9999
                    else
                    bolt=popl*ints
                    bolt=log10(bolt)!Herbers2023                    
                    endif
                    if (bolt.ge.ctlpar(C_INTLM)) then
                      ino=ino+1
                      if (ino.gt.3*DIMTOT) stop 'INTALL:error:ino'
                      qsp(1,ino)=sup
                      qsp(2,ino)=tlo
                      qsp(3,ino)=qlo(tlo,Q_K)
                      qsp(4,ino)=ilo(qlo(tlo,Q_V1))/2
                      qsp(5,ino)=(2*jlo+2-ilo(qlo(tlo,Q_V1)))/2
                      qsp(6,ino)=tup
                      qsp(7,ino)=qup(tup,Q_K)
                      qsp(8,ino)=iup(qup(tup,Q_V1))/2
                      qsp(9,ino)=(2*jup+2-iup(qup(tup,Q_V1)))/2
                      qsp(10,ino)=0
                      dsp(1,ino)=freq
                      dsp(2,ino)=ints
                      dsp(3,ino)=min(elo(tlo),eup(tup))
                    end if
 10                 continue
                  end do        ! tup
                end do          ! tlo
                ino=ino+1
                qsp(1,ino)=sup
                qsp(2,ino)=-1
                qsp(6,ino)=-1
              end do            ! Sup
              nno=ino
              do iino=1, nno
                tlo=qsp(2,iino)
                tup=qsp(6,iino)
C                write(*,*)'tlo,tup',tlo,tup

                do sup=0, size(S_G) ! a new S lop: energies will have to be reloaded 
        call rdmat(jup,sup,fup,bup,f1up,sizeup,qvkup,hrup,hiup,eup,qup)  ! energies have to be reloaded for each S species for population calculations to be precise.
        call rdmat(jlo,sup,flo,blo,f1lo,sizelo,qvklo,hrlo,hilo,elo,qlo)  
        
        
c                  do itry=1,2
                  do io=1, nno
                    ino=sp(sup)+io
C     check for end of sup
                    if ((qsp(2,ino).eq.-1)
     $                   .and.(qsp(6,ino).eq.-1)) goto 30
                    if ((tlo.eq.-1).and.(tup.eq.-1)
     $                   .and.(qsp(10,ino).eq.0)) goto 50
C     check for matching tup and tlo
                    if ((tlo.ne.qsp(2,ino))
     $                   .or.(tup.ne.qsp(6,ino))) goto 40
 50                 continue

                    freq=dsp(1,ino)
                    ints=dsp(2,ino)
                    ener=dsp(3,ino)
                    qsp(10,ino)=1
                    if (freq.gt.0.0) then
                      j1 =jup
                      t1 =qsp(6,ino)
                      k1 =qsp(7,ino)
                      km1=qsp(8,ino)
                      kp1=qsp(9,ino)
                      f1 = fup
                      f11 = f1up
                      
                      j2 =jlo
                      t2 =qsp(2,ino)
                      k2 =qsp(3,ino)
                      km2=qsp(4,ino)
                      kp2=qsp(5,ino)
                      f2 = flo
                      f12 = f1lo
                      
                      de2=elo(t2)-enull(sup+1)!Herbers2024
                      de1=eup(t1)-enull(sup+1)!Herbers2024
                      elow=de2*1000.
                      
                      if (sup.eq.isf) then
                        sfrq=freq
                        t01=t1
                        t02=t2
                      end if
                    else
                      j2 =jup
                      t2 =qsp(6,ino)
                      k2 =qsp(7,ino)
                      km2=qsp(8,ino)
                      kp2=qsp(9,ino)
                      f2 =fup
                      f12 = f1up
                      
                      j1 =jlo
                      t1 =qsp(2,ino)
                      k1 =qsp(3,ino)
                      km1=qsp(4,ino)
                      kp1=qsp(5,ino)
                      f1 =flo
                      f11 = f1lo
                      
                      freq=-freq
                      
                      de1=elo(t1)-enull(sup+1)!Herbers2024
                      de2=eup(t2)-enull(sup+1)!Herbers2024
                      elow=de2*1000.
                      
                      if (sup.eq.isf) then
                        sfrq=freq
                        t01=t1
                        t02=t2
                      end if
                    end if
                    
                    

                    
                    strr=ints                ! following the notation in line 669 of calcat from picketts calpgm
                    fac = (4.16231e-5)/ctlpar(C_QSUM)  ! see https://spec.jpl.nasa.gov//ftp//pub/calpgm/calcat.c
                    frq = freq*1000.          ! all of this is Herbers2023
                    tmq = tmc/(temp*cmc)
                    tmql = tmq*0.43429448
                    dgn = 1.0! (2*jlo+1) ! 1.0 seems to work
                    
                    str=dgn*strr*fac*frq*(1. - exp(tmq * frq))
                    popl=(statw)*exp(-de1/(kghz*temp)) !always jlo, lower of two js for delta m equals zero
                    popl=popl-(statw)*exp(-de2/(kghz*temp))!delta m=0 , degeneracy should be entered the same here.
                    popl=abs(popl)/ctlpar(C_QSUM)                     
                    if (str.eq.0) then
                    strlg=-9999
                    bolt=-9999
                    else
                    strlg=log10(str+0.)+tmql*elow
                    bolt=popl*ints
                    bolt=log10(bolt)!Herbers2023                    
                    endif
                    spcat=strlg                    


C                    if ((jlo.eq.2).and.(jup.eq.3).and.
C     $                   (flo.eq.1).and.(fup.eq.3)) then
C                     write(0,*) tup,tlo,sup,elow,tmql,sj,str
C                    end if


                    if (ctlint(C_SPIN).ne.0) then !START-SPIN-1
                    if (ctlint(C_SPIN2).ne.0) then !START-SPIN-2
                    if ((t02.eq.t2).and.(t01.eq.t1).and.(sup.gt.isf))!START-t !Both SPINS are present.
     $                   then
             write(*,'(39X,(A,I2),10X,F13.6,F10.4,5F12.4,2(A,2I3))')
     $               '  S',sup
     $               ,freq,(freq-sfrq)*1000.0,log10(ints),bolt,statw !Herbers2023 this is int 3
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    else !ELSE-t
                      if (sup.gt.0) then !START-S
             write(*,'(2(3I3,1X),A4,2I3,A3,2I3,3(A,I2),
     $                F13.6,10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,' F1 ', f11, f12
     $               ,' F ', f1, f2
     $               ,'  S',sup,'  V',vup,'  B',bup
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                      end if !END-S
                    end if!END-t
                    if (sup.eq.0) then!START-S
             write(*,'(2(3I3,1X),A4,2I3,A3,2I3,A15,F13.6,
     $                  10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,' F1 ', f11, f12,' F '
     $               , f1, f2,
     $               '  rigid  '
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    end if!END-S
                    
                    else ! ELSE-SPIN2
                    
                    if ((t02.eq.t2).and.(t01.eq.t1).and.(sup.gt.isf))!START-t
     $                   then
             write(*,'(30X,(A,I2),10X,F13.6,F10.4,5F12.4,2(A,2I3))')
     $               '  S',sup
     $               ,freq,(freq-sfrq)*1000.0,log10(ints),bolt,statw !Herbers2023 this is int 3
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    else !ELSE-t
                      if (sup.gt.0) then !START-S
             write(*,'(2(3I3,1X),A3,2I3,3(A,I2),
     $                F13.6,10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,' F ', f11, f12,
     $               '  S',sup,'  V',vup,'  B',bup
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                      end if !END-S
                    end if!END-t
                    if (sup.eq.0) then!START-S
             write(*,'(2(3I3,1X),A3,2I3,A15,F13.6,
     $                  10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,' F ', f11, f12, 
     $               '  rigid  '
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    end if!END-S
                    
                    end if !END-SPIN2
                    
                    else !ELSE-SPIN1
                    
                     
                    if ((t02.eq.t2).and.(t01.eq.t1).and.(sup.gt.isf)) !No hyperfine
     $                   then
             write(*,'(20X,(A,I2),10X,F13.6,F10.4,5F12.4,2(A,2I3))')
     $               '  S',sup
     $               ,freq,(freq-sfrq)*1000.0,log10(ints),bolt,statw !Herbers2023 this is int 3
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    else
                      if (sup.gt.0) then
             write(*,'(2(3I3,1X),3(A,I2),F13.6,10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,
     $               '  S',sup,'  V',vup,'  B',bup
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                      end if
                    end if
                    if (sup.eq.0) then
             write(*,'(2(3I3,1X),A15    ,F13.6,10X,5F12.4,2(A,2I3))')
     $                 j1 ,km1, kp1,j2 ,km2, kp2,
     $               '  rigid  '
     $               ,freq                   ,log10(ints),bolt,statw !Herbers2023
     $               ,log10(popl),spcat,'  K',k1,k2,'  t',t1,t2
                    end if
                    
                    
                    endif
 40                 continue
                  end do        ! io
c                  end do
 30               continue
                end do          ! sup
                if ((tlo.eq.-1).and.(tup.eq.-1)) goto 20
              end do            ! ino
 20           continue
            end if
            end if              ! end at very end.
            end do
            end do
            end do              ! Fup
            end do              ! Flo
            end do              ! Jup
          end do                ! Jlo
          write(*,*)
        end do                  ! Vlo
      end do                    ! Blo
      write(*,*) ' total is product of pop dif (d_pop) and reduced'
      write(*,*) 'transition dipole moment squared (str**2)'
      write(*,*)
      return
      end
C                fno=char(vup+64)//char(j1+64)//char(t1+64)
C     $               //char(j2+64)//char(t2+64)//char(sup+64)
C ----------------------------------------------------------------------
C     -------------------------------------------------------------
      subroutine intall(mux,muy,muz,temp)
      implicit none
      include 'iam.fi'
      real*8  muz, mux, muy
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer sizeup(DIMSIZ),sizelo(DIMSIZ)
      integer qvkup(DIMTOT,Q_K:Q_V+DIMTOP),qvklo(DIMTOT,Q_K:Q_V+DIMTOP)
C     real*8  vrup(DIMTOT),vrlo(DIMTOT),viup(DIMTOT),vilo(DIMTOT)
      integer qup(DIMTOT,DIMQLP),qlo(DIMTOT,DIMQLP)
      real*8  elo(DIMTOT),eup(DIMTOT),ints,freq
      real*8  hrup(DIMTOT,DIMTOT),hrlo(DIMTOT,DIMTOT)
      real*8  hiup(DIMTOT,DIMTOT),hilo(DIMTOT,DIMTOT)
      integer jup,jlo
      integer checkjjff(0:DIMJ,0:DIMJ,0:DIMJ+(DIMQ+1)/2
     $  ,0:DIMJ+(DIMQ+1)/2)
      integer j1,j2,km1,km2,kp1,kp2,k1,k2,t1,t2,i,f1,f2
      integer f11,f12 !herbers2024
      integer fup,flo,tup,tlo,sup,slo,bup,blo,vup,vlo,isf
      integer f1up,f1lo !herbers2024
      integer skip
      integer iup(DIMVV),ilo(DIMVV)
      real*8  kghz,bolt,temp,popl,statw,hv,spcat
      real*8  tmql,tmq,strr,strlg,str,fac,frq,dgn,cmc,tmc,elow !Herbers2023
      real*8  delo,deup,sj,sj2
      real*8  enull(size(S_G)+1) !Herbers 2023
      parameter (kghz=20.8366191)!Herbers2023
      parameter (cmc=29979.2458)!Herbers2023
      parameter (tmc=-1.43878)!Herbers2023
      integer endidx
      integer stepflo,stepf1lo
      integer stepfup,stepf1up
      integer runf,runf1
      integer hfup, hflo
      integer hf1up, hf1lo
      real*8 dump(DIMTOT), dump2(DIMTOT)
      real*8 sjsum
      real*8 sj2sum
      integer F1slo(DIMQ2)
      real*8 storesj(DIMQ2,DIMQ2)
      real*8 nF1slo(DIMTOT,DIMQ2) ! here only DIMQ, because only vector is read not matrix 
      integer F1sup(DIMQ2)
      real*8 nF1sup(DIMTOT,DIMQ2)
      integer jq,iq! runnin indicies for sj sj2 treatment
      
      
      fup=-1
      flo=-1
      f1lo=-1
      f1up=-1
      isf=0
      enull=0.0
      if (ctlint(C_NTOP).gt.0) isf=1
      
      do blo=1,size(S_NB)
          write(*,*) 'Linestrengths**2 str**2 are reduced'
        write(*,*) 'transition dipole moments squared'
        write(*,*)
        write(*,*) 'd_pop is log base 10 of calculated population'
        write(*,*) 'difference, based on qsum'
        write(*,*)        
        write(*,*) ' spcat is intensity as expected from spcat'
        if ((ctlint(C_SPIN).ne.0).and.(ctlint(C_SPIN2).ne.0)) then  ! Depending on the number of spins the space before Freq is printed is larger
        write(*,'(A,I2,44X,A12,5A12)') '-- B',blo,'Freq ','log10_str**2'
     $       ,'log10_total','stat.w.','d_pop','spcat'
        else 
         if ((ctlint(C_SPIN).ne.0).and.(ctlint(C_DW).eq.3)) then
        write(*,'(A,I2,34X,A12,5A12)') '-- B',blo,'Freq ','log10_str**2'
     $       ,'log10_total','stat.w.','d_pop','spcat'
         else
        write(*,'(A,I2,25X,A12,5A12)') '-- B',blo,'Freq ','log10_str**2'
     $       ,'log10_total','stat.w.','d_pop','spcat'
         end if
        end if
        bup=blo
        call rdtori(tori,blo)
        do vlo=1,size(S_VV)
          vup=vlo
          do slo=isf,size(S_G)
            sup=slo
            checkjjff=-1
            
              if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
               runf1=ctlint(C_SPIN)
               endidx=size(S_MAXK)-ctlint(C_SPIN)-ctlint(C_SPIN2)
               if ((ctlint(C_SPIN2).ne.0)) then
                runf=ctlint(C_SPIN2)
               else
                runf=0
               end if
              else
              endidx=size(S_MAXK)
              runf1=0
              runf=0
              end if
              
              
            do jlo=0, endidx
            do jup=jlo,min(jlo+1,endidx)
            do stepf1lo=0,runf1
            do stepf1up=0,runf1
            do stepflo=0,runf
            do stepfup=0,runf
             if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
              f1lo=2*jlo-ctlint(C_SPIN)+2*(stepf1lo)
              f1up=2*jup-ctlint(C_SPIN)+2*(stepf1up)
             endif
             if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN2).ne.0)) then
             
              flo=f1lo-ctlint(C_SPIN2)+2*(stepflo)
              fup=f1up-ctlint(C_SPIN2)+2*(stepfup)
             else
              flo=f1lo
              fup=f1up
             endif
            
            if ((ctlint(C_DW).ne.3).or.((abs(fup-flo).le.2).and.
     $      (fup.ge.0).and.(flo.ge.0).and.(f1up.ge.0).and.(f1lo.ge.0))
     $      .or.((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).eq.0))) then !Will end if at very end
            skip=0
            
            if ((ctlint(C_SPIN2).ne.0).and.(ctlint(C_DW).eq.3)) then !IF-CHECK-EXIST-SPIN2
            
            if (flo.lt.ctlint(C_SPIN2))then ! checking for states that do not exist
             if ((abs(flo+f1lo)).lt.ctlint(C_SPIN2))then
              skip=1
             end if
            end if
            
            if (fup.lt.ctlint(C_SPIN2))then ! checking for states that do not exist
             if ((abs(fup+f1up)).lt.ctlint(C_SPIN2))then
              skip=1
             end if
            end if
            
            end if!END-CHECK-EXIST-SPIN2
            
            if ((ctlint(C_SPIN).ne.0).and.(ctlint(C_DW).eq.3)) then !IF-CHECK-EXIST-SPIN
            
            if (f1lo.lt.ctlint(C_SPIN))then ! checking for states that do not exist
             if ((abs(f1lo+2*jlo)).lt.ctlint(C_SPIN))then
              skip=1
             end if
            end if
            
            if (f1up.lt.ctlint(C_SPIN))then ! checking for states that do not exist
             if ((abs(f1up+2*jup)).lt.ctlint(C_SPIN))then
              skip=1
             end if
            end if
            
            end if !END-CHECK-EXIST-SPIN
            
            if (skip.ne.1) then!IF-SKIP

            if (((ctlint(C_SPIN).ne.0).and.(ctlint(C_SPIN2).ne.0))
     $          .and.((flo.ne.-1).and.(fup.ne.-1))) then
            call rdmatQ(jup,sup,fup,bup,f1up,sizeup
     $                ,qvkup,hrup,hiup,eup,qup,F1sup,nF1sup)
            call rdmatQ(jlo,slo,flo,blo,f1lo
     $                 ,sizelo,qvklo,hrlo,hilo,elo,qlo,F1slo,nF1slo)
            else
            call rdmat(jup,sup,fup,bup,f1up,sizeup
     $                ,qvkup,hrup,hiup,eup,qup)
            call rdmat(jlo,slo,flo,blo,f1lo,sizelo
     $                ,qvklo,hrlo,hilo,elo,qlo)
            end if
             if ((jlo.eq.0).and.(jup.eq.0)) then ! Herbers2023 ... reference energy level for pop calculation
             if (enull(sup+1).ne.0.0) then !this condition was added, to have only one hyperfine J=0 as reference.
             ! Do nothing
             else
             enull(sup+1) = elo(1)!Herbers2023 sup+1 is required since fortran arrays start at 1... so 1 is S 0...
             end if
            end if
                do i=1, DIMVV
                  ilo(i)=0
                end do
                storesj=-100.0
                do tlo=1, sizelo(S_H)
                  ilo(qlo(tlo,Q_V1))=ilo(qlo(tlo,Q_V1))+1
                  do i=1, DIMVV
                    iup(i)=0
                  end do
                  do tup=1, sizeup(S_H)
                    iup(qup(tup,Q_V1))=iup(qup(tup,Q_V1))+1
                    freq=eup(tup)-elo(tlo)
                    if ((dabs(freq).lt.ctlpar(C_FRQLO))
     $                   .or.(dabs(freq).gt.ctlpar(C_FRQUP))) goto 10
                    if ((tup.le.tlo).and.(jup.eq.jlo)) goto 10
                    if ((qup(tup,Q_V1).ne.vup).or.
     $                   (qlo(tlo,Q_V1).ne.vlo)) goto 10
                    if (jup.gt.jlo) then                     
                      call calir(ints,tori,jup,jlo,sup,muz, mux, muy
     $                     ,sizeup,qvkup,hrup(1,tup),hiup(1,tup)
     $                     ,sizelo,qvklo,hrlo(1,tlo),hilo(1,tlo))
                    end if                                   
                    if (jup.eq.jlo) then                     
                      call caliq(ints,tori,jup,jlo,sup,muz, mux, muy
     $                     ,sizeup,qvkup,hrup(1,tup),hiup(1,tup)
     $                     ,sizelo,qvklo,hrlo(1,tlo),hilo(1,tlo))
                    end if
                    if ((flo.ge.0).and.(ctlint(C_SPIN).ne.0)) then
                    statw=min(flo,fup)+1
                    else
                    statw=2*jlo+1
                    endif
                    delo=elo(tlo)-enull(sup+1)!Herbers2023 Offset is corrected when internal rotation is involved. J=0 is set to E=0.
                    deup=eup(tup)-enull(sup+1)!Herbers2023 Offset is corrected when internal rotation is involved. J=0 is set to E=0.
                    
                    
C                    !START. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !START. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !START. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !START. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !START. Inserted here is the sixj treatment for hyperfine approximate intensities.
      if ((ctlint(C_SPIN2).ne.0).and.(ctlint(C_SPIN).ne.0)
     $                                   .and.(flo.ne.-1)) then ! approximate to nuclei intensities, by weighting the wigner symbol prorducts by the contribution of the various F1s.
      sjsum=0.0
      sj2sum=0.0
      
C            call rdvecQ(jlo,tlo,slo,flo,blo,f1lo,sizelo ! not needed anymore since already read by rdmatQ
C     $             ,qvklo,dump,dump2,elo(tlo),qlo(tlo,:)
C     $             ,F1slo,nF1slo(tlo,:))
C            call rdvecQ(jup,tup,sup,fup,bup,f1up,sizeup
C     $             ,qvkup,dump,dump2,eup(tup),qup(tup,:)
C     $             ,F1sup,nF1sup(tup,:))
      
      do iq=1,DIMQ2 
        do jq=1, DIMQ2 
          if ((F1sup(iq).ne.(-1)).and.(F1slo(jq).ne.(-1))) then
           if (storesj(iq,jq).lt.(-90.0)) then
           call sixj(2*jup,F1sup(iq),ctlint(C_SPIN),F1slo(jq)
     $                                              ,2*jlo,2,sj)
           call sixj(F1sup(iq),fup,ctlint(C_SPIN2),flo,F1slo(jq)
     $                                                   ,2,sj2)
           sj=sj*((F1sup(iq)+1)*(F1slo(jq)+1))**0.5
           sj2=sj2*((fup+1)*(flo+1))**0.5
           storesj(iq,jq)=abs(sj*sj2)
           end if
           sjsum=sjsum+abs(storesj(iq,jq)
     $                    *(nF1slo(tlo,jq)*nF1sup(tup,iq)))
          end if
        end do
      end do
C     ints=ints*sjsum
      ints=ints*sjsum**2
      
      else
       if ((ctlint(C_DW).eq.3).and.(ctlint(C_SPIN).ne.0)) then
         if ((ctlint(C_SPIN2).eq.0).or.(flo.eq.-1)) then
          sj2=1.0
         else
          call sixj(f1up,fup,ctlint(C_SPIN2),flo,f1lo,2,sj2)
          sj2=sj2**2
          sj2=abs(sj2*(fup+1)*(flo+1))
         end if
         if (f1lo.eq.-1) then
          sj=1.0
         else
          call sixj(2*jup,f1up,ctlint(C_SPIN),f1lo,2*jlo,2,sj)
          sj=sj**2
          sj=abs(sj*(f1up+1)*(f1lo+1))!/(ctlint(C_SPIN)+1)) !Pickett does not use devision by 2I+1
         end if
         ints=ints*sj*sj2
       end if
      end if
      
C                    !END. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !END. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !END. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !END. Inserted here is the sixj treatment for hyperfine approximate intensities.
C                    !END. Inserted here is the sixj treatment for hyperfine approximate intensities.
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    if (freq.gt.0.0) then
                        frq = freq*1000.
                        !popl=exp(-elo(tlo)/(kghz*temp))!Herber2023 moved
                        
                        popl=statw*exp(-delo/(kghz*temp))
                        popl=popl-statw*exp(-deup/(kghz*temp))!Do no use jup here in degeneracy. Since we do deltaM=0 the degeneracies of upper and lowe should be put the same.
                        popl=abs(popl)/ctlpar(C_QSUM) 
                        
                    else
                        !popl=exp(-eup(tup)/(kghz*temp))
                        frq = -freq*1000.
                        popl=statw*exp(-delo/(kghz*temp))
                        popl=popl-statw*exp(-deup/(kghz*temp))
                        popl=abs(popl)/ctlpar(C_QSUM)     
                    end if
                    
                    

                    !hv  =1.0d0-exp(-abs(freq)/(kghz*temp))!Herber2023 moved
                    strr=ints                ! following the notation in line 669 of calcat from picketts calpgm
                    fac = (4.16231e-5)/ctlpar(C_QSUM)  ! see https://spec.jpl.nasa.gov//ftp//pub/calpgm/calcat.c
                    frq=frq          ! all of this is Herbers2023
                    tmq = tmc/(temp*cmc)
                    tmql = tmq*0.43429448
                    dgn = 1 !(2*jlo+1) ! 1.0 seems to work
                    elow=min(delo,deup)*1000
                    str=dgn*strr*fac*frq*(1. - exp(tmq * frq))
                    
                                        
                    
                    !bolt=statw*popl*hv*ints!Herber2023 moved

                    if (str.eq.0) then
                    strlg=-9999
                    spcat=strlg
                    bolt=-9999
                    else
                    strlg=log10(str+0.)+tmql*elow
                    spcat=strlg
                    bolt=popl*ints
                    bolt=log10(bolt)!Herbers2023                    
                    endif
                    
C                    if ((jlo.eq.2).and.(jup.eq.3).and.
C     $                   (flo.eq.1).and.(fup.eq.3)) then
C                    end if
                    
                    
                    
                    if (bolt.ge.ctlpar(C_INTLM)) then ! START INTENSITY-CONDITION
                      if (freq.gt.0.0) then
                        j1 =jup
                        f1 = fup
                        f11 = f1up
                        km1=iup(qup(tup,Q_V1))/2
                        kp1=(2*jup+2-iup(qup(tup,Q_V1)))/2
                        k1 =qup(tup,Q_K)
                        t1 =tup
                        j2 =jlo
                        f2 = flo
                        f12 = f1lo
                        km2=ilo(qlo(tlo,Q_V1))/2
                        kp2=(2*jlo+2-ilo(qlo(tlo,Q_V1)))/2
                        k2 =qlo(tlo,Q_K)
                        t2 =tlo
                      else
                        freq=-freq !Herbers2023
                        j2 =jup
                        f2 =fup
                        f12= f1up
                        km2=iup(qup(tup,Q_V1))/2
                        kp2=(2*jup+2-iup(qup(tup,Q_V1)))/2
                        k2 =qup(tup,Q_K)
                        t2 =tup
                        j1 =jlo
                        f1 =flo
                        f11 =f1lo
                        km1=ilo(qlo(tlo,Q_V1))/2
                        kp1=(2*jlo+2-ilo(qlo(tlo,Q_V1)))/2
                        k1 =qlo(tlo,Q_K)
                        t1 =tlo
                      end if
         if (ctlint(C_SPIN).ne.0) then !Start SPIN1
          if (ctlint(C_SPIN2).ne.0) then !START SPIN2
           if (.false.) then !START-GAM (gamma(sup,0).eq.0) then   ! deactivated the -K prints because I cant figure out what to use them for.
         write(*,'(2(I3,A3,I3,1X),A4,I3,I3,A3,I3,I3,2(A,I2)
     $   ,F13.6,5F12.4,A,I2,2(A,2I3))')
     $                  j1 ,' K ', k1, j2 ,' K ',k2,' F1 ', f11 , f12 
     $                       ,' F ', f1 , f2 
     $                     ,'  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl)!Herbers2023
     $                     ,spcat,'  B',bup,'  K',k1,k2,'  t',t1,t2
           else! ELSE-GAM
          write(*,'(2(3I3,1X),A4,I3,I3,A3,I3,I3,2(A,I2)
     $    ,F13.6,5F12.4,A,I2,2(A,2I3))')
     $                       j1 ,km1, kp1, j2 ,km2, kp2,' F1 ', f11 
     $                      , f12, ' F ', f1 , f2      
     $                      ,'  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl),spcat!Herbers2023
     $                     ,'  B',bup,'  K',k1,k2,'  t',t1,t2
           end if !END-GAM
          else           !ELSE SPIN 2
           if (.false.) then !START_GAM(gamma(sup,0).eq.0) then   ! deactivated the -K prints because I cant figure out what to use them for.
         write(*,'(2(I3,A3,I3,1X),A3,I3,I3,2(A,I2)
     $   ,F13.6,5F12.4,A,I2,2(A,2I3))')
     $                  j1 ,' K ', k1, j2 ,' K ',k2,' F ', f11 , f12 
     $                      ,'  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl),spcat!Herbers2023
     $                     ,'  B',bup,'  K',k1,k2,'  t',t1,t2
           else!ELSE-GAM
          write(*,'(2(3I3,1X),A3,I3,I3,2(A,I2)
     $    ,F13.6,5F12.4,A,I2,2(A,2I3))')
     $                       j1 ,km1, kp1, j2 ,km2, kp2,' F ', f11  
     $                      , f12,'  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl),spcat!Herbers2023
     $                     ,'  B',bup,'  K',k1,k2,'  t',t1,t2
          end if !END-IF-GAM
         end if !END-SPIN2
         else !ELSE-SPIN1
          if (.false.) then !START-GAM(gamma(sup,0).eq.0) then ! deactivated the -K prints because I cant figure out what to use them for.
          write(*,'(2(I3,A3,I3,1X),2(A,I2),F13.6,5F12.4,A,I2,2(A,2I3))')
     $                       j1 ,' K ', k1,j2 ,' K ',k2,
     $                       '  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl),spcat!Herbers2023
     $                     ,'  B',bup,'  K',k1,k2,'  t',t1,t2
          else!ELSE-GAM
          write(*,'(2(3I3,1X),2(A,I2),F13.6,5F12.4,A,I2,2(A,2I3))')
     $                       j1 ,km1, kp1,j2 ,km2, kp2,
     $                       '  S',sup,'  V',vup
     $                    ,freq,log10(ints),bolt,statw,log10(popl),spcat!Herbers2023
     $                     ,'  B',bup,'  K',k1,k2,'  t',t1,t2
          end if!END-GAM
          end if!END-SPIN1
                    end if!END-INTENSITY CONDITION
 10                 continue
                  end do
                end do
C             if (flo.ne.-1)checkjjff(jlo,jup,hflo,hfup)=0
C             end if
              
              end if ! END-IF-SKIP
              end if ! The new F conditions are ended here.
              
              end do
            end do
            end do
            end do
            end do
            end do
            write(*,*)
          end do
          write(*,*)
        end do
      end do
      write(*,*) 
      return
      end
C ----------------------------------------------------------------------
C     -------------------------------------------------------------
C                                                                      !Wigner 3j and 6j functions are taken from https://github.com/ogorton/wigner/blob/main/wigner.f90 and rewritten into subroutines.
      subroutine threej(two_j1, two_j2, two_j3,                        !Function from 2022 Oliver C. Gorton (ogorton@sdsu.edu) translated to subroutine
     $                  two_m1, two_m2, two_m3, tj) 

            ! Computes the Wigner 3-J symbol with arguments
            !   two_j1/2 two_j1/2 two_j3/2
            !   two_m1/2 two_m2/2 two_m3/2
            ! using clebsh-gordon vector-coupling coefficients.

            implicit none
            integer :: two_j1,two_j2,two_j3,two_m1,two_m2,two_m3
            real(kind=8) :: tj
            call vector_couple(two_j1, two_m1, two_j2
     $                   , two_m2, two_j3, -two_m3, tj)
            tj = (-1) **((two_j1 - two_j2 - two_m3)/2)
     $           /sqrt(dble(two_j3)+1d0)  * tj

      return
      end  

      subroutine vector_couple(two_j1, two_m1, two_j2,  
     $                         two_m2, two_jc, two_mc, cg)  

            ! Computes the Clebsh-Gordon vector-coupling coefficient
            !  (j1 m1 j2 m2 | j1 j2 jc mc)
            ! using algebraic expressions in factorials from Edmonds.

            implicit none
            integer :: two_j1,two_j2,two_jc,two_m1,two_m2,two_mc
            real(kind=8) :: cg, fac_prod, fac_sum, den
            real(kind=8) lf, ls
            integer :: t1, t2, t3, t4
            integer :: d1, d2, d3, d4, d5, d6
            integer :: dd1, dd2
            integer :: zmin, zmax, z

            cg = 0d0
            if (two_m1+two_m2 /= two_mc) return
            if (two_j1<0) return
            if (two_j2<0) return
            if (two_jc<0) return
            if (abs(two_m1)>two_j1) return
            if (abs(two_m2)>two_j2) return
            if (abs(two_jc)>two_jc) return
            if (mod(two_j1+two_m1,2) /= 0) return
            if (mod(two_j2+two_m2,2) /= 0) return
            if (mod(two_jc+two_mc,2) /= 0) return

            t1 = ( two_j1 + two_j2 - two_jc)/2
            t2 = ( two_j1 - two_j2 + two_jc)/2
            t3 = (-two_j1 + two_j2 + two_jc)/2
            t4 = ( two_j1 + two_j2 + two_jc)/2
            if (t1<0) return
            if (t2<0) return
            if (t3<0) return

            d1 = (two_j1 + two_m1)/2
            d2 = (two_j1 - two_m1)/2
            d3 = (two_j2 + two_m2)/2
            d4 = (two_j2 - two_m2)/2
            d5 = (two_jc + two_mc)/2
            d6 = (two_jc - two_mc)/2

            dd1 = (two_jc - two_j2 + two_m1)/2
            dd2 = (two_jc - two_j1 - two_m2)/2
            
            ls=0.0d0
            call logfac(d1,lf)
            ls=ls+lf
            call logfac(d2,lf)
            ls=ls+lf
            call logfac(d3,lf)
            ls=ls+lf
            call logfac(d4,lf)
            ls=ls+lf
            call logfac(d5,lf)
            ls=ls+lf
            call logfac(d6,lf)
            ls=ls+lf
            
            fac_prod=0.0
            call triangle(two_j1, two_j2, two_jc, fac_prod) 
            fac_prod = sqrt(dble(two_jc)+1d0) * fac_prod
            fac_prod = exp(0.5d0*ls)*fac_prod
            
            zmin = max(0, -dd1, -dd2)
            zmax = min(t1, d2, d3)

            fac_sum = 0d0
            do z = zmin, zmax
                ls=0.0d0
                 call logfac(z,lf)
                 ls=ls+lf
                 call logfac(t1-z,lf)
                 ls=ls+lf
                 call logfac(d2-z,lf)
                 ls=ls+lf
                 call logfac(d3-z,lf)
                 ls=ls+lf
                 call logfac(dd1+z,lf)
                 ls=ls+lf
                 call logfac(dd2+z,lf)
                 ls=ls+lf
                den = -ls
                fac_sum = fac_sum + (-1) ** (z) * exp(den)
            end do

            cg = fac_prod * fac_sum

      return
      end 







      subroutine sixj(two_j1,two_j2,two_j3,two_l1,two_l2,two_l3,sj) !Function from 2022 Oliver C. Gorton (ogorton@sdsu.edu) translated to subroutine
 
      ! Computes the wigner six-j symbol with arguments
      !    two_j1/2 two_j2/2 two_j3/2
      !    two_l1/2 two_l2/2 two_l3/2
      ! using explicit algebraic expressions from Edmonds (1955/7).
      
      implicit none
      integer :: two_j1,two_j2,two_j3,two_l1,two_l2,two_l3
      integer :: z, zmin, zmax
      integer  :: n2, n3, n4
      integer  :: d1, d2, d3, d4
      real(kind=8) :: sj
      real(kind=8) :: fac_sum, triangle_prod, num, den
      real(kind=8) :: delta, lf
      
      sj = 0d0
      if ( two_j1 < 0) return
      if ( two_j2 < 0) return
      if ( two_j3 < 0) return
      if ( two_l1 < 0) return
      if ( two_l2 < 0) return
      if ( two_l3 < 0) return
      if ( two_j1 < abs(two_j2-two_j3)) return
      if ( two_j1 > two_j2+two_j3 ) return
      if ( two_j1 < abs(two_l2-two_l3)) return
      if ( two_j1 > two_l2+two_l3 ) return
      if ( two_l1 < abs(two_j2-two_l3)) return
      if ( two_l1 > two_j2+two_l3 ) return
      if ( two_l1 < abs(two_l2-two_j3)) return
      if ( two_l1 > two_l2+two_j3 ) return
      
      n2 = (two_j1+two_j2+two_l1+two_l2)/2
      n3 = (two_j2+two_j3+two_l2+two_l3)/2
      n4 = (two_j3+two_j1+two_l3+two_l1)/2
      
      d1 = (two_j1+two_j2+two_j3)/2
      d2 = (two_j1+two_l2+two_l3)/2
      d3 = (two_l1+two_j2+two_l3)/2
      d4 = (two_l1+two_l2+two_j3)/2
      
      triangle_prod = 0.0
      call triangle(two_j1,two_j2,two_j3,delta)  
      triangle_prod = delta
      call triangle(two_j1,two_l2,two_l3,delta)
      triangle_prod = triangle_prod*delta
      call triangle(two_l1,two_j2,two_l3,delta)
      triangle_prod = triangle_prod*delta
      call triangle(two_l1,two_l2,two_j3,delta)
      triangle_prod = triangle_prod*delta
      

      zmin = max(d1, d2, d3, d4)
      zmax = min(n2, n3, n4)
      
      fac_sum = 0d0
      
      do z = zmin, zmax
          call logfac(z+1,lf)
          num = lf
C         num = logfac(z+1)
          call logfac(n2-z,lf)
          den = lf
          call logfac(n3-z, lf)
          den=den+lf
          call logfac(n4-z,lf)
          den=den+lf
          call logfac(z-d1, lf)
          den=den+lf
          call logfac(z-d2, lf)
          den=den+lf
          call logfac(z-d3, lf)
          den=den+lf
          call logfac(z-d4, lf)
          den=den+lf
C         den = logfac(n2-z) + logfac(n3-z) + logfac(n4-z) 
C    $        + logfac(z-d1) + logfac(z-d2) + logfac(z-d3) 
C    $        + logfac(z-d4)
          fac_sum = fac_sum + (-1) ** (z) * exp(num - den)
      end do
      
      sj = triangle_prod * fac_sum
      
      return
      end  
C     -------------------------------------------------------------
      subroutine triangle(two_j1, two_j2, two_j3, delta)  !Function from 2022 Oliver C. Gorton (ogorton@sdsu.edu) translated to subroutine
      
      ! Computes the triangle functions, typically denoted as
      !     \delta(abc) = (a+b-c)!(a-b+c)!(b+c-a)!/(a+b+c+1)!
      
      implicit none
      integer :: two_j1, two_j2, two_j3
      integer :: c1, c2, c3, c4
      real(kind=8) :: delta, lf
      
      delta = 0.0D0
      c1 = two_j1 + two_j2 - two_j3
      c2 = two_j1 - two_j2 + two_j3
      c3 =-two_j1 + two_j2 + two_j3
      c4 = two_j1 + two_j2 + two_j3
      
      if (c1.lt.0) return
      if (c2.lt.0) return
      if (c3.lt.0) return
      if (c4+1.lt.0) return
      
      if (mod(c1,2).ne.0) return ! \= means .ne.
      if (mod(c2,2).ne.0) return
      if (mod(c3,2).ne.0) return
      if (mod(c4,2).ne.0) return
      delta =0.0
      lf=0.0
      call logfac(c1/2,lf)
      delta = delta + lf
      call logfac(c2/2, lf)
      delta = delta + lf
      call logfac(c3/2, lf)
      delta = delta +lf
      call logfac(c4/2+1, lf)
      delta = delta - lf
      delta = 0.5d0*delta
      delta = exp(delta)
C     delta = exp(0.5d0*(logfac(c1/2)+logfac(c2/2)
C    $    +logfac(c3/2)-logfac(c4/2+1)))
      
      return
      
      end  
C     -------------------------------------------------------------
      subroutine logfac(n,lf) !Function from 2022 Oliver C. Gorton (ogorton@sdsu.edu) translated to subroutine
      
      ! Computes log(n!)
      ! Log-factorial is used to delay numerical overflow.
      
      implicit none
      integer :: n
      real(kind=8) :: lf
      integer :: i
      if (n<=0) then
          lf = 0d0
          return
      end if
      lf = 1d0
      do i = 1, n
          lf = lf * dble(i)
      end do
      lf = log(lf)
      return
      
      end  