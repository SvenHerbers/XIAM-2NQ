C----------------------------------------------------------------------
      subroutine calvjk_qdj2f12(jselect,gam,f,ib,f1select,evalv
     $ ,ovv,rotm,rott,tori,atot,qmv,ifittot,dfit,palc
     $ ,pali,npar,fistat,evhs)
C     calculation of the eigenvalues of one matrix with specified j,f,gam
C     the evalues are put in the field of dnv(1..ndata,Q_ENG,Q_UP/LO)
C     the deviations DE/DPi in dnv(1..ndata,DQ_ENG,Q_UP/LO(i))
C     fistat = 0 for regular calculation of Eigenvalues
C     fistat > 0  Eigenvalues for differential quotient
      implicit none
      include 'iam.fi'
      integer j, gam, f, ib, npar, fistat, is, f1, k, p, w, Temp                    
      integer startf1, endf1
      integer jselect
      integer f1select, D2
      integer sizef1,sizepre
      integer sj,ej,cm,cmb
      integer cm2
      integer minJ,maxJ
      integer occupied
      integer n
      real*8  hs(DIMQ2,DIMQ+DIMQ2,DIMTOT,DIMTOT)
      real*8  evhs(DIMQ2,DIMQ+DIMQ2,DIMTOT)                    
      real*8  h_2(DIMQ*DIMTOT,DIMQ*DIMTOT) 
      real*8  h_2NQ1(DIMQ*DIMTOT,DIMQ*DIMTOT)
      integer h_2_sizes(DIMQ2) 
      real*8  h_3(DIMQ2*DIMQ*DIMTOT,DIMQ2*DIMQ*DIMTOT) !h2024
      real*8  h_3NQ2(DIMQ2*DIMQ*DIMTOT,DIMQ2*DIMQ*DIMTOT) !For elements off-diagonal in F1 for HQ2
      real*8  evh_3(DIMQ2*DIMQ*DIMTOT)
      real*8  zr_3(DIMTOT*DIMQ*DIMQ2,DIMTOT*DIMQ*DIMQ2)
      real*8  zi_3(DIMTOT*DIMQ*DIMQ2,DIMTOT*DIMQ*DIMQ2)
      real*8  e_3(DIMQ2*DIMQ*DIMTOT),e2_3(DIMQ2*DIMQ*DIMTOT)
      real*8  tau_3(2,DIMQ2*DIMQ*DIMTOT)   
      real*8  signs(DIMTOT)
      real*8  vector(DIMTOT)
      integer qcase2
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP) 
      real*8  atot(DIMPAR,DIMVB)
      real*8  palc(DIMFIT,-1:DIMPLC) ! not used, but maybe later
      real*8  normis(DIMQ+DIMQ2)
      real*8  normis2(DIMQ2)
      real*8  normisF1(DIMTOT*DIMQ2*DIMQ,DIMQ2)     ! F1 quantum number assignment is not as easy as I thought. 
      real*8  nF1s(DIMQ2,DIMQ+DIMQ2,DIMTOT,DIMQ2)         ! same as normisF1 but after F1 block assignment.
      integer wF1s(DIMQ2)         ! I need to save the vector normisF1 as well as wF1s to write them to file for improved intensity prediciton.
              
      integer noJsinF1s(DIMQ+DIMQ2,DIMQ2)
      integer nJPPMinF1(DIMQ+DIMQ2,0:1,0:1,DIMQ2)
      real*8  collectnorms((DIMQ+DIMQ2)*(2*DIMJ+1),DIMQ2) ! A
      real*8  newnorms((DIMQ+DIMQ2)*(2*DIMJ+1))           ! B
      real*8  NormF1(DIMQ2)
      integer qcasesF1(DIMTOT*DIMQ2*DIMQ)
      real*8  normisEO(0:1)                                !Even or odd
      real*8  normisPM(0:1)                                !PM for wang assignment
      integer EOofI(DIMQ2*DIMQ*DIMTOT)                            ! separation by odd and even K quantum numbers.
      integer PMofI(DIMQ2*DIMQ*DIMTOT)                            ! separation by wangs gamma + or -
      integer indicesJ(DIMQ+DIMQ2,(DIMQ+DIMQ2)*(2*DIMJ+1))
      integer indicesJPPM(DIMQ+DIMQ2,0:1,0:1,(DIMQ+DIMQ2)*(DIMJ+1))
      integer indis((DIMQ+DIMQ2)*(2*DIMJ+1))
      integer hitsj(DIMQ+DIMQ2)
      integer hitsJPPM(DIMQ+DIMQ2,0:1,0:1)
      integer qcasesJ(DIMTOT*DIMQ2*DIMQ)
      integer pali(DIMFIT, 0:DIMPLC,2)       ! not used but maybe later.
      integer qmv(DIMV),ifittot(DIMPAR,DIMVB),dfit(DIMFIT)
      integer qmvs(DIMQ2,DIMQ,DIMV)
C     quantum numbers
      integer qmks(DIMQ2,DIMQ,DIMTOT,DIMQLP)
      integer qvks(DIMQ2,DIMQ,DIMTOT,Q_K:Q_V+DIMTOP) 
      integer qvs(DIMQ2,DIMQ,DIMTOT)
      integer qcase
C     work
      real*8  zrs(DIMQ2,DIMQ+DIMQ2,DIMTOT,DIMTOT)
      real*8  zis(DIMQ2,DIMQ+DIMQ2,DIMTOT,DIMTOT)
C      real*8  dedp(DIMPAR) ! not used but maybe later
      integer check2
      integer results
      integer results2
      integer ie,i,itop,ivr,ivc,it1,it2
      integer eused(DIMTOT), ierr
      integer ruse(DIMVV,DIMVV,DIMTOP)
      integer mycounters(DIMQ2,DIMQ+DIMQ2)
      integer hitsJP((DIMQ+DIMQ2),0:1)                         ! remove if not needed.
      integer indicesJP(DIMQ+DIMQ2,0:1,(DIMQ+DIMQ2)*(DIMJ+1)) ! remove if not needed.
      character*4 fnpre
C      character*6 fnpost
      logical masave
      logical complex
      integer myand
      external myand
      integer  mclock
      external mclock
      if (DIMJ.lt.(jselect+ctlint(C_SPIN)+ctlint(C_SPIN2))) then  
       write(*,*) 'Error DIMJ < J+2I'
       write(0,*) 'Error DIMJ < J+2I'
       stop
      end if
      masave=.false.
      
      if (ctlint(C_SPIN2).eq.1) then                                 !Some program termination conditions for input errors.
       write(0,*) "spin2 read as 1, but I_2 = 0.5",                  !Some program termination conditions for input errors.
     $  "is not implemented, stopping the program."                  !Some program termination conditions for input errors.
       stop                                                          !Some program termination conditions for input errors.
      end if                                                         !Some program termination conditions for input errors.
      if (ctlint(C_SPIN).eq.1) then                                  !Some program termination conditions for input errors.
       write(0,*) "spin read as 1, but I_2 = 0.5",                   !Some program termination conditions for input errors.
     $ "is not implemented with exact quadrupole coupling,",         !Some program termination conditions for input errors.
     $ "stopping the program."                                       !Some program termination conditions for input errors.
       if (ctlint(C_SPIN2).eq.0) then                                !Some program termination conditions for input errors.
        write(0,*) "You can try XIAM's spin rotation constants,",    !Some program termination conditions for input errors.
     $   "with a 0.5 nucleus using DWctrl 0 ."                       !Some program termination conditions for input errors.
       end if                                                        !Some program termination conditions for input errors.
       stop                                                          !Some program termination conditions for input errors.
      end if                                                         !Some program termination conditions for input errors.
      if (mod(f1select+ctlint(C_SPIN),2).ne.0) then                  !Some program termination conditions for input errors.
       write(0,*) "INPUT Error in 2f1=",f1select," vs 2I1=",         !Some program termination conditions for input errors.
     $ ctlint(C_SPIN),"both must be even or both must be odd",       !Some program termination conditions for input errors.
     $ "stopping the program."                                       !Some program termination conditions for input errors.
       stop                                                          !Some program termination conditions for input errors.
      end if                                                         !Some program termination conditions for input errors.
      if ((f.ne.-1).and.(mod(f+ctlint(C_SPIN)
     $                       +ctlint(C_SPIN2),2).ne.0)) then         !Some program termination conditions for input errors.
       write(0,*) "INPUT Error in 2f=",f," vs 2I1+2I2=",             !Some program termination conditions for input errors.
     $ ctlint(C_SPIN),"+",ctlint(C_SPIN),"both quantities must be",  !Some program termination conditions for input errors.
     $ "even or both must be odd, stopping the program."             !Some program termination conditions for input errors.
       stop                                                          !Some program termination conditions for input errors.
      end if                                                         !Some program termination conditions for input errors.
      if ((2*DIMJ).le.(f+ctlint(C_SPIN2))) then                      !An extra condition checking for DIMJ.
       if (ctlint(C_SPIN2).gt.0) then                                !An extra condition checking for DIMJ.
       write(0,*) "DIMJ too small for input f, stopping program"     !An extra condition checking for DIMJ.
        stop                                                         !An extra condition checking for DIMJ.
       end if                                                        !An extra condition checking for DIMJ.
      end if                                                         !An extra condition checking for DIMJ.
      
      
C      h_2=0.0
      hs=0.0! (:ctlint(C_SPIN),:(ctlint(C_SPIN)+ctlint(C_SPIN2)),:,:) dimension restriction didnt cause speedup in intialization.
      evhs=0.0
      normis=0.0
      mycounters=1
      wF1s=-1
      nF1s=0.0

      startf1=f-ctlint(C_SPIN2)
      endf1=(f+ctlint(C_SPIN2))
      if (startf1.lt.0) then
        startf1=mod(ctlint(C_SPIN),2) ! if less than 0, starting f1 is put to be 0 or 0.5 (*2=1) depending on wether I1 is a half integer or full integer
      end if
      !if spin2 is equal to 0 this should return f1select.
      !if spin2 is 0 then f=f1.
      if (f.lt.ctlint(C_SPIN2))then ! if F is smaller than I, then there are F1 states that can not exist
       if ((abs(f+startf1)).lt.ctlint(C_SPIN2))then
       startf1=ctlint(C_SPIN2)-f
       end if
      end if

      if (f.eq.-1) then !only place where f1select should be needed
       startf1=f1select !only place where f1select should be needed
       endf1=f1select   !only place where f1select should be needed
      end if            !only place where f1select should be needed
      
      minJ=(startf1-ctlint(C_SPIN))/2 ! these will be needed to get the vector components right later.
      if (minJ.lt.0) then
       minJ=0
      end if
       maxJ=(endf1+ctlint(C_SPIN))/2 ! these will be needed to get the vector components right later.
      

      !!! for non existing J/F states with J>jselect no matrix has to be built up
      if (ctlint(C_EVAL).gt.3)   masave=.true.

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1))
     $     write(*,'(A,5I3)')
     $     'starting with J,S,B,F Fit_stat=',jselect,gam,ib,
     $     f,fistat
      fnpre='xiam'
      if (gam.eq.0) ctlint(C_NTOP)=0
      complex=.true.!.false.

       do itop=1, ctlint(C_NTOP)
        do ivr=1, size(S_VV)
          do ivc=1, size(S_VV)
            ruse(ivr,ivc,itop)=1
          end do
        end do
      end do

      do it1=1, ctlint(C_NTOP)
        do it2=1, ctlint(C_NTOP)
          if (it2.ne.it1) then
            do ivr=1, size(S_VV)
              do ivc=1, size(S_VV)
                if (qvv(ivr,it2,ib).ne.qvv(ivc,it2,ib)) then    
                  ruse(ivr,ivc,it1)=0
                end if
              end do
            end do
          end if
        end do
      end do


C------Construction of Htot Starts
C------Construction of Htot Starts  
C------Construction of Htot Starts
      noJsinF1s=0
      nJPPMinF1=0
      h_2_sizes=0
      h_3NQ2=0.0
      h_3 =0.0

      do f1=startf1,endf1,2 ! uses a step size of 2
       
       cm=(f1-startf1)/2!count matrices
C       h_2=0.0 ! is reinitilaized in nqvmat_ir anyway
       wF1s(cm+1)=f1 ! need to save the f1 to written vectors and matrices together with F1 contributions for later evaluation in intensity predictions
       sj=(2*jselect-(f1-ctlint(C_SPIN)))/2!startj offset
       ej=((f1+ctlint(C_SPIN))-2*jselect)/2!endj offset
       if ((jselect-sj).lt.0) sj=jselect   ! the start incidces and end indices for the F1 matrix depend on the f1 used

      if (f1.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs(f1+2*(jselect-sj)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       sj=sj-(abs(abs(f1-2*(jselect-sj))-ctlint(C_SPIN))/2)
       end if 
      end if

       
       if (f1select.eq.-1) then
        sj=0
        ej=0
       end if
       
       call nqcmat_ir(jselect,gam,f1,ib,evalv,ovv,rotm,rott !This sets up the semi-rigid rotor and internal rotation matrix elements for the various J blocks in the F matrix
     $  ,tori,atot,qmv,h_2,qvks(cm+1,:,:,:),qmks(cm+1,:,:,:)
     $  ,qmvs(cm+1,:,:),qvs(cm+1,:,:),ruse)  
     
     
       call nqcmat_NQ1(jselect,f1,ib,h_2NQ1,atot,ctlint(C_SPIN)) ! This builds the Matrix elements for HQ1, these are all diagonal in F1.
       h_2=h_2+h_2NQ1 ! adding first nucleus hamiltonian 
       size(s_h)=0  
       do j=jselect-sj,jselect+ej ! the sizes depend on the specific F1 also.
         cm2=j-minJ
         noJsinF1s(cm2+1,cm+1)=noJsinF1s(cm2+1,cm+1)+2*j+1  !will later be used in quantum number assignments...
         nJPPMinF1(cm2+1,0,0,cm+1)=nJPPMinF1(cm2+1,0,0,cm+1)+1
         do k=1,j
          if (mod(k,2).eq.0) then
           nJPPMinF1(cm2+1,0,0,cm+1)=nJPPMinF1(cm2+1,0,0,cm+1)+1 ! one with odd one with even wang factor
           nJPPMinF1(cm2+1,0,1,cm+1)=nJPPMinF1(cm2+1,0,1,cm+1)+1 ! one with odd one with even wang factor
          else
           nJPPMinF1(cm2+1,1,0,cm+1)=nJPPMinF1(cm2+1,1,0,cm+1)+1 ! one with odd one with even wang factor
           nJPPMinF1(cm2+1,1,1,cm+1)=nJPPMinF1(cm2+1,1,1,cm+1)+1 ! one with odd one with even wang factor
          end if
         end do
         if (j.ge.0) size(s_h)=size(s_h)+2*j+1             ! 
       end do
       occupied=sum(h_2_sizes)
       h_2_sizes(cm+1)=size(s_h)
      
C      !Matrix elements diagonal in J for second nucleus
       if (ctlint(C_SPIN2).ne.0) then
        call nqcmat_NQ2df1eqn(jselect,f1,f,ib,atot,h_3NQ2,occupied,0) ! Matrix HQ2, diagonal in F1
        if (f1.le.endf1-1) then
         call nqcmat_NQ2df1eqn(jselect,f1,f,ib,atot,h_3NQ2,occupied,1) ! Matrix HQ2, offdiagonal in F1 by +-1
        end if
        if (f1.le.endf1-2) then
         call nqcmat_NQ2df1eqn(jselect,f1,f,ib,atot,h_3NQ2,occupied,2) ! Matrix HQ2, offdiagonal in F1 by +-2
        end if
       end if
C      !Off diagonal elements delta F1 pm 1

      h_3(1+occupied:occupied+size(s_h)
     $             ,1+occupied:occupied+size(s_h))=h_2(1:size(s_h)
     $             ,1:size(s_h)) ! filling up the diagonals.
      end do

      
C     Updating size

      size(s_h)=occupied+size(s_h)
      
C     Adding second nucleus to hamiltonian
      D2=DIMQ2*DIMQ*DIMTOT 
C      D2=size(S_H) ! replacing D2 with size(S_H) leads to quite some speed up but will cause stack overflow for iodotoluene for some reason?
      h_3(1:size(S_H),1:size(S_H)) = 
     $   h_3(1:size(S_H),1:size(S_H)) 
     $ + h_3NQ2(1:size(S_H),1:size(S_H)) ! Adding NQ2 matrix to total matrix.

C       The diagonalization can not handle zero rows + columns, so the following offset was added.
      do i=1,size(S_H)                     !
      if (abs(h_3(i,i)).le.1.0E-14) then   !The diagonalization can not handle zero rows + columns, so this offset was added.
       h_3(i,i) = 5.0E-14                  !The diagonalization can not handle zero rows + columns, so this offset was added.
      end if                               !The diagonalization can not handle zero rows + columns, so this offset was added.
      end do                               !The diagonalization can not handle zero rows + columns, so this offset was added.
C------Construction of Htot finished
C------Construction of Htot finished   
C------Construction of Htot finished   
     
      

C       e_3(1:D2)=0.0      ! These are outputs of htrid3, intialization should be not required.
C       e2_3(1:D2)=0.0     ! These are outputs of htrid3, intialization should be not required.
C       evh_3(1:D2)=0.0    ! These are outputs of htrid3, intialization should be not required.
C       tau_3=0.0          ! These are outputs of htrid3, intialization should be not required.
       zr_3(1:D2,1:D2)=0.0
      do i=1, D2 !initialize zr_2 for diagonalization routine
       zr_3(i,i)=1.0
      end do

      call htrid3 (D2,size(S_H),h_3(:D2,:D2),evh_3(:D2),e_3(:D2), ! requires zr to be a unit matrix as input.
     $  e2_3(:D2),tau_3(:,:D2))
      ierr=0.0
      call tql2 (D2,size(S_H),evh_3(:D2),e_3(:D2),zr_3(:D2,:D2),ierr)

      if (ierr.ne.0) then
          write (*,'(a,i5)') 'Error in tql2 ',ierr
          stop
      endif
      call htrib3 (D2,size(S_H),h_3(:D2,:D2),tau_3(:,:D2)
     $           ,size(S_H),zr_3(:D2,:D2),zi_3(:D2,:D2))
cc     sort eigenvalues in **ascending** order
      call heigsrt(evh_3(:D2),zr_3(:D2,:D2),zi_3(:D2,:D2),size(S_H),D2) 

      do i=1, size(S_H)
        do ie=1, size(S_H)
            h_3(i,ie)=sign(dsqrt(zr_3(i,ie)**2+zi_3(i,ie)**2)
     $           ,zr_3(i,ie)+zi_3(i,ie))              
        end do
      end do
      

C------Diagnoalization finished
C------Diagnoalization finished
C------Diagnoalization finished

        cm2=0                                     
        hitsj=0
        hitsJPPM=0
        hitsJP=0
        do i = 1, size(S_H)   ! here goes one size... 
          normis=0   !normis for j identification.
          normis2=0  ! normis2 for f1 identification
          normisEO=0.0
          normisPM=0.0   
          do f1= startf1, endf1,2         
           cm2=(f1-startf1)/2
                 sizef1=h_2_sizes(cm2+1)  
                 sizepre=sum(h_2_sizes(1:cm2+1))-sizef1                  
                 if (f1.ge.0) then

                    normis2(cm2+1)=
     $                sum((h_3(sizepre+1:sizepre+sizef1,i)**2)) ! for f1 quantum number assignment.
                 end if
           normisF1(i,cm2+1)=normis2(cm2+1) !saving F1 components for all lines.
           wF1s(cm2+1)=f1
           
           sj=(2*jselect-(f1-ctlint(C_SPIN)))/2!startj offset
           ej=((f1+ctlint(C_SPIN))-2*jselect)/2!endj offset
           if ((jselect-sj).lt.0) sj=jselect
           if (f1.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
            if (abs(f1+2*(jselect-sj)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
            sj=sj-(abs(abs(f1-2*(jselect-sj))-ctlint(C_SPIN))/2)
            end if 
           end if
           
           check2=0
           
           do j=jselect-sj,jselect+ej
             cm=j-minJ ! here sj and ej differ for the different f1. To add all the components properly they have to be referenced relative to a constant value, which here is the minimum J
C             write(0,*) h_3(sizepre+check2+1+j-3:sizepre+check2+1+j+3,i)
             if (j.ge.0) then
                    
                 normis(cm+1)=normis(cm+1)
     $                +sum((h_3(sizepre+check2+1:sizepre+check2+2*j+1
     $                 ,i)**2))                                      !sizepre+check2+1:sizepre+check2+2*j+1 it climbs through all j contributions by moving check2 by 2*j+1 each time.
                 if (normisPM(0).lt.0) then !Sign flip so PM(0) always starts positive.
                 normisPM(1)=-normisPM(1)
                 end if
                 normisPM(0)=abs(normisPM(0))
     $                                 +(h_3(sizepre+check2+1+j,i))**2  !adding J component at K=0, K=0 means Plus and Even so PM(0) and EO(0) get this added.
                 normisEO(0)=normisEO(0)+h_3(sizepre+check2+1+j,i)**2
                 do n=1, j
                 if ((h_3(sizepre+check2+1+j,i))**2             ! case Sign taken from vectors.
     $                          .ge.abs(normisPM(0))) then
                 
                 normisPM(0)=sign(normisPM(0)
     $                 ,h_3(sizepre+check2+1+j+n,i))
     $                 +sign(h_3(sizepre+check2+1+j+n,i)**2
     $                 ,h_3(sizepre+check2+1+j+n,i))
                 normisPM(1)=normisPM(1)
     $                     +sign(h_3(sizepre+check2+1+j-n,i)**2
     $                     ,h_3(sizepre+check2+1+j-n,i))

                 else
                  if (normisPM(0)*h_3(sizepre+check2+1+j+n,i).lt.0) then  !if sign change on one expected, also sign change on other
                  normisPM(1)=normisPM(1)
     $                       -sign(h_3(sizepre+check2+1+j-n,i)**2
     $                       ,h_3(sizepre+check2+1+j-n,i))
                  else
                  normisPM(1)=normisPM(1)
     $                       +sign(h_3(sizepre+check2+1+j-n,i)**2
     $                       ,h_3(sizepre+check2+1+j-n,i))
                  end if
                 
                 normisPM(0)=normisPM(0)
     $                  +sign(h_3(sizepre+check2+1+j+n,i)**2
     $                  ,normisPM(0))    ! sign taken from normisPM
                 end if
                  if (mod(n,2).eq.0) then
                  normisEO(0)=normisEO(0)+h_3(sizepre+check2+1+j+n,i)**2 !Adding +K
     $                      +h_3(sizepre+check2+1+j-n,i)**2              !Adding -K
                  else
                  normisEO(1)=normisEO(1)+h_3(sizepre+check2+1+j+n,i)**2 !Adding +K
     $                      +h_3(sizepre+check2+1+j-n,i)**2              !Adding -K
                  end if
                 end do
                 
                 check2=check2+2*j+1
             end if
           end do !j loop
          end do !f1 loop
          
          !now the maximums are found by maxloc, this is after the loops of j and f1, but the i loop is still open.
         if (normisEO(0).gt.normisEO(1)) then
          EOofI(i)=0
         else
          EOofI(i)=1
         end if
         
         
         if (abs(normisPM(1)).gt.0.9) then ! exception where wang does not work anymore.
           PMofI(i)=1
         else
         if ((normisPM(0)*normisPM(1)).lt.0) then !wangs gamma separation, see also Gordy equation 7.22

          if (abs(normisPM(1)).gt.0.1) then ! for all but K=0 the normisPM(0) and (1) should be about 0.5 in their absolutes., for K=0 normisPM(1) must be very close to zero.
           PMofI(i)=1
          else
           PMofI(i)=0
          end if
         else
          PMofI(i)=0
         end if
         end if
         
C         if (gam.eq.1) then
C         write(0,*) normisPM(0), normisPM(1) ! 
C         end if
C         if (abs(normisPM(1)).ge.0.9) then
C          qcase=maxloc(normis,1)-1+minJ
C          write(0,*) qcase
C          write(0,*) startf1,endf1
C          check2=0
C          do n=0,qcase-1
C          check2=check2+2*i+1
C          end do
C          do f1= startf1, endf1,2         
C           cm2=(f1-startf1)/2
C                 sizef1=h_2_sizes(cm2+1)  
C                 sizepre=sum(h_2_sizes(1:cm2+1))-sizef1                  
CC          write(0,*) h_3(sizepre+check2+1:
CC     $         sizepre+check2+2*qcase+1,i)
C         end do
C         if (qcase.eq.7) then
C         write(0,*) h_3(49+11+13:,i)
C         end if
C         end if

         qcase2=maxloc(normis2,1) !F1 for testing, not assignment
         sizef1=h_2_sizes(qcase2) 
         sizepre=sum(h_2_sizes(1:qcase2))-sizef1 
         
         qcase=maxloc(normis,1) !J - actual assignment, should work as long as J is a nearly good quantum number.
         qcasesJ(i)=qcase              !records the q case of line i

         hitsj(qcase)=hitsj(qcase)+1   !counts the number of lines assigned to each j quantum number
         hitsJPPM(qcase,EOofI(i),PMofI(i))=
     $                hitsJPPM(qcase,EOofI(i),PMofI(i))+1
         hitsJP(qcase,EOofI(i))=
     $                hitsJP(qcase,EOofI(i))+1
         
         indicesJ(qcase,hitsj(qcase))=i!saves the line index of the nth line of the J block with corresponding to q case.
         indicesJPPM(qcase,EOofI(i),PMofI(i), 
     $         hitsJPPM(qcase,EOofI(i),PMofI(i)))=i
         indicesJP(qcase,EOofI(i),hitsJP(qcase,EOofI(i)))=i
        
          ! qcase2 gives the f1 quantum numbers, qcase gives the J quantum number.
         if (normis2(qcase2).le.0.05) then
          write(0,*) 'Vector f1 not found'
          qcase2=-1          !e.g. for empty matrix entries where all normis2 are zero
         end if
          if (normis(qcase).le.0.05) then !e.g. for empty matrix entries where all normis are zero
           write(0,*) 'Vector j not found'
           qcase=-1
          end if
        end do   !end for column i of h_3  
        
        
        
        
        
        
        
        if ((ctlint(C_SORT).eq.1).or.
     $          ((ctlint(C_SORT).eq.3).and.
     $          ((gam.eq.1).or.(gam.eq.0)))) then !gam=0 und gam=1 sind S0 oder S1, also rigid rotor or S1 !Sort! and 3 cause A species to be assigned with E/O filter and Wangs gamma filter, For Sort 1 also E species is assigned this way.
C        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IF CASES FOR F1 ASSIGNMENT SCHEME STARTS
        do j=minJ,maxJ       ! for each J block < this works perfectly for A species
          cm=j-minJ
          NormF1=0.0
          collectnorms=0.0
         do p=0,1
         do w=0,1
          do i=1, hitsJPPM(cm+1,p,w)!for each line in the block
           NormF1=normisF1(indicesJPPM(cm+1,p,w,i),:) !read the F1 vector sums saved to normisF1
           collectnorms(i,:)=NormF1            !X
          end do
          do f1=startf1,endf1,2                 ! for all F1 quantum numbers
            cm2=(f1-startf1)/2
            newnorms=collectnorms(:,cm2+1)              !load vector contributions to F1
            call argsort(newnorms,
     $         indis,2*(DIMJ+DIMQ+DIMQ2)+1,hitsJPPM(cm+1,p,w)) ! sort
            do n=1,nJPPMinF1(cm+1,p,w,cm2+1)
               qcasesF1(indicesJPPM(cm+1,p,w,indis(n)))=cm2+1 ! intial F1 assignment
               collectnorms(indis(n),:)=0! this line is already assigned, and is removed from subsequent assignment processes.
            end do    ! end n        
          end do ! end f1
         end do ! end w
         end do ! end p
        end do  ! end J
        
        
        
        do i=1,size(S_H) ! Final F1 assignment in double do loop
          j=0
          do while (j.lt.size(S_H))
           j=j+1
           if (j.ne.i) then
            if (qcasesJ(i).eq.qcasesJ(j)) then
            if (EOofI(i).eq.EOofI(j)) then
            if (PMofI(i).eq.PMofI(j)) then
             if ((normisF1(i,qcasesF1(j))+normisF1(j,qcasesF1(i))).gt.
     $           (normisF1(i,qcasesF1(i))+normisF1(j,qcasesF1(j)))) then
             Temp=qcasesF1(i)
             qcasesF1(i)=qcasesF1(j)
             qcasesF1(j)=Temp
             j=0 !resets the loop if a swap was made to check again.
             end if
            end if
            end if
            end if
           end if
          end do
        end do
C        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ELSE CASES FOR F1 ASSIGNMENT SCHEME STARTS
        else ! if sort 2, or sort 3 and E species
C        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ELSE CASES FOR F1 ASSIGNMENT SCHEME STARTS
        do j=minJ,maxJ       ! for each J block
          cm=j-minJ
          NormF1=0.0
          collectnorms=0.0
         do p=0,1
          do i=1, hitsJP(cm+1,p)!for each line in the block
           NormF1=normisF1(indicesJP(cm+1,p,i),:) !read the F1 vector sums saved to normisF1
           collectnorms(i,:)=NormF1            !X
          end do
          do f1=startf1,endf1,2                 ! for all F1 quantum numbers
            cm2=(f1-startf1)/2
            newnorms=collectnorms(:,cm2+1)              !load vector contributions to F1
            call argsort(newnorms,
     $         indis,2*(DIMJ+DIMQ+DIMQ2)+1,hitsJP(cm+1,p)) ! sort
            do n=1,nJPPMinF1(cm+1,p,0,cm2+1)+nJPPMinF1(cm+1,p,1,cm2+1)
               qcasesF1(indicesJP(cm+1,p,indis(n)))=cm2+1
               collectnorms(indis(n),:)=0! this line is already assigned, and is removed from subsequent assignment processes.
            end do    ! end n        
          end do ! end f1
         end do ! end p
        end do  ! end J
        
        
        do i=1,size(S_H) ! Final F1 assignment in double do loop
          j=0
          do while (j.lt.size(S_H))
           j=j+1
           if (j.ne.i) then
            if (qcasesJ(i).eq.qcasesJ(j)) then
            if (EOofI(i).eq.EOofI(j)) then
            if (.true.) then !(PMofI(i).eq.PMofI(j)) then ! here the pm filter is turned of.
             if ((normisF1(i,qcasesF1(j))+normisF1(j,qcasesF1(i))).gt.
     $           (normisF1(i,qcasesF1(i))+normisF1(j,qcasesF1(j)))) then
             Temp=qcasesF1(i)
             qcasesF1(i)=qcasesF1(j)
             qcasesF1(j)=Temp
             j=0 !resets the loop if a swap was made to check again.
             end if
            end if
            end if
            end if
           end if
          end do
        end do
C        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END IF CASES FOR F1 ASSIGNMENT SCHEME STARTS
        end if
C        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END IF CASES FOR F1 ASSIGNMENT SCHEME STARTS



        do i = 1, size(S_H)   ! 
         qcase2=qcasesF1(i) !'loading' qcase for F1 from new assignment process
         qcase=qcasesJ(i) !'loading' qcase for J, still assuming that J remains a nearly good quantum number.
         
         do f1= startf1, endf1,2  !another f loop, to add the corresponding j blocks on one each other, and to save the energy as well as the final vectors.   
           
           cm2=(f1-startf1)/2
           sizef1=h_2_sizes(cm2+1)  
           sizepre=sum(h_2_sizes(1:cm2+1))-sizef1
           
           
           sj=(2*jselect-(f1-ctlint(C_SPIN)))/2!startj offset 
           ej=((f1+ctlint(C_SPIN))-2*jselect)/2!endj offset
           if ((jselect-sj).lt.0) sj=jselect   
           if (f1.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
            if (abs(f1+2*(jselect-sj)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
            sj=sj-(abs(abs(f1-2*(jselect-sj))-ctlint(C_SPIN))/2)
            end if 
           end if
           check2=0

           do j=jselect-sj,jselect+ej              !another j loop to find the propper indices using check2.
             
             
             cm=j-minJ 
             if (qcase.eq.cm+1) then 
             vector=0.0               
             
             vector(1:2*j+1)=h_3(sizepre+check2+1
     $                       :sizepre+check2+2*j+1,i) ! last mention of i.
              hs(qcase2,qcase,1:2*j+1,mycounters(qcase2,qcase))=
     $        hs(qcase2,qcase,1:2*j+1,mycounters(qcase2,qcase))+  !these will add for each j found in the various matrices. Now we have to make sure the signs match.
     $        vector(1:2*j+1)**2                                       !Here it is still the vectors squared.
              if (qcase2.eq.cm2+1) then 
                results=check2+1                                  ! I will also use this case to get upper and lower index for sign vector reading.
                results2=check2+2*j+1                             ! I will also use this case to get upper and lower index for sign vector reading.
                signs=h_3(sizepre+1:sizepre+sizef1,i)             ! signs are saved.
              end if
             end if 
            check2=check2+2*j+1 ! to step the j indices in the vector component adding.
           end do ! close j loop for vector components.
          end do ! close f1 loop. We have now added all the vector components and saved the energy
          check2=results2-results+1! hamiltonian size
         nF1s(qcase2,qcase,mycounters(qcase2,qcase),:)=normisF1(i,:)
         
         evhs(qcase2,qcase,mycounters(qcase2,qcase))=evh_3(i)  ! if also the f1 quantum number is matched, save energy.
          hs(qcase2,qcase,1:check2,mycounters(qcase2,qcase))=
     $      SQRT(hs(qcase2,qcase,1:check2,mycounters(qcase2,qcase)))! taking the square root of the square sum of vectors. mycounters starts at one, after the previous do loops it should be stepped by one.
          hs(qcase2,qcase,1:check2,mycounters(qcase2,qcase))=
     $     SIGN(hs(qcase2,qcase,1:check2,mycounters(qcase2,qcase)) !  
     $           ,signs(results:results2)) !now the signs were restored following the signs of the largest contributing f1 !normalization should not be required since all components are used.
         mycounters(qcase2,qcase)=mycounters(qcase2,qcase)+1   !counter start at 1
         end do   !end for column i of h_3   
         
      if ((ctlint(C_INTS).gt.0).and.(fistat.eq.0)) then
      zis=0.0
      zrs=hs   
      end if
      
      do f1= startf1, endf1,2 
      cm2=(f1-startf1)/2
      sj=(2*jselect-(f1-ctlint(C_SPIN)))/2!startj offset 
      ej=((f1+ctlint(C_SPIN))-2*jselect)/2!endj offset
      if ((jselect-sj).lt.0) sj=jselect
      if (f1.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs(f1+2*(jselect-sj)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       sj=sj-(abs(abs(f1-2*(jselect-sj))-ctlint(C_SPIN))/2)
       end if 
      end if

      do j=jselect-sj,jselect+ej
      cm=j-minJ
      cmb=j-(jselect-sj)
      
      size(S_H)= 2*(j)+1  ! added because it is used in wvrec,
      size(S_K)= 2*(j)+1 ! saved S_K added because S_K will be used in other routines
      
      
      
      call assgn(j,gam,f,ib,f1,hs(cm2+1,cm+1,:,:),evhs(cm2+1,cm+1,:)
     $           ,qvks(cm2+1,cmb+1,:,:),qmvs(cm2+1,cmb+1,:),
     $           qmks(cm2+1,cmb+1,:,:),qvs(cm2+1,cmb+1,:),fistat)
      call esave(j,gam,f,ib,f1,qmks(cm2+1,cmb+1,:,:)
     $           ,evhs(cm2+1,cm+1,:),eused)

      if ((ctlint(C_INTS).gt.0).and.(fistat.eq.0)) then
      
      call wrvecQ(zrs(cm2+1,cm+1,:,:),zis(cm2+1,cm+1,:,:)
     $           ,evhs(cm2+1,cm+1,:),j,gam,f,ib,f1
     $           ,qvks(cm2+1,cmb+1,:,:),qmks(cm2+1,cmb+1,:,:)
     $           ,wF1s,nF1s(cm2+1,cm+1,:,:))
      
      
      
C        call wrvec(zrs(cm2+1,cm+1,:,:),zis(cm2+1,cm+1,:,:)
C     $           ,evhs(cm2+1,cm+1,:),j,gam,f,ib,f1
C     $           ,qvks(cm2+1,cmb+1,:,:),qmks(cm2+1,cmb+1,:,:))
      end if

       end do !js
       end do !f1s
        
      return
       end 
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      subroutine calvjk(j,gam,f,ib,f1,h,evalv,ovv,rotm,rott,tori
     $     ,a,qmv,ifit,dfit,palc,pali,npar,fistat,evh)
C     calculation of the eigenvalues of one matrix with specified j,f,gam
C     the evalues are put in the field of dnv(1..ndata,Q_ENG,Q_UP/LO)
C     the deviations DE/DPi in dnv(1..ndata,DQ_ENG,Q_UP/LO(i))
C     fistat = 0 for regular calculation of Eigenvalues
C     fistat > 0  Eigenvalues for differential quotient

      implicit none
      include 'iam.fi'
      integer j, gam, f, ib, npar, fistat, is, f1
      real*8  h(DIMTOT,DIMTOT),evh(DIMTOT)
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      real*8  a(DIMPAR)
      real*8  palc(DIMFIT,-1:DIMPLC)
      
      integer pali(DIMFIT, 0:DIMPLC,2)
      integer qmv(DIMV),ifit(DIMPAR),dfit(DIMFIT)

C     quantum numbers
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP), qmk(DIMTOT,DIMQLP)
      integer qv(DIMTOT)
C     work
      real*8  e(DIMTOT),e2(DIMTOT),tau(2,DIMTOT)
      real*8  zr(DIMTOT,DIMTOT),zi(DIMTOT,DIMTOT)
      real*8  dedp(DIMPAR)
      integer id,ie,i,iv,ik,itop,ivr,ivc,ir,ic,it1,it2,ikr,ikc
      integer eused(DIMTOT), usert,ierr
      integer ruse(DIMVV,DIMVV,DIMTOP)
      character*30 fmtstr
      character*40 fmtstr2
      character*4 fnpre
      character*6 fnpost
      real*8 tt
      logical masave
      logical complex
      integer myand
	  
      external myand

      integer  mclock,t1,t2
      external mclock
      masave=.false.
      if (ctlint(C_EVAL).gt.3)   masave=.true.

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1))
     $     write(*,'(A,5I3)')
     $     'starting with J,S,B,F Fit_stat=',J,gam,ib,f,fistat
      fnpre='xiam'
      if (gam.eq.0) ctlint(C_NTOP)=0
      complex=.false.
      do i=1,ctlint(C_NTOP)
        if (a(PI_GAMA+(i-1)*DIMPIR+DIMPRR).ne.0.0) complex=.true.
      end do
      usert = 1
      if ((a(P_QYZ).ne.0.0).or.(a(P_QXY).ne.0.0)) complex=.true.

      size(S_K)=2*j+1
C     initialize the quantum no.s qvk
      i=0
      do iv=1, size(S_VV)
        do ik=1, size(S_K)
          i=i+1
          qvk(i,Q_K) =ik-j-1
          qv(i)=iv
          do itop=1, ctlint(C_NTOP)
            qvk(i,Q_V+itop)=qvv(iv,itop,ib)
          end do
        end do
      end do
      size(S_H)=i
      do ivr=1, size(S_H)
        do ivc=1, size(S_H)
          h(ivr,ivc)=0.0
        end do
      end do

      do itop=1, ctlint(C_NTOP)
        do ivr=1, size(S_VV)
          do ivc=1, size(S_VV)
            ruse(ivr,ivc,itop)=1
          end do
        end do
      end do

      do it1=1, ctlint(C_NTOP)
        do it2=1, ctlint(C_NTOP)
          if (it2.ne.it1) then
            do ivr=1, size(S_VV)
              do ivc=1, size(S_VV)
                if (qvv(ivr,it2,ib).ne.qvv(ivc,it2,ib)) then    
                  ruse(ivr,ivc,it1)=0
                end if
              end do
            end do
          end if
        end do
      end do

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1)) then
        write(fmtstr,'(A,I1,A,I2,A)') '(',ctlint(C_NTOP),'I2,A,'
     $       ,size(S_VV),'I3)'
        write(*,*) fmtstr
        do itop=1, ctlint(C_NTOP)
          write(*,'(A,I3)') ' ruse array for top', itop 
          do ivr=1, size(S_VV)
            write(*,fmtstr)(qvv(ivr,it2,ib),it2=1,ctlint(C_NTOP))
     $           ,' | ',(ruse(ivr,ivc,itop),ivc=1, size(S_VV))
          end do
        end do
      end if 

      if (masave) then
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'t'//fnpost,status='unknown')
          write(55,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              write(55,*)
     $             ((rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,itop)
     $             *tori(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc
     $             ,gamma(gam,itop),itop)
     $             ,ikc=1,size(S_K))
     $             ,ivc=1, size(S_V+itop))
            end do
          end do
          close(55)
        end do
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(56,file=fnpre//'p'//fnpost,status='unknown')
          write(56,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              do ivc=1, size(S_V+itop)
                do ikc=1, size(S_K)
                  if (qvk(ikr,Q_K).eq.qvk(ikc,Q_K)) then
                    write(56,'(D22.14,$)')
     $                   ovv(ivr,ivc,PM_PI,gamma(gam,itop)
     $                   ,qvk(ikr,Q_K),itop)
                  else
                    write(56,'(D22.14,$)') 0.0
                  end if
                end do
              end do
              write(56,*)
            end do
          end do
          close(56)
        end do
        
        do it1=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         it1,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'r'//fnpost,status='unknown')
          write(55,*) size(S_H),0,0
          do ir=1, size(S_H)
            do ic=1, size(S_H)
              tt=1.0d0
              do it2=1,ctlint(C_NTOP)
                tt=tt*tori(qvk(ir,Q_K),    qvk(ic,Q_K),
     $                   qvk(ir,Q_V+it2),qvk(ic,Q_V+it2),
     $                   gamma(gam,it2),it2)
              end do
              write(55,'(D22.14,$)')
     $             rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,it1)
     $             *tt
            end do
            write(55,*)
          end do
          close(55)
        end do
      end if

      if (gam.ne.0) then
C      t1=mclock()
C     build the rotated D^{T} E_{K v \sigma} D
C        if (usert.eq.0) then
C          call bld1vjk(j,gam,f,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori)
C        else
        call bld2vjk(j,gam,f1
     $         ,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,complex)
C        end if

c      write(*,*) 'bld2vjk',mclock()-t1

      else
        if (size(S_VV).gt.1) stop ' size vv > 1 for rigid rotor!'
      end if
C      t1=mclock()
      call addrig_old(j,gam,f1
     $     ,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,complex)
        
c      write(*,*) 'addrig',mclock()-t1
      if (myand(ctlint(C_PRI),AP_MH).ne.0) then
        write(*,*) ' H_tot '
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
C---
C      t1=mclock()
      if (complex) then
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            zr(ir,ic)=0.0
          end do
          zr(ir,ir)=1.0
        end do
        ierr=0
        call htrid3 (DIMTOT,size(s_h),h,evh,e,e2,tau)
        call tql2 (DIMTOT,size(s_h),evh,e,zr,ierr)
        if (ierr.ne.0) then
          write (*,'(a,i5)') 'Error in tql2 ',ierr
          stop
        endif
        call htrib3 (DIMTOT,size(S_H),h,tau,size(S_H),zr,zi)

C     sort eigenvalues in **ascending** order
        call heigsrt(evh,zr,zi,size(S_H),DIMTOT)

        do i=1, size(S_H)
          do ie=1, size(S_H)
            h(i,ie)=sign(dsqrt(zr(i,ie)**2+zi(i,ie)**2)
     $           ,zr(i,ie)+zi(i,ie))              
          end do
        end do
      else
        write(fnpost,'(A,I1,A,I1,I1)')
     $       's',gam,'.j',int(j/10),mod(j,10)
        if (masave) then
          open(55,file=fnpre//'h'//fnpost,status='unknown')
          write(55,*) size(S_H),j,gam
          do i=1, size(S_H)
            write(55,*)(h(i,ie),ie=1,size(S_H))
          end do
          close(55)
        end if
 
        call hdiag(DIMTOT,size(S_H),h,evh,e,ierr)
        if (ierr.ne.0) then
          write (*,'(a,i5)') 'Error in hdiag ',ierr
          stop
        endif
C     sort eigenvalues in **ascending** order
        call eigsrt(evh,h,size(S_H),DIMTOT)
C---
        if (masave) then
          open(55,file=fnpre//'e'//fnpost,status='unknown')
          write(55,*) size(S_H),j,gam
          do i=1, size(S_H)
            write(55,*) evh(i)
          end do
          close(55)
C---  
          open(55,file=fnpre//'v'//fnpost,status='unknown')
          write(55,*) size(S_H),j,gam
          do i=1, size(S_H)
            write(55,*)(h(i,ie),ie=1,size(S_H))
          end do
          close(55)
        end if
        do i=1, size(S_H)
          do ie=1, size(S_H)
            zi(i,ie)=0.0
            zr(i,ie)=h(i,ie)
          end do
        end do        
      end if
c      write(*,*) 'diag   ',mclock()-t1

      if (myand(ctlint(C_PRI),AP_EH).ne.0) then
        write(*,*) ' Eigenvectors'
        do ir=1, size(S_H)
          do ic=1, size(S_H)
            if (h(ir,ic).lt.1000.0) then
              write(*,'(2F10.6,A,$)') zr(ir,ic),zi(ir,ic),'; '
            else
              write(*,'(2F10.2,A,$)') zr(ir,ic),zi(ir,ic),'; '
            end if
          end do
          write(*,*)
        end do
        write(*,*)
        write(*,*) ' Eigenvalues'
        do ir=1, size(S_H)
          write(*,'(F20.8,A,$)') evh(ir),'; '
        end do
        write(*,*)
      end if
 

      call assgn(j,gam,f,ib,f1,h,evh,qvk,qmv,qmk,qv,fistat)
        
C     save the eigenvalues  
      call esave(j,gam,f,ib,f1,qmk,evh,eused)

C     write the eigenvalues and vectors to disk
      if ((ctlint(C_INTS).gt.0).and.(fistat.eq.0)) then
        call wrvec(zr,zi,evh,j,gam,f,ib,f1,qvk,qmk)
      end if
            
C     calculate the deviation dedp
C      if (usert.eq.0) then
      complex=.true.
      if (complex) then
        do id=1, size(S_H)
          if ((eused(id).ne.0).or.(ctlint(C_DFRQ).ne.0)) then
            do ie=1, size(S_H)
              e(ie)=zr(ie,id)
              e2(ie)=zi(ie,id)
            end do
            call hcaldev_old(e,e2,j,gam,f1,ib,ifit,npar
     $           ,qvk,ruse,a,dedp,evalv,ovv,rotm,tori)
            call devsave(j,gam,f1,ib,qmk(id,Q_T)
     $           ,ifit,dfit,dedp,palc,pali)
          end if
        end do
c      else
c        do id=1, size(S_H)
c          if ((eused(id).ne.0).or.(ctlint(C_DFRQ).ne.0)) then
c            do ie=1, size(S_H)
c              e(ie)=zr(ie,id)
c            end do
c            call caldev(e,j,gam,f,ifit,npar
c     $           ,qvk,ruse,a,dedp,evalv,ovv,rotm,tori)
c            call devsave(j,gam,f,qmk(id,Q_T),ifit,dfit,dedp,palc)
c          end if
c        end do
      end if
      return
      end 
C----------------------------------------------------------------------
      subroutine calvjk_d(j,gam,f,ib,f1,evalv,ovv,rotm,rott,tori
     $     ,atot,qmv,ifittot,dfit,palc,pali,npar,fistat,hs,evhs)
C     calculation of the eigenvalues of one matrix with specified j,f,gam
C     the evalues are put in the field of dnv(1..ndata,Q_ENG,Q_UP/LO)
C     the deviations DE/DPi in dnv(1..ndata,DQ_ENG,Q_UP/LO(i))
C     fistat = 0 for regular calculation of Eigenvalues
C     fistat > 0  Eigenvalues for differential quotient

      implicit none
      include 'iam.fi'
      integer j, gam, f, ib, npar, fistat, is, f1
      real*8  h(DIMTOT,DIMTOT),evh(DIMTOT)
      real*8  hs(2,DIMTOT,DIMTOT),evhs(2,DIMTOT)
C      real*8  h1out(DIMTOT,DIMTOT),evh1out(DIMTOT)! outputs will be written to these
C      real*8  h2out(DIMTOT,DIMTOT),evh2out(DIMTOT)! For each of the two couble states.
      real*8  h_2(2*DIMTOT,2*DIMTOT), evh_2(2*DIMTOT)
C      real*8  hresort(2*DIMTOT,2*DIMTOT) !Turns out the diagonalization likely works only if I blocksort first.
C      real*8  evhresort(2*DIMTOT)
      integer sorti(2*DIMTOT)
      integer counti,countis
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      real*8  atot(DIMPAR,DIMVB)
      real*8  al(DIMPAR), au(DIMPAR)
      real*8  palc(DIMFIT,-1:DIMPLC)
      real*8  normi1, normi2 !renormalization in spearted hamiltonian matrices
      integer pali(DIMFIT, 0:DIMPLC,2)
      integer qmv(DIMV),ifittot(DIMPAR,DIMVB),dfit(DIMFIT)
      integer qmvs(2,DIMV)
C     quantum numbers
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP), qmk(DIMTOT,DIMQLP)
      integer qmks(2,DIMTOT,DIMQLP)
      integer qvks(2,DIMTOT,Q_K:Q_V+DIMTOP) 
      integer qv(DIMTOT)
      integer qvs(2,DIMTOT)
C     work
      real*8  e(DIMTOT),e2(DIMTOT),tau(2,DIMTOT)
      real*8  e_2(2*DIMTOT),e2_2(2*DIMTOT),tau_2(2,2*DIMTOT)
      real*8  zr(DIMTOT,DIMTOT),zi(DIMTOT,DIMTOT)
      real*8  zr_2(DIMTOT*2,DIMTOT*2),zi_2(DIMTOT*2,DIMTOT*2)
      real*8  zrs(2,DIMTOT,DIMTOT),zis(2,DIMTOT,DIMTOT)
      real*8  dedp(DIMPAR)
      real*8  Gx,Gy,Gz,Fxy,Fxz,Fyz,Chixy,Chiyz,Chixz
      real*8  check1,check2               
      real*8  fjn
      integer id,ie,i,iv,ik,itop,ivr,ivc,ir,ic,it1,it2,ikr,ikc
      integer ibl, ibu, ibselect
      integer eused(DIMTOT), usert,ierr
      integer ruse(DIMVV,DIMVV,DIMTOP)
      integer counter1
      integer counter2
      real*8 dj,djj1,e1,djjc,di,dii1,df,dff1,dg
      character*30 fmtstr
      character*30 fmtstr2
      character*4 fnpre
      character*6 fnpost
      real*8 tt
      logical masave
      logical complex
      integer myand
      external myand

      integer  mclock,t1,t2
      external mclock
      masave=.false.
      h=0.0
      h_2=0.0
C      h1out=0.0
C      h2out=0.0
C      hresort=0.0
C      evhresort=0.0
      sorti=0
      counti= 0
      countis=0
C      evh1out=0.0
C      evh2out=0.0
      counter1=0
      counter2=0

      dj=dble(j)
      djj1=dj*(dj+1.0)
      e1=0.0
C     djjc is used for spin rotation coupling to prevent a division by zero
C     for J=0
      djjc=1.0
      if ((ctlint(C_SPIN).ne.0).and.(j.gt.0).and.(f1.ge.0)) then
        di=dble(ctlint(C_SPIN))/2.0d0
        dii1=di*(di+1.0)
        df=dble(f1)/2.0d0!using f1 here
        dff1=df*(df+1.0)
        djjc=djj1
        dg=dff1-dii1-djj1
        if (ctlint(C_SPIN).gt.1) then
          e1= (0.75*dg*(dg+1.0)-dii1*djj1)
     $         /(2.0*di*(2.0*di-1.0)*djj1*(2.0*dj-1.0)*(2.0*dj+3.0))
        else
          e1=0.0
        end if
      end if
      
      if (ib.le.2) then
        ibl=1
        ibu=2
        if (ib.eq.1) then
          ibselect=1
        else
          ibselect=2
        end if
      else
        ibl=3
        ibu=4  
        if (ib.eq.3) then
         ibselect=1
        else
         ibselect=2
        endif
      end if     
      al=atot(:,ibl)
      au=atot(:,ibu)  

      if (ctlint(C_EVAL).gt.3)   masave=.true.

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1))
     $     write(*,'(A,5I3)')
     $     'starting with J,S,B,F Fit_stat=',J,gam,ibl,f1,fistat
      fnpre='xiam'
      if (gam.eq.0) ctlint(C_NTOP)=0
      complex=.true. ! this marks the hamiltonian as containing complex matrix elements in several subroutines. Not sure if putting it false for certain parametersets really speeds things up, I will keep it as true for now

      usert = 1


      size(S_K)=2*j+1
C     initialize the quantum no.s qvk
      i=0
      do iv=1, size(S_VV)
        do ik=1, size(S_K)
          i=i+1
          qvk(i,Q_K) =ik-j-1
          qv(i)=iv
          do itop=1, ctlint(C_NTOP)
            qvk(i,Q_V+itop)=qvv(iv,itop,ibl)
          end do
        end do
      end do
      size(S_H)=i

      do itop=1, ctlint(C_NTOP)
        do ivr=1, size(S_VV)
          do ivc=1, size(S_VV)
            ruse(ivr,ivc,itop)=1
          end do
        end do
      end do

      do it1=1, ctlint(C_NTOP)
        do it2=1, ctlint(C_NTOP)
          if (it2.ne.it1) then
            do ivr=1, size(S_VV)
              do ivc=1, size(S_VV)
                if (qvv(ivr,it2,ibl).ne.qvv(ivc,it2,ibl)) then    
                  ruse(ivr,ivc,it1)=0
                end if
              end do
            end do
          end if
        end do
      end do

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1)) then
        write(fmtstr,'(A,I1,A,I2,A)') '(',ctlint(C_NTOP),'I2,A,'
     $       ,size(S_VV),'I3)'
        write(*,*) fmtstr
        do itop=1, ctlint(C_NTOP)
          write(*,'(A,I3)') ' ruse array for top', itop 
          do ivr=1, size(S_VV)
            write(*,fmtstr)(qvv(ivr,it2,ibl),it2=1,ctlint(C_NTOP))
     $           ,' | ',(ruse(ivr,ivc,itop),ivc=1, size(S_VV))
          end do
        end do
      end if 

      if (masave) then
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'t'//fnpost,status='unknown')
          write(55,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              write(55,*)
     $             ((rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,itop)
     $             *tori(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc
     $             ,gamma(gam,itop),itop)
     $             ,ikc=1,size(S_K))
     $             ,ivc=1, size(S_V+itop))
            end do
          end do
          close(55)
        end do
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(56,file=fnpre//'p'//fnpost,status='unknown')
          write(56,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              do ivc=1, size(S_V+itop)
                do ikc=1, size(S_K)
                  if (qvk(ikr,Q_K).eq.qvk(ikc,Q_K)) then
                    write(56,'(D22.14,$)')
     $                   ovv(ivr,ivc,PM_PI,gamma(gam,itop)
     $                   ,qvk(ikr,Q_K),itop)
                  else
                    write(56,'(D22.14,$)') 0.0
                  end if
                end do
              end do
              write(56,*)
            end do
          end do
          close(56)
        end do
        
        do it1=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         it1,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'r'//fnpost,status='unknown')
          write(55,*) size(S_H),0,0
          do ir=1, size(S_H)
            do ic=1, size(S_H)
              tt=1.0d0
              do it2=1,ctlint(C_NTOP)
                tt=tt*tori(qvk(ir,Q_K),    qvk(ic,Q_K),
     $                   qvk(ir,Q_V+it2),qvk(ic,Q_V+it2),
     $                   gamma(gam,it2),it2)
              end do
              write(55,'(D22.14,$)')
     $             rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,it1)
     $             *tt
            end do
            write(55,*)
          end do
          close(55)
        end do
      end if

      if (gam.ne.0) then 
        call bld2vjk(j,gam,f1
     $         ,qvk,ruse,h,al,evalv,ovv,rotm,rott,tori,complex)

C        end if

c      write(*,*) 'bld2vjk',mclock()-t1

      else
        if (size(S_VV).gt.1) stop ' size vv > 1 for rigid rotor!'
      end if
C      t1=mclock()
      call addrig_old(j,gam,f1
     $     ,qvk,ruse,h,al,evalv,ovv,rotm,rott,tori,complex)
C---
C      t1=mclock()
       do i=1, DIMTOT
         do ie=1, DIMTOT
           h_2(i,ie)=h(i,ie) !copying h to h2
         end do
       end do 
      qvks(ibl,:,:)=qvk
      qmks(ibl,:,:)=qmk
      qmvs(ibl,:)=qmv
      qvs(ibl,:)=qv


C     Starting a copy with ibu
C     Starting a copy with ibu
C     Starting a copy with ibu
C     Starting a copy with ibu
C     Starting a copy with ibu
C     Starting a copy with ibu
      i=0
      do iv=1, size(S_VV)
        do ik=1, size(S_K)
          i=i+1
          qvk(i,Q_K) =ik-j-1
          qv(i)=iv
          do itop=1, ctlint(C_NTOP)
            qvk(i,Q_V+itop)=qvv(iv,itop,ibu)
          end do
        end do
      end do
      size(S_H)=i
      do ivr=1, size(S_H)
        do ivc=1, size(S_H)
          h(ivr,ivc)=0.0
        end do
      end do

      do itop=1, ctlint(C_NTOP)
        do ivr=1, size(S_VV)
          do ivc=1, size(S_VV)
            ruse(ivr,ivc,itop)=1
          end do
        end do
      end do

      do it1=1, ctlint(C_NTOP)
        do it2=1, ctlint(C_NTOP)
          if (it2.ne.it1) then
            do ivr=1, size(S_VV)
              do ivc=1, size(S_VV)
                if (qvv(ivr,it2,ibu).ne.qvv(ivc,it2,ibu)) then    
                  ruse(ivr,ivc,it1)=0
                end if
              end do
            end do
          end if
        end do
      end do

      if ((myand(ctlint(C_PRI),AP_ST).ne.0).and.(xde.ge.1)) then
        write(fmtstr,'(A,I1,A,I2,A)') '(',ctlint(C_NTOP),'I2,A,'
     $       ,size(S_VV),'I3)'
        write(*,*) fmtstr
        do itop=1, ctlint(C_NTOP)
          write(*,'(A,I3)') ' ruse array for top', itop 
          do ivr=1, size(S_VV)
            write(*,fmtstr)(qvv(ivr,it2,ibu),it2=1,ctlint(C_NTOP))
     $           ,' | ',(ruse(ivr,ivc,itop),ivc=1, size(S_VV))
          end do
        end do
      end if 

      if (masave) then
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'t'//fnpost,status='unknown')
          write(55,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              write(55,*)
     $             ((rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,itop)
     $             *tori(qvk(ikr,Q_K),qvk(ikc,Q_K),ivr,ivc
     $             ,gamma(gam,itop),itop)
     $             ,ikc=1,size(S_K))
     $             ,ivc=1, size(S_V+itop))
            end do
          end do
          close(55)
        end do
        do itop=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         itop,gam,'.j',int(j/10),mod(j,10)
          open(56,file=fnpre//'p'//fnpost,status='unknown')
          write(56,*) size(S_V+itop)*size(S_K),0,0
          do ivr=1, size(S_V+itop)
            do ikr=1, size(S_K)
              do ivc=1, size(S_V+itop)
                do ikc=1, size(S_K)
                  if (qvk(ikr,Q_K).eq.qvk(ikc,Q_K)) then
                    write(56,'(D22.14,$)')
     $                   ovv(ivr,ivc,PM_PI,gamma(gam,itop)
     $                   ,qvk(ikr,Q_K),itop)
                  else
                    write(56,'(D22.14,$)') 0.0
                  end if
                end do
              end do
              write(56,*)
            end do
          end do
          close(56)
        end do
        
        do it1=1,ctlint(C_NTOP)
          write(fnpost,'(I1,I1,A,I1,I1)')
     $         it1,gam,'.j',int(j/10),mod(j,10)
          open(55,file=fnpre//'r'//fnpost,status='unknown')
          write(55,*) size(S_H),0,0
          do ir=1, size(S_H)
            do ic=1, size(S_H)
              tt=1.0d0
              do it2=1,ctlint(C_NTOP)
                tt=tt*tori(qvk(ir,Q_K),    qvk(ic,Q_K),
     $                   qvk(ir,Q_V+it2),qvk(ic,Q_V+it2),
     $                   gamma(gam,it2),it2)
              end do
              write(55,'(D22.14,$)')
     $             rotm(qvk(ikr,Q_K),qvk(ikc,Q_K),1,it1)
     $             *tt
            end do
            write(55,*)
          end do
          close(55)
        end do
      end if

      if (gam.ne.0) then 

        call bld2vjk(j,gam,f1
     $         ,qvk,ruse,h,au,evalv,ovv,rotm,rott,tori,complex)

      else
        if (size(S_VV).gt.1) stop ' size vv > 1 for rigid rotor!'
      end if
C      t1=mclock()
      call addrig_old(j,gam,f1
     $     ,qvk,ruse,h,au,evalv,ovv,rotm,rott,tori,complex)
        

C---
C      t1=mclock()
        do i=1, size(S_H)
          do ie=1, size(S_H)
            h_2(i+size(S_H),ie+size(S_H))=h(i,ie) !copying 2nd h to h2
          end do
        end do
        
       qvks(ibu,:,:)=qvk
       qmks(ibu,:,:)=qmk
       qmvs(ibu,:)=qmv
        qvs(ibu,:)=qv
      
C      End of copy
C      End of copy        

        
C      The following implementation of Gx,y,z and Fxy, Fyz and Fxz, differes from that in SPFIT/SPCAT       
C      While if all Fs and all Gs are fit on their own, the Data simulated using SPFIT/SPCAT is reproduced
C      Mixing the Fs and Gs will lead to disagreement.       
C      This is likely due to inconsistensies in phase convention in my implementation here.
C      For now, I discourage mixing G and F parameters until this disagreement is fixed.
C      
C             --- Sven 25-07-2024
       Gz    = 0.0
       Gy    = 0.0
       Gx    = 0.0
       Fxy   = 0.0
       Fxz   = 0.0
       Fyz   = 0.0
       Chixy = 0.0
       Chixz = 0.0
       Chiyz = 0.0
       if (ib.le.2) then
         Gz =al(P_GZ12)
         Gy =al(P_GY12)
         Gx =al(P_GX12)
         Fxy=al(P_FXY1)
         Fxz=al(P_FXZ1)
         Fyz=al(P_FYZ1)
         Chixy=al(P_WQXY1)*(-1.0) !Sign change to match relative signs in spfit output !these are offdiagonal nqcc matrix elements but used offdiagonal in v. Matrix elements offdiagonal in J neglected.
         Chiyz=al(P_WQYZ1)*(-1.0) !
         Chixz=al(P_WQXZ1)*(-1.0) !
       else
         Gz =al(P_GZ34)
         Gy =al(P_GY34)
         Gx =al(P_GX34)
         Fxy=al(P_FXY3)  
         Fxz=al(P_FXZ3)!
         Fyz=al(P_FYZ3)!      
         Chixy=al(P_WQXY3)*(-1.0)
         Chiyz=al(P_WQYZ3)*(-1.0)
         Chixz=al(P_WQXZ3)*(-1.0)
       end if
        
C      ADDING OFF DIAGONAL ELEMENTS FOR GX, GY, GZ
       if (Gz .ne. 0.0) then
         do ik=1, 2*j+1
           h_2(ik,size(S_H)+ik)=Gz*1.0*(ik-1-j)! on complex off diag Gordy_Cook eq 7.135 p290 
         end do !  
       end if 
       if (Gy .ne. 0.0) then
         do ik=1, 2*j
             h_2(ik,size(S_H)+ik+1)=h_2(ik,size(S_H)+ik+1)+0.5*Gy
     $                *sqrt(1.0*(j-(ik-1-j))*(j+(ik-1-j)+1)) ! on complex off diag
             h_2(ik+1,size(S_H)+ik)=h_2(ik+1,size(S_H)+ik)+0.5*Gy
     $                *-1.0*sqrt(1.0*(j-(ik-1-j))*(j+(ik-1-j)+1)) ! on complex off diag
         end do
       end if 
       if (Gx .ne. 0.0) then
         do ik=2, 2*j+1
             h_2(size(S_H)+ik,ik-1)=h_2(size(S_H)+ik,ik-1)+0.5*Gx
     $                *sqrt(1.0*(j+(ik-1-j))*(j-(ik-1-j)+1)) ! on real off diag
              h_2(size(S_H)+ik-1,ik)=h_2(size(S_H)+ik-1,ik)+0.5*Gx
     $                *sqrt(1.0*(j+(ik-1-j))*(j-(ik-1-j)+1)) ! on real off diag
         end do 
       end if
       
       !FXY is imaginary!Based on Evaluation and optimal computation of angular momentum matrix elements: An information theory approach
       if ((Fxy .ne. 0.0).or.(Chixy .ne. 0.0)) then !and based on doi.org/10.1063/1.1677430
         do ik=1, 2*j-1
           if(j.eq.1)then
             fjn=(0.5*j*(j+1))**2
           else
             fjn=0.25*(j*(j+1)-(ik-j)*((ik-j)+1))
     $          *(j*(j+1)-(ik-j)*((ik-j)-1))

           end if 
             h_2(ik,size(S_H)+ik+2)=h_2(ik,size(S_H)+ik+2) 
     $        +sqrt(fjn)*(Fxy-2*Chixy*e1)
             h_2(ik+2,size(S_H)+ik)=h_2(ik+2,size(S_H)+ik)
     $         -1.0*sqrt(fjn)*(Fxy-2*Chixy*e1)
         end do
       end if        
       
       ! FYZ is put on the imaginary
       if ((Fyz .ne. 0.0).or.(Chiyz .ne. 0.0)) then !Based on Evaluation and optimal computation of angular momentum matrix elements: An information theory approach
         do ik=1, 2*j
             h_2(ik,size(S_H)+ik+1)=
     $  h_2(ik,size(S_H)+ik+1)+0.5*(2*(ik-1-j)+1)
     $  *(j**2+j-(ik-1-j)**2-(ik-1-j))**0.5*(Fyz-2.0*Chiyz*e1)
             h_2(ik+1,size(S_H)+ik)=
     $  h_2(ik+1,size(S_H)+ik)-0.5*(2*(ik-1-j)+1)
     $  *(j**2+j-(ik-1-j)**2-(ik-1-j))**0.5*(Fyz-2.0*Chiyz*e1)
         end do
       end if        
       !FXZ is put on the real
       if ((Fxz .ne. 0.0).or.(Chixz .ne. 0.0)) then !Based on Evaluation and optimal computation of angular momentum matrix elements: An information theory approach
         do ik=1, 2*j
            h_2(size(S_H)+ik+1,ik)=h_2(size(S_H)+ik+1,ik)+(0.5* 
     $         (2*(ik-1-j)+1)*(j**2+j-(ik-1-j)**2-(ik-1-j))**0.5)
     $          *(Fxz-2.0*Chixz*e1)
             h_2(size(S_H)+ik,ik+1)=h_2(size(S_H)+ik,ik+1)+(0.5*
     $          (2*(ik-1-j)+1)*(j**2+j-(ik-1-j)**2-(ik-1-j))**0.5)
     $          *(Fxz-2.0*Chixz*e1)
         end do
       end if        
       
        do ir=1, 2*size(S_H)
          do ic=1, 2*size(S_H)
            zr_2(ir,ic)=0.0
          end do
          zr_2(ir,ir)=1.0
        end do
        ierr=0
              
        
C        The problem is that there can be 0 energy levels e.g. J=0. 
C        however the diagonalization routine does not accept zero entries. 
C        This means I have to check, if for both rotational states the entries are zero.
C        If this is not the case, I will add a small offset to the one that is non zero
        do i = 1, 2*size(S_H) !
         check1=abs(h_2(i,i))    
         if (check1.le.1.0e-14) then
             h_2(i,i)=5.0e-14 
         end if
        end do     
             
         
C      hresort(1:size(S_H),1:size(S_H))= !resorting should not be needed anymore since DIMTOT was replaced with size(S_H) in construction
C     $ h_2(1:size(S_H),1:size(S_H))  ! Diagonal Block lower state
C      hresort(size(S_H)+1:2*size(S_H),size(S_H)+1:2*size(S_H))=
C     $ h_2(DIMTOT+1:DIMTOT+size(S_H),DIMTOT+1:DIMTOT+size(S_H))  ! Diagonal Block upper state
C      hresort(1:size(S_H),size(S_H)+1:2*size(S_H))=
C     $ h_2(1:size(S_H),DIMTOT+1:DIMTOT+size(S_H))     ! Offdiagonal 1
C      hresort(size(S_H)+1:2*size(S_H),1:size(S_H))=  
C     $ h_2(DIMTOT+1:DIMTOT+size(S_H),1:size(S_H))     ! Offdiagonal 2

C        h_2=hresort
        
        call htrid3 (2*DIMTOT,2*size(s_h),h_2,evh_2,e_2,e2_2,tau_2)
        call tql2 (2*DIMTOT,2*size(s_h),evh_2,e_2,zr_2,ierr)
        
        
        if (ierr.ne.0) then
          write (*,'(a,i5)') 'Error in tql2 ',ierr
          stop
        endif
        call htrib3 (2*DIMTOT,2*size(S_H),h_2,tau_2
     &           ,2*size(S_H),zr_2,zi_2)

C     sort eigenvalues in **ascending** order
        call heigsrt(evh_2,zr_2,zi_2,2*size(S_H),2*DIMTOT)

        do i=1, 2*size(S_H)
          do ie=1, 2*size(S_H)
            h_2(i,ie)=sign(dsqrt(zr_2(i,ie)**2+zi_2(i,ie)**2)
     $           ,zr_2(i,ie)+zi_2(i,ie))              
          end do
        end do
        
        counter1=0
        counter2=0
        do i = 1, 2*size(S_H)   !adding and renormalizing the vector components....
          normi1=sum((h_2(1:size(S_H),i)**2)) 
          normi2=sum((h_2(size(S_H)+1:2*size(S_H),i)**2))
         if (normi1.ge.0.5) then !e.g. for empty matrix entries where all normis are zero
          counter1=counter1+1
          evhs(1,counter1)=evh_2(i)
          hs(1,1:size(S_H),counter1)=SIGN(h_2(1:size(S_H),i)**2+
     $         h_2(size(S_H)+1:2*size(S_H),i)**2,h_2(1:size(S_H),i))  ! This is different for normalization - I use the K info from both states.
         else
          counter2=counter2+1
          evhs(2,counter2)=evh_2(i)
          hs(2,1:size(S_H),counter2)=SIGN(h_2(1:size(S_H),i)**2+
     $         h_2(size(S_H)+1:2*size(S_H),i)**2,
     $         h_2(size(S_H)+1:2*size(S_H),i))  ! This is different for normalization - I use the K info from both states.
         end if
        end do
        
C       do i=1, 2*size(S_H)
C          do ie=1, 2*size(S_H)
C            zi_2(i,ie)=0.0
C            zr_2(i,ie)=h_2(i,ie)
C          end do
C       end do        
       do i=1, size(S_H)
          do ie=1, size(S_H)
            zis(1,i,ie)=0.0
            zrs(1,i,ie)=hs(1,i,ie)
          end do
       end do        
       do i=1, size(S_H)
          do ie=1, size(S_H)
            zis(2,i,ie)=0.0
            zrs(2,i,ie)=hs(2,i,ie)
          end do
       end do        
      call assgn(j,gam,f,ib,f1,hs(ibselect,:,:),evhs(ibselect,:)
     $           ,qvks(ibselect,:,:),qmvs(ibselect,:),
     $           qmks(ibselect,:,:),qvs(ibselect,:),fistat)
        
      call esave(j,gam,f,ib,f1,qmks(ibselect,:,:)
     $           ,evhs(ibselect,:),eused)
C
C     write the eigenvalues and vectors to disk
      if ((ctlint(C_INTS).gt.0).and.(fistat.eq.0)) then
       call wrvec(zrs(ibselect,:,:),zis(ibselect,:,:),evhs(ibselect,:)
     $           ,j,gam,f,ib,f1,qvks(ibselect,:,:),qmks(ibselect,:,:))
      end if
            
C     calculate the deviation dedp ! removed analytical gradients for now.
C     complex=.true.
C     if (complex) then
C       do id=1, size(S_H)
C         if ((eused(id).ne.0).or.(ctlint(C_DFRQ).ne.0)) then
C           do ie=1, size(S_H)
C             e(ie)=zrs(ib,ie,id)
C             e2(ie)=zis(ib,ie,id)
C           end do
C           call hcaldev(e,e2,j,gam,f,ib,ifittot(:,ib),npar
C    $           ,qvks(ibselect,:,:),ruse,atot(:,ib),dedp
C    $           ,evalv,ovv,rotm,tori)
C           call devsave(j,gam,f,ib,qmks(ib,id,Q_T)
C    $           ,ifittot(:,ib),dfit,dedp,palc,pali)
C         end if
C       end do
C      end if
C      qmk= qmks(ibselect,:,:)
C      qvk=qvks(ibselect,:,:)
C      zr=zrs(ibselect,:,:)
C      zi=zis(ibselect,:,:)
C      evh=evhs(ibselect,:)
C      h=hs(ibselect,:,:)
       return
       end 

C----------------------------------------------------------------------
      subroutine addrig(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,complex)
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      logical complex
C     work
      real*8 vr(DIMTOT), vi(DIMTOT), vor(DIMTOT), voi(DIMTOT)
      integer i,iv

      do i =1, size(S_H)
        do iv=1, size(S_H)
          vr(iv)=0.0
          vor(iv)=0.0
          vi(iv)=0.0
          voi(iv)=0.0
        end do
        vr(i )=1.0d0
        call hmulthrr(j,gam,f
     $       ,qvk,ruse,a,vr,vi,vor,voi,evalv,ovv,rotm,tori,0,0)
        do iv=1, i
          h(i,iv)=h(i,iv)+vor(iv)
        end do
        do iv=i+1, size(S_H)
          h(i,iv)=h(i,iv)+voi(iv)
        end do
      end do
      return
      end
C----------------------------------------------------------------------
      subroutine addrig_old(j,gam,f,qvk,ruse
     $     ,h,a,evalv,ovv,rotm,rott,tori,complex)
      implicit none
      include 'iam.fi'
      integer j,gam,f
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  h(DIMTOT,DIMTOT),a(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      logical complex
C     work
      real*8 vr(DIMTOT), vi(DIMTOT), vor(DIMTOT), voi(DIMTOT)
      integer i,iv

      do i =1, size(S_H)
        do iv=1, size(S_H)
          vr(iv)=0.0
          vor(iv)=0.0
          vi(iv)=0.0
          voi(iv)=0.0
        end do
        vr(i )=1.0d0
        call hmulthrr_old(j,gam,f
     $       ,qvk,ruse,a,vr,vi,vor,voi,evalv,ovv,rotm,tori,0,0)
        do iv=1, i
          h(i,iv)=h(i,iv)+vor(iv)
        end do
        do iv=i+1, size(S_H)
          h(i,iv)=h(i,iv)+voi(iv)
        end do
      end do
      return
      end

C----------------------------------------------------------------------
      subroutine hcaldev(vr,vi,j,gam,f,ib,ifit,npar
     $     ,qvk,ruse,a,dedp,evalv,ovv,rotm,tori)
C     (hermitian) calculation of the < J tau | Op | J tau' >  
C     Hellmann-Feynman-Theorem?
      implicit none
      include 'iam.fi'
      real*8  vi(DIMTOT),vr(DIMTOT)
      integer j,gam,f,npar,ib
      integer qvk(DIMTOT,Q_K:Q_V1+DIMTOP-1), ifit(DIMPAR)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  a(DIMPAR),dedp(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      real*8  da(DIMPAR)
      real*8  vor(DIMTOT),voi(DIMTOT)
      integer i,ifs

      do ifs=1, DIMPRR
        dedp(ifs)=0.0
      end do
      do ifs=1, DIMPRR
        if (ifit(ifs).ne.0) then
          do i=1,DIMPRR
            da(i)=0.0
          end do
          da(ifs)=1.0
          do i=1, size(S_H)
            vor(i)=0.0
            voi(i)=0.0
          end do
          call hmulthrr(j,gam,f
     $         ,qvk,ruse,da,vr,vi,vor,voi,evalv,ovv,rotm,tori,ifs,0)
          do i=1, size(S_H)
            dedp(ifs)=dedp(ifs)+vor(i)*vr(i)+voi(i)*vi(i)
          end do
        end if
      end do
      return
      end
C----------------------------------------------------------------------
      subroutine hcaldev_old(vr,vi,j,gam,f,ib,ifit,npar
     $     ,qvk,ruse,a,dedp,evalv,ovv,rotm,tori)
C     (hermitian) calculation of the < J tau | Op | J tau' >  
C     Hellmann-Feynman-Theorem? 
      implicit none
      include 'iam.fi'
      real*8  vi(DIMTOT),vr(DIMTOT)
      integer j,gam,f,npar,ib
      integer qvk(DIMTOT,Q_K:Q_V1+DIMTOP-1), ifit(DIMPAR)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  a(DIMPAR),dedp(DIMPAR)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
C     work
      real*8  da(DIMPAR)
      real*8  vor(DIMTOT),voi(DIMTOT)
      integer i,ifs

      do ifs=1, DIMPRR
        dedp(ifs)=0.0
      end do
      do ifs=1, DIMPRR
        if (ifit(ifs).ne.0) then
          do i=1,DIMPRR
            da(i)=0.0
          end do
          da(ifs)=1.0
          do i=1, size(S_H)
            vor(i)=0.0
            voi(i)=0.0
          end do
          call hmulthrr_old(j,gam,f
     $         ,qvk,ruse,da,vr,vi,vor,voi,evalv,ovv,rotm,tori,ifs,0)
          do i=1, size(S_H)
            dedp(ifs)=dedp(ifs)+vor(i)*vr(i)+voi(i)*vi(i)
          end do
        end if
      end do
      return
      end
      
C----------------------------------------------------------------------
      subroutine hmulthrr(j,gam,f
     $     ,qvk,ruse,a,vr,vi,vor,voi,evalv,ovv,rotm,tori,ifs,it)
C     multiply the complex vector vr,vi by the rigid part of the
C     Hamilton matrix, yielding vor,voi
      implicit none
      include 'iam.fi'
      integer j,gam,f,ifs,it
      integer qvk(DIMTOT,Q_K:Q_V1+DIMTOP-1)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  a(DIMPAR)
      real*8  vr(DIMTOT),vor(DIMTOT),vi(DIMTOT),voi(DIMTOT)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)

C     work
      real*8  dk,dj,djj1,dff,t,t1,t2,djjc
      real*8  e1,df,dff1,di,dii1,dg
      real*8  adelk,adelj,ah2,ah3,ahk,ar6,ahjk,ahj
      real*8  alj,aljk,alkj,alk !Herbers2026
      real*8  al2,al3,al4 !Herbers2026
      real*8  DXTERM, DYTERM, DZTERM
      integer ik,ir,ic,ivr,ivc,itop
      integer off,voff
      integer myand
      external myand
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in HMULTHRR'

      dj=dble(j)
      djj1=dj*(dj+1.0)
      e1=0.0
C     djjc is used for spin rotation coupling to prevent a division by zero
C     for J=0
      djjc=1.0
      if ((ctlint(C_SPIN).ne.0).and.(j.gt.0).and.(f.ge.0)) then
        di=dble(ctlint(C_SPIN))/2.0d0
        dii1=di*(di+1.0)
        df=dble(f)/2.0d0
        dff1=df*(df+1.0)
        djjc=djj1
        dg=dff1-dii1-djj1
        if (ctlint(C_SPIN).gt.1) then
          e1= (0.75*dg*(dg+1.0)-dii1*djj1)
     $         /(2.0*di*(2.0*di-1.0)*djj1*(2.0*dj-1.0)*(2.0*dj+3.0))
        else
          e1=0.0
        end if
      end if

C     the centrifugal distortion parameters
C     Watson A 
      if (ctlint(C_RED).eq.0) then
        adelj=a(P_DJD)
        adelk=a(P_DKD)
        ar6  =0.0d0
        ahj  =a(P_HJD)
        ahjk =a(P_HJKD)
        ah2  =0.0d0
        ahk  =a(P_HKD)
        ah3  =0.0d0
        alj = a(P_LJD)
        aljk = a(P_LJKD)
        alkj = a(P_LKJD)
        alk  = a(P_LKD)
      end if
C     Watson S
      if (ctlint(C_RED).eq.1) then
        adelj=-a(P_DJD)
        adelk=0.0d0
        ar6  =a(P_DKD)
        ahj  =a(P_HJD)
        ahjk =0.0d0
        ah2  =a(P_HJKD)
        ahk  =0.0d0
        ah3  =a(P_HKD)
        alj = a(P_LJD)
        al2 = a(P_LJKD)
        al3 = a(P_LKJD)
        al4 = a(P_LKD)
      end if
C     van Eijck / Typke
      if (ctlint(C_RED).eq.2) then
        adelj=a(P_DJD)
        adelk=0.0d0
        ar6  =a(P_DKD)
        ahj  =0.5d0*a(P_HJD)
        ahjk =0.0d0
        ah2  =0.25d0*a(P_HJKD)
        ahk  =0.0d0
        ah3  =0.125d0*a(P_HKD)
        alj  =0.5*a(P_LJD) ! no further implementation for octic.
      end if

      if ( (a(P_BJ) .ne.0.0).or.(a(P_BK) .ne.0.0).or.
     $     (a(P_DJ) .ne.0.0).or.(a(P_DJK).ne.0.0).or.
     $     (a(P_DK) .ne.0.0).or.(a(P_HJ ).ne.0.0).or.
     $     (a(P_HJK).ne.0.0).or.(a(P_HKJ).ne.0.0).or.
     $     (a(P_HK) .ne.0.0).or.(a(P_QZ) .ne.0.0).or.
     $     (a(P_LJ) .ne.0.0).or.(a(P_LK) .ne.0.0).or.
     $     (a(P_LJJK) .ne.0.0).or.(a(P_LJK) .ne.0.0).or.
     $     (a(P_LKKJ) .ne.0.0).or.
     $     (a(P_E)  .ne.0.0)) then !Herbers2023
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          t=     a(P_BJ) *djj1
     $         + a(P_BK) *dk*dk
     $         - a(P_DJ) *djj1**2
     $         - a(P_DJK)*djj1*dk**2
     $         - a(P_DK) *dk**4
     $         + a(P_HJ) *djj1**3
     $         + a(P_HJK)*(djj1**2)*(dk**2)
     $         + a(P_HKJ)*djj1*(dk**4)
     $         + a(P_HK) *dk**6                      !     $         + a(P_QZ) * e1 * (3.0*(dk**2)-djj1) remvoed
     $         + a(P_CP) * 0.5*dg*(1.0-dk*dk/djjc)
     $         + a(P_CZ) * 0.5*dg*dk*dk/djjc
     $         + a(P_E)  !Herbers2023
     $         + a(P_LJ) *djj1**4                 !Herbers2026
     $         + a(P_LK) *dk**8                   !Herbers2026
     $         + a(P_LJK)*(djj1**2)*(dk**4)       !Herbers2026
     $         + a(P_LJJK)*(djj1**3)*(dk**2)      !Herbers2026
     $         + a(P_LKKJ)*(djj1**1)*(dk**6)      !Herbers2026
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
      
      if (a(P_PZ) .ne.0.0) then
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          t=     a(P_PZ) *dk
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
      
      DZTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)   DZTERM=a(P_DZ1)  !herbers2024
      if (gam.eq.2)   DZTERM=a(P_DZ2)  !herbers2024
      if (gam.eq.3)   DZTERM=a(P_DZ3)  !herbers2024
      if (gam.eq.4)   DZTERM=a(P_DZ4)  !herbers2024
      if (gam.eq.5)   DZTERM=a(P_DZ5)  !herbers2024
      if (gam.eq.6)   DZTERM=a(P_DZ6)  !herbers2024
      if (gam.eq.7)   DZTERM=a(P_DZ7)  !herbers2024
      if (gam.eq.8)   DZTERM=a(P_DZ8)  !herbers2024
      if (gam.eq.9)   DZTERM=a(P_DZ9)  !herbers2024
      if (gam.eq.10)  DZTERM=a(P_DZ10) !herbers2024
      if (gam.eq.11)  DZTERM=a(P_DZ11) !herbers2024
      if (DZTERM.ne.0.0) then          !herbers2024
        do ik=1, size(S_K)             !herbers2024
          dk=dble(qvk(ik,Q_K))         !herbers2024
          t=     DZTERM *dk            !herbers2024
          do ivc=1, size(S_VV)         !herbers2024
            ic=ik+size(S_K)*(ivc-1)    !herbers2024
            vor(ic)=vor(ic)+vr(ic)*t   !herbers2024
            voi(ic)=voi(ic)+vi(ic)*t   !herbers2024
          end do                       !herbers2024
        end do                         !herbers2024
      end if                           !herbers2024
      
      
      
      if ((ctlint(C_RED).eq.2).and. 
     $     ((ar6.ne.0.0d0).or.(ah2.ne.0.0d0))) then
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          dff=  (djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         +(djj1-dk*(dk-1.0))*(djj1-(dk-1.0)*(dk-2.0))
     $         -2.0d0*(djj1-dk**2)**2
          t=  (ar6 + ah2*djj1)*dff
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
        
CC     real off diagonal k/k+1   ! removed 
C      if ((a(P_QXZ).ne.0.0)) then
C        off=1
C        do ik=1, size(S_K)-off
C          dk=dble(qvk(ik,Q_K))
C          dff=(1.0+2.0*dk)*dsqrt(djj1-dk*(dk+1.0))
C          t=    e1* a(P_QXZ) *dff  
C          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
C        end do
C      end if

C     real off diagonal k/k+1 ! in XIAM P_PX is real and P_PY gives imaginary, this is flipped convention from Gordy 2.63
      if ((a(P_PX).ne.0.0).or.(a(P_DZX).ne.0.0).or.
     $     (a(P_DZXJ).ne.0.0).or.(a(P_DZXK).ne.0.0)) then!Adding Dab, DabJ, DabK here
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=  a(P_PX)*dff
     $       +a(P_DZX)*(dk+(dk+1.0))*dff
     $       +a(P_DZXJ)*(dk+(dk+1.0))*djj1*dff
     $       +a(P_DZXK)*(dk**3+(dk+1.0)**3)*dff     
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if
      
      
      
      DXTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)  DXTERM=a(P_DX1)   !herbers2024
      if (gam.eq.2)  DXTERM=a(P_DX2)   !herbers2024
      if (gam.eq.3)  DXTERM=a(P_DX3)   !herbers2024
      if (gam.eq.4)  DXTERM=a(P_DX4)   !herbers2024
      if (gam.eq.5)  DXTERM=a(P_DX5)   !herbers2024
      if (gam.eq.6)  DXTERM=a(P_DX6)   !herbers2024
      if (gam.eq.7)  DXTERM=a(P_DX7)   !herbers2024
      if (gam.eq.8)  DXTERM=a(P_DX8)   !herbers2024
      if (gam.eq.9)  DXTERM=a(P_DX9)   !herbers2024
      if (gam.eq.10) DXTERM=a(P_DX10)   !herbers2024
      if (gam.eq.11) DXTERM=a(P_DX11)   !herbers2024
      if (DXTERM.ne.0.0) then          !herbers2024
C       real off diagonal k/k+1 
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=  DXTERM *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if      

CC     imaginaer off diagonal k/k+1  ! removed
C      if ((a(P_QYZ).ne.0.0)) then
C        off=1
C        do ik=1, size(S_K)-off
C          dk=dble(qvk(ik,Q_K))
C          dff=(1.0+2.0*dk)*dsqrt(djj1-dk*(dk+1.0))
C          t=    e1* a(P_QYZ) *dff *(-1) ! Sign change added by Sven 2024
C          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
C        end do
C      end if

C     imaginaer off diagonal k/k+1 
      if ((a(P_PY).ne.0.0)) then
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=   a(P_PY) *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if
      
      DYTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)  DYTERM=a(P_DY1)   !herbers2024
      if (gam.eq.2)  DYTERM=a(P_DY2)   !herbers2024
      if (gam.eq.3)  DYTERM=a(P_DY3)   !herbers2024
      if (gam.eq.4)  DYTERM=a(P_DY4)   !herbers2024
      if (gam.eq.5)  DYTERM=a(P_DY5)   !herbers2024
      if (gam.eq.6)  DYTERM=a(P_DY6)   !herbers2024
      if (gam.eq.7)  DYTERM=a(P_DY7)   !herbers2024
      if (gam.eq.8)  DYTERM=a(P_DY8)   !herbers2024
      if (gam.eq.9)  DYTERM=a(P_DY9)   !herbers2024
      if (gam.eq.10) DYTERM=a(P_DY10)   !herbers2024
      if (gam.eq.11) DYTERM=a(P_DY11)   !herbers2024
      if (DYTERM.ne.0.0) then          !herbers2024
C       imaginaer off diagonal k/k+1 
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=   DYTERM *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if
      
C     real off diagonal k/k+2 
      if ((a(P_BD).ne.0.0).or.(adelj.ne.0.0).or.
     $    (adelk.ne.0.0).or.(a(P_QD).ne.0.0).or.
     $    (ahj.ne.0.0).or.(a(P_CD).ne.0.0).or.
     $    (ahjk.ne.0.0).or.(ahk.ne.0.0).or.(alj.ne.0.0).or.
     $    (aljk.ne.0.0).or.(alkj.ne.0.0).or.(alk.ne.0.0)
     $     ) then
        off=2
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
          t=     a(P_BD) *dff
     $         - adelj   *2.0d0*dff*djj1
     $         - adelk   *dff*((dk+2.0d0)**2+dk**2)
     $         + ahj*2.0d0*dff*djj1**2
     $         + ahjk*dff*((dk+2.0d0)**2+dk**2)*djj1
     $         + ahk*dff*((dk+2.0d0)**4+dk**4)   !     $         + a(P_QD) *dff*e1 removed
     $         + a(P_CD) *0.5*dg*dff/djjc
     $         + alj*2.0d0*dff*djj1**3                        !Herbers2026
     $         + aljk*dff*((dk+2.0d0)**2+dk**2)*djj1**2       !Herbers2026
     $         + alkj*dff*((dk+2.0d0)**4+dk**4)*djj1          !Herbers2026
     $         + alk*dff*((dk+2.0d0)**6+dk**6)                !Herbers2026
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

      if ((ctlint(C_RED).eq.2).and.(ah3.ne.0.0d0)) then
        off=2
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
     $         *((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         +(djj1-dk*(dk-1.0))*(djj1-(dk-1.0)*(dk-2.0))
     $         +(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0)))
          t=  ah3*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

CC     imaginaer off diagonal k/k+2  ! removed
C      if ((a(P_QXY).ne.0.0)) then
C        off=2
C        do ik=1, size(S_K)-off
C          dk=dble(qvk(ik,Q_K))
C          dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
C          t=    
C     $         + a(P_QXY)*dff*e1*2.0*(-1.0) ! Sign change Sven 2024
C          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
C        end do
C      end if
 
C     Watson S off diagonal k/k+4 (evtl. change dff)
      if ((ar6.ne.0.0d0).or.(ah2.ne.0.0d0).or.(al2.ne.0.0d0)) then
        off=4
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))) 
          t=   
     $          ar6*dff + ah2*djj1*dff + al2*djj1**2*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

C     Watson S off diagonal k/k+6 
      if ((ah3.ne.0.0d0).or.(al3.ne.0.0)) then
        off=6
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))
     $         *(djj1-(dk+4.0)*(dk+5.0))*(djj1-(dk+5.0)*(dk+6.0))) 
          t=   
     $          ah3*dff+al3*djj1*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if
      
C     Watson S off diagonal k/k+8                                     !Herbers2026
      if (al4.ne.0.0) then                                            !Herbers2026
        off=8                                                         !Herbers2026
        do ik=1, size(S_K)-off                                        !Herbers2026
          dk=dble(qvk(ik,Q_K))                                        !Herbers2026
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))       !Herbers2026
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))     !Herbers2026
     $         *(djj1-(dk+4.0)*(dk+5.0))*(djj1-(dk+5.0)*(dk+6.0))     !Herbers2026
     $         *(djj1-(dk+6.0)*(dk+7.0))*(djj1-(dk+7.0)*(dk+8.0)))    !Herbers2026
          t=                                                          !Herbers2026
     $          al4*dff                                               !Herbers2026
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)            !Herbers2026
        end do
      end if
      
      return
      end

C----------------------------------------------------------------------
      subroutine hmulthrr_old(j,gam,f
     $     ,qvk,ruse,a,vr,vi,vor,voi,evalv,ovv,rotm,tori,ifs,it)
C     multiply the complex vector vr,vi by the rigid part of the
C     Hamilton matrix, yielding vor,voi
      implicit none
      include 'iam.fi'
      integer j,gam,f,ifs,it
      integer qvk(DIMTOT,Q_K:Q_V1+DIMTOP-1)
      integer ruse(DIMVV,DIMVV,DIMTOP)
      real*8  a(DIMPAR)
      real*8  vr(DIMTOT),vor(DIMTOT),vi(DIMTOT),voi(DIMTOT)
      real*8            evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)

C     work
      real*8  dk,dj,djj1,dff,t,t1,t2,djjc
      real*8  e1,df,dff1,di,dii1,dg
      real*8  adelk,adelj,ah2,ah3,ahk,ar6,ahjk,ahj
      real*8  alj,aljk,alkj,alk !Herbers2026
      real*8  al2,al3,al4 !Herbers2026
      real*8  DXTERM, DYTERM, DZTERM
      integer ik,ir,ic,ivr,ivc,itop
      integer off,voff
      integer myand
      external myand
      if (size(S_H).gt.DIMTOT) stop 'Dimension Error in HMULTHRR'

      dj=dble(j)
      djj1=dj*(dj+1.0)
      e1=0.0
C     djjc is used for spin rotation coupling to prevent a division by zero
C     for J=0
      djjc=1.0
      if ((ctlint(C_SPIN).ne.0).and.(j.gt.0).and.(f.ge.0)) then
        di=dble(ctlint(C_SPIN))/2.0d0
        dii1=di*(di+1.0)
        df=dble(f)/2.0d0
        dff1=df*(df+1.0)
        djjc=djj1
        dg=dff1-dii1-djj1
        if (ctlint(C_SPIN).gt.1) then
          e1= (0.75*dg*(dg+1.0)-dii1*djj1)
     $         /(2.0*di*(2.0*di-1.0)*djj1*(2.0*dj-1.0)*(2.0*dj+3.0))
        else
          e1=0.0
        end if
      end if

C     the centrifugal distortion parameters
C     Watson A 
      if (ctlint(C_RED).eq.0) then
        adelj=a(P_DJD)
        adelk=a(P_DKD)
        ar6  =0.0d0
        ahj  =a(P_HJD)
        ahjk =a(P_HJKD)
        ah2  =0.0d0
        ahk  =a(P_HKD)
        ah3  =0.0d0
        alj = a(P_LJD)
        aljk = a(P_LJKD)
        alkj = a(P_LKJD)
        alk  = a(P_LKD)
      end if
C     Watson S
      if (ctlint(C_RED).eq.1) then
        adelj=-a(P_DJD)
        adelk=0.0d0
        ar6  =a(P_DKD)
        ahj  =a(P_HJD)
        ahjk =0.0d0
        ah2  =a(P_HJKD)
        ahk  =0.0d0
        ah3  =a(P_HKD)
        alj = a(P_LJD)
        al2 = a(P_LJKD)
        al3 = a(P_LKJD)
        al4 = a(P_LKD)
      end if
C     van Eijck / Typke
      if (ctlint(C_RED).eq.2) then
        adelj=a(P_DJD)
        adelk=0.0d0
        ar6  =a(P_DKD)
        ahj  =0.5d0*a(P_HJD)
        ahjk =0.0d0
        ah2  =0.25d0*a(P_HJKD)
        ahk  =0.0d0
        ah3  =0.125d0*a(P_HKD)
        alj  =0.5*a(P_LJD) ! no further implementation for octic.
      end if

      if ( (a(P_BJ) .ne.0.0).or.(a(P_BK) .ne.0.0).or.
     $     (a(P_DJ) .ne.0.0).or.(a(P_DJK).ne.0.0).or.
     $     (a(P_DK) .ne.0.0).or.(a(P_HJ ).ne.0.0).or.
     $     (a(P_HJK).ne.0.0).or.(a(P_HKJ).ne.0.0).or.
     $     (a(P_HK) .ne.0.0).or.(a(P_QZ) .ne.0.0).or.
     $     (a(P_LJ) .ne.0.0).or.(a(P_LK) .ne.0.0).or.
     $     (a(P_LJJK) .ne.0.0).or.(a(P_LJK) .ne.0.0).or.
     $     (a(P_LKKJ) .ne.0.0).or.
     $     (a(P_E)  .ne.0.0)) then !Herbers2023
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          t=     a(P_BJ) *djj1
     $         + a(P_BK) *dk*dk
     $         - a(P_DJ) *djj1**2
     $         - a(P_DJK)*djj1*dk**2
     $         - a(P_DK) *dk**4
     $         + a(P_HJ) *djj1**3
     $         + a(P_HJK)*(djj1**2)*(dk**2)
     $         + a(P_HKJ)*djj1*(dk**4)
     $         + a(P_HK) *dk**6
     $         + a(P_QZ) * e1 * (3.0*(dk**2)-djj1)
     $         + a(P_CP) * 0.5*dg*(1.0-dk*dk/djjc)
     $         + a(P_CZ) * 0.5*dg*dk*dk/djjc
     $         + a(P_E)  !Herbers2023
     $         + a(P_LJ) *djj1**4                 !Herbers2026
     $         + a(P_LK) *dk**8                   !Herbers2026
     $         + a(P_LJK)*(djj1**2)*(dk**4)       !Herbers2026
     $         + a(P_LJJK)*(djj1**3)*(dk**2)      !Herbers2026
     $         + a(P_LKKJ)*(djj1**1)*(dk**6)      !Herbers2026
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
      
      if (a(P_PZ) .ne.0.0) then
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          t=     a(P_PZ) *dk
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
      
      DZTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)   DZTERM=a(P_DZ1)  !herbers2024
      if (gam.eq.2)   DZTERM=a(P_DZ2)  !herbers2024
      if (gam.eq.3)   DZTERM=a(P_DZ3)  !herbers2024
      if (gam.eq.4)   DZTERM=a(P_DZ4)  !herbers2024
      if (gam.eq.5)   DZTERM=a(P_DZ5)  !herbers2024
      if (gam.eq.6)   DZTERM=a(P_DZ6)  !herbers2024
      if (gam.eq.7)   DZTERM=a(P_DZ7)  !herbers2024
      if (gam.eq.8)   DZTERM=a(P_DZ8)  !herbers2024
      if (gam.eq.9)   DZTERM=a(P_DZ9)  !herbers2024
      if (gam.eq.10)  DZTERM=a(P_DZ10) !herbers2024
      if (gam.eq.11)  DZTERM=a(P_DZ11) !herbers2024
      if (DZTERM.ne.0.0) then          !herbers2024
        do ik=1, size(S_K)             !herbers2024
          dk=dble(qvk(ik,Q_K))         !herbers2024
          t=     DZTERM *dk            !herbers2024
          do ivc=1, size(S_VV)         !herbers2024
            ic=ik+size(S_K)*(ivc-1)    !herbers2024
            vor(ic)=vor(ic)+vr(ic)*t   !herbers2024
            voi(ic)=voi(ic)+vi(ic)*t   !herbers2024
          end do                       !herbers2024
        end do                         !herbers2024
      end if                           !herbers2024
      
      
      
      if ((ctlint(C_RED).eq.2).and. 
     $     ((ar6.ne.0.0d0).or.(ah2.ne.0.0d0))) then
        do ik=1, size(S_K)
          dk=dble(qvk(ik,Q_K))
          dff=  (djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         +(djj1-dk*(dk-1.0))*(djj1-(dk-1.0)*(dk-2.0))
     $         -2.0d0*(djj1-dk**2)**2
          t=  (ar6 + ah2*djj1)*dff
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)=vor(ic)+vr(ic)*t
            voi(ic)=voi(ic)+vi(ic)*t
          end do
        end do
      end if
        
C     real off diagonal k/k+1 
      if ((a(P_QXZ).ne.0.0)) then
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=(1.0+2.0*dk)*dsqrt(djj1-dk*(dk+1.0))
          t=    e1* a(P_QXZ) *dff  
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

C     real off diagonal k/k+1 ! in XIAM P_PX is real and P_PY gives imaginary, this is flipped convention from Gordy 2.63
      if ((a(P_PX).ne.0.0).or.(a(P_DZX).ne.0.0).or.
     $     (a(P_DZXJ).ne.0.0).or.(a(P_DZXK).ne.0.0)) then!Adding Dab, DabJ, DabK here
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=  a(P_PX)*dff
     $       +a(P_DZX)*(dk+(dk+1.0))*dff
     $       +a(P_DZXJ)*(dk+(dk+1.0))*djj1*dff
     $       +a(P_DZXK)*(dk**3+(dk+1.0)**3)*dff     
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if
      
      DXTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)  DXTERM=a(P_DX1)   !herbers2024
      if (gam.eq.2)  DXTERM=a(P_DX2)   !herbers2024
      if (gam.eq.3)  DXTERM=a(P_DX3)   !herbers2024
      if (gam.eq.4)  DXTERM=a(P_DX4)   !herbers2024
      if (gam.eq.5)  DXTERM=a(P_DX5)   !herbers2024
      if (gam.eq.6)  DXTERM=a(P_DX6)   !herbers2024
      if (gam.eq.7)  DXTERM=a(P_DX7)   !herbers2024
      if (gam.eq.8)  DXTERM=a(P_DX8)   !herbers2024
      if (gam.eq.9)  DXTERM=a(P_DX9)   !herbers2024
      if (gam.eq.10) DXTERM=a(P_DX10)   !herbers2024
      if (gam.eq.11) DXTERM=a(P_DX11)   !herbers2024
      if (DXTERM.ne.0.0) then          !herbers2024
C       real off diagonal k/k+1 
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=  DXTERM *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if      

C     imaginaer off diagonal k/k+1 
      if ((a(P_QYZ).ne.0.0)) then
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=(1.0+2.0*dk)*dsqrt(djj1-dk*(dk+1.0))
          t=    e1* a(P_QYZ) *dff *(-1) ! Sign change added by Sven 2024
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if

C     imaginaer off diagonal k/k+1 
      if ((a(P_PY).ne.0.0)) then
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=   a(P_PY) *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if
      
      DYTERM=0.0                       !herbers2024 I added these to use PxPyPz sperate for the various S1,S2,S3,S4,S5 without the need of defining a new set of constants.
      if (gam.eq.1)  DYTERM=a(P_DY1)   !herbers2024
      if (gam.eq.2)  DYTERM=a(P_DY2)   !herbers2024
      if (gam.eq.3)  DYTERM=a(P_DY3)   !herbers2024
      if (gam.eq.4)  DYTERM=a(P_DY4)   !herbers2024
      if (gam.eq.5)  DYTERM=a(P_DY5)   !herbers2024
      if (gam.eq.6)  DYTERM=a(P_DY6)   !herbers2024
      if (gam.eq.7)  DYTERM=a(P_DY7)   !herbers2024
      if (gam.eq.8)  DYTERM=a(P_DY8)   !herbers2024
      if (gam.eq.9)  DYTERM=a(P_DY9)   !herbers2024
      if (gam.eq.10) DYTERM=a(P_DY10)   !herbers2024
      if (gam.eq.11) DYTERM=a(P_DY11)   !herbers2024
      if (DYTERM.ne.0.0) then          !herbers2024
C       imaginaer off diagonal k/k+1 
        off=1
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt(djj1-dk*(dk+1.0))
          t=   DYTERM *dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if
      
C     real off diagonal k/k+2 
      if ((a(P_BD).ne.0.0).or.(adelj.ne.0.0).or.
     $    (adelk.ne.0.0).or.(a(P_QD).ne.0.0).or.
     $    (ahj.ne.0.0).or.(a(P_CD).ne.0.0).or.
     $    (ahjk.ne.0.0).or.(ahk.ne.0.0).or.(alj.ne.0.0).or.
     $    (aljk.ne.0.0).or.(alkj.ne.0.0).or.(alk.ne.0.0)
     $     ) then
        off=2
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
          t=     a(P_BD) *dff
     $         - adelj   *2.0d0*dff*djj1
     $         - adelk   *dff*((dk+2.0d0)**2+dk**2)
     $         + ahj*2.0d0*dff*djj1**2
     $         + ahjk*dff*((dk+2.0d0)**2+dk**2)*djj1
     $         + ahk*dff*((dk+2.0d0)**4+dk**4)
     $         + a(P_QD) *dff*e1
     $         + a(P_CD) *0.5*dg*dff/djjc
     $         + alj*2.0d0*dff*djj1**3                        !Herbers2026
     $         + aljk*dff*((dk+2.0d0)**2+dk**2)*djj1**2       !Herbers2026
     $         + alkj*dff*((dk+2.0d0)**4+dk**4)*djj1          !Herbers2026
     $         + alk*dff*((dk+2.0d0)**6+dk**6)                !Herbers2026
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

      if ((ctlint(C_RED).eq.2).and.(ah3.ne.0.0d0)) then
        off=2
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
     $         *((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         +(djj1-dk*(dk-1.0))*(djj1-(dk-1.0)*(dk-2.0))
     $         +(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0)))
          t=  ah3*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

C     imaginaer off diagonal k/k+2 
      if ((a(P_QXY).ne.0.0)) then
        off=2
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=0.5d0*dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0)))
          t=    
     $         + a(P_QXY)*dff*e1*2.0*(-1.0) ! Sign change Sven 2024
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,1)
        end do
      end if
 
C     Watson S off diagonal k/k+4 (evtl. change dff)
      if ((ar6.ne.0.0d0).or.(ah2.ne.0.0d0).or.(al2.ne.0.0d0)) then
        off=4
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))) 
          t=   
     $          ar6*dff + ah2*djj1*dff + al2*djj1**2*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if

C     Watson S off diagonal k/k+6 
      if ((ah3.ne.0.0d0).or.(al3.ne.0.0)) then
        off=6
        do ik=1, size(S_K)-off
          dk=dble(qvk(ik,Q_K))
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))
     $         *(djj1-(dk+4.0)*(dk+5.0))*(djj1-(dk+5.0)*(dk+6.0))) 
          t=   
     $          ah3*dff+al3*djj1*dff
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)
        end do
      end if
      
C     Watson S off diagonal k/k+8                                     !Herbers2026
      if (al4.ne.0.0) then                                            !Herbers2026
        off=8                                                         !Herbers2026
        do ik=1, size(S_K)-off                                        !Herbers2026
          dk=dble(qvk(ik,Q_K))                                        !Herbers2026
          dff=dsqrt((djj1-dk*(dk+1.0))*(djj1-(dk+1.0)*(dk+2.0))       !Herbers2026
     $         *(djj1-(dk+2.0)*(dk+3.0))*(djj1-(dk+3.0)*(dk+4.0))     !Herbers2026
     $         *(djj1-(dk+4.0)*(dk+5.0))*(djj1-(dk+5.0)*(dk+6.0))     !Herbers2026
     $         *(djj1-(dk+6.0)*(dk+7.0))*(djj1-(dk+7.0)*(dk+8.0)))    !Herbers2026
          t=                                                          !Herbers2026
     $          al4*dff                                               !Herbers2026
          call vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,0)            !Herbers2026
        end do
      end if
      
      
      return
      end

C----------------------------------------------------------------------



      subroutine vadd(ik,off,gam,t,tori,qvk,vr,vi,vor,voi,ri)
      implicit none
      include 'iam.fi'
      integer ik,off,gam,ri
      real*8  t
      integer qvk(DIMTOT,Q_K:Q_V1+DIMTOP-1)
      real*8  vr(DIMTOT),vor(DIMTOT),vi(DIMTOT),voi(DIMTOT)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer ivr,ivc,ir,ic,itop,voff
      real*8  t1,t2
      integer myand
      external myand
      
      if (myand(ctlint(C_WOODS),32).ne.0) then
        do ivr=1, size(S_VV)
          do ivc=1, size(S_VV)
            ir=(ivr-1)*size(S_K)+ik
            ic=(ivc-1)*size(S_K)+ik
            t1=1.0d0
            t2=1.0d0
            do itop=1, ctlint(C_NTOP)
              voff=size(S_MINV+itop)-1
              t1=t1*tori
     $             (qvk(ir,Q_K),qvk(ic+off,Q_K)
     $             ,qvk(ir,Q_V+itop)-voff
     $             ,qvk(ic+off,Q_V+itop)-voff
     $             ,gamma(gam,itop),itop)
              t2=t2*tori
     $             (qvk(ir+off,Q_K),qvk(ic,Q_K)
     $             ,qvk(ir+off,Q_V+itop)-voff
     $             ,qvk(ic,Q_V+itop)-voff
     $             ,gamma(gam,itop),itop)
            end do  
            if (ri.eq.0) then
              vor(ir)    =vor(ir    )+vr(ic+off)*t*t2
              vor(ir+off)=vor(ir+off)+vr(ic    )*t*t1
              voi(ir)    =voi(ir    )+vi(ic+off)*t*t2
              voi(ir+off)=voi(ir+off)+vi(ic    )*t*t1
            else
              vor(ir)    =vor(ir    )+vi(ic+off)*t*t1
              vor(ir+off)=vor(ir+off)-vi(ic    )*t*t2
              voi(ir)    =voi(ir    )-vr(ic+off)*t*t1
              voi(ir+off)=voi(ir+off)+vr(ic    )*t*t2
            end if
          end do
        end do
      end if
      if (myand(ctlint(C_WOODS),32).eq.0) then
        if (ri.eq.0) then
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)    =vor(ic    )+vr(ic+off)*t
            vor(ic+off)=vor(ic+off)+vr(ic    )*t
            voi(ic)    =voi(ic    )+vi(ic+off)*t
            voi(ic+off)=voi(ic+off)+vi(ic    )*t
          end do
        else
          do ivc=1, size(S_VV)
            ic=ik+size(S_K)*(ivc-1)
            vor(ic)    =vor(ic    )+vi(ic+off)*t
            vor(ic+off)=vor(ic+off)-vi(ic    )*t
            voi(ic)    =voi(ic    )-vr(ic+off)*t
            voi(ic+off)=voi(ic+off)+vr(ic    )*t
          end do
        end if
      end if
      return
      end
C----------------------------------------------------------------------
      subroutine esave(j,gam,f,ib,f1,qmk,eval,eused)
      implicit none
      include 'iam.fi'
      real*8  eval(DIMTOT)
      integer j,gam,f,ib,qmk(DIMTOT,DIMQLP), eused(DIMTOT)
      integer f1 ! f1 is used here, to check which qlin it belongs to.
      integer il,ie,iq,vok
      integer use_t,use_tj
      integer myand,myor
      external myand,myor

      do ie=1, size(S_H)
        eused(ie)=0
      end do
        
      do iq=1,2
        do il=1, ctlint(C_NDATA) !compare with all data.
          if ((j.eq.qlin(il,Q_J,iq))
     $           .and.  (f.eq.qlin(il,Q_F,iq))!upper and lower is iq.
     $           .and.  (f1.eq.qlin(il,Q_F1,iq))
     $           .and. (ib.eq.qlin(il,Q_B,iq))
     $           .and.(gam.eq.qlin(il,Q_S,iq))) then 
            use_tj=myand(qlin(il,Q_STAT,iq),2)
            use_t =myand(qlin(il,Q_STAT,iq),4)
            do ie=1, size(S_H)
              if (
     $             ((qmk(ie,Q_T).eq.qlin(il,Q_T,iq))
     $             .and.(use_t.ne.0))
     $             .or.
     $             ((qmk(ie,Q_TJ).eq.qlin(il,Q_TJ,iq))
     $             .and.(use_tj.ne.0))
     $             .or.
     $             ((qmk(ie,Q_K).eq.qlin(il,Q_K,iq))
     $             .and.(use_t.eq.0).and.(use_tj.eq.0))
     $             ) then
                vok=0
                if ((qmk(ie,Q_V1).eq.qlin(il,Q_V1,iq)).or.(use_t.ne.0))
     $               vok=ctlint(C_NTOP)
                if ((myand(qlin(il,Q_STAT,iq),1).eq.0)
     $               .and.(vok.eq.ctlint(C_NTOP))) then 
                  dnv(il,NV_ENG,iq)=eval(ie)
                  qlin(il,Q_STAT,iq)=myor(qlin(il,Q_STAT,iq),1)
                  qlin(il,Q_TJ,iq)=qmk(ie,Q_TJ)
                  qlin(il,Q_T ,iq)=qmk(ie,Q_T )
                  qlin(il,Q_K, iq)=qmk(ie,Q_K )
                  qlin(il,Q_K2,iq)=qmk(ie,Q_K2)
                  qlin(il,Q_GK,iq)=qmk(ie,Q_GK)
                  qlin(il,Q_V1,iq)=qmk(ie,Q_V1)
C                  do itop=1, ctlint(C_NTOP)
C                    qlin(il,Q_V+itop,iq)=qmk(ie,Q_V+itop)
C                  end do
                  if (dln(il,LN_ERR).ne.NOFIT)
     $                 eused(ie)=myor(eused(ie),1)
                end if
              end if
            end do
          end if
        end do
      end do
      return
      end
C----------------------------------------------------------------------
C----------------------------------------------------------------------
      subroutine devsave(j,gam,f,ib,t,ifit,dfit,dedp,palc,pali)
      implicit none
      include 'iam.fi'
      integer j,gam,f,t,ib,ifit(DIMPAR),dfit(DIMFIT)
      real*8  dedp(DIMPAR)
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)

      integer il,ip,iq,ifp
      
      do ifp=1, DIMFIT
        if (dfit(ifp).gt.0) then
          do il=1, ctlint(C_NDATA) 
            do iq=1,2
              if (      (j  .eq.qlin(il,Q_J,iq))
     $             .and.(f  .eq.qlin(il,Q_F,iq))
     $             .and.(gam.eq.qlin(il,Q_S,iq))
     $             .and.(t  .eq.qlin(il,Q_T,iq))
     $             .and.(ib .eq.qlin(il,Q_B,iq))) then
c                do ip=1, DIMPAR
c                  dnv(il,ifp+1,iq)=dnv(il,ifp+1,iq)
c     $                 +palc(ifp,ip)*dedp(ip)
c                end do
                do ip=1, pali(ifp,0,1)
                  if (pali(ifp,ip,2).eq.ib)
     $                 dnv(il,ifp+NV_DEF,iq)=dnv(il,ifp+NV_DEF,iq)
     $                 +palc(ifp,ip)*dedp(pali(ifp,ip,1))
                end do
              end if
            end do
          end do
        end if
      end do
      return
      end

C----------------------------------------------------------------------
      subroutine assgn(j,gam,f,ib,f1,h,eval,qvk,qmv,qmk,qv,fistat)
C     qv: Quantum No.s from buildvjk ( 1  1  1  1  1  2  2  2  2  2 ..)
C     qk: Quantum No.s from buildvjk (-2 -1  0  1  2 -2 -1  0  1  2 ..)
C     qmv: Qunatum No.s from buildm  ( 0  3 -3  6  ...) or ( 1 -2  4 -5 ...)

      implicit none
      include 'iam.fi'
      integer j,gam,f,ib,fistat,f1!f1 is not used for anything in the moment, also not printing.
      real*8  h(DIMTOT,DIMTOT),eval(DIMTOT)
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP)
      integer qmv(DIMV),qmk(DIMTOT,DIMQLP),qv(DIMTOT)

C     work
      integer i,ik,icc,iv,goodn,itop,ofistat
      integer i1,i2,iv1,iv2,v1,v2,ik1,ik2
C      integer bestk(DIMTOT),bestv(DIMTOT,DIMTOP),bestvv(DIMTOT)
      integer bestk(DIMTOT),bestvv(DIMTOT)
      integer bestvk(1)
C      integer scndk(DIMTOT),scndv(DIMTOT,DIMTOP),scndvv(DIMTOT)
      integer scndk(DIMTOT),scndvv(DIMTOT)
      integer scndvk(1)
      integer qtp(DIMVV),qvp(DIMVV)
      real*8  besth(DIMTOT),scndh(DIMTOT)
      real*8  besth1(1),scndh1(1)
      real*8  ksum(DIM2J1,DIMTOT)
C      real*8  vsum(DIMV,DIMTOT,DIMTOP)
      real*8  vsum(DIMVV,DIMTOT)
      real*8  vvsum(DIMVV,DIMTOT)
      real*8  difk,difv(DIMTOP),degn
      character*1 sgnchr
      logical kdegen
      data  ofistat /0/
C      save /xx/ofistat
C     
C      integer found(DIMTOT),fp,bestgk,bestfp,k
C
      if (size(S_K).gt.DIM2J1) stop ' ERROR: DIM 2J1 exceeded'
      degn=2.0d-4
      kdegen=.false.
      do icc=1, size(S_H)
        qmk(icc,Q_K)=NaQN
        qmk(icc,Q_T)=NaQN
        qmk(icc,Q_K2)=NaQN
        qmk(icc,Q_GK)=-200
      end do
      do iv=1, size(S_VV)
        qtp(iv)=0
        qvp(iv)=0
      end do
  
C     init the ksum and vsum matrix
      do ik=1, size(S_K)
        do i=1, size(S_H)
          ksum(ik,i)=0.0
        end do
      end do
c      do itop=1, ctlint(C_NTOP)
c        do i=1, size(S_H)
c          do iv=1, size(S_V+itop)
c            vsum(iv,i,itop)=0.0
c          end do
c        end do 
c      end do
      do i=1, size(S_H)
        do iv=1, size(S_VV)
          vvsum(iv,i)=0.0
          vsum(iv,i)=0.0
        end do
      end do

      if (ctlint(C_NTOP).gt.0) then
        do i=1, size(S_H)
          do iv=1, size(S_VV)
            do ik=1, size(S_K)
              i1=(iv-1)*size(S_K)+ik
              vvsum(iv,i)=vvsum(iv,i)+h(i1,i)**2
            end do
          end do
        end do
      end if
c      write(*,'(A,2I3)') ' vvsum: j,gam',j,gam
c      do i=1, size(S_H)
c        write(*,'(F10.5,$)') (vvsum(iv,i),iv=1, size(S_VV))
c         write(*,*)
c      end do

      if (ctlint(C_NTOP).eq.2) then
        do iv1=1, size(S_VV)
          v1=qvv(iv1,1,ib)
          v2=qvv(iv1,2,ib)
          if (v1.gt.v2) then
c     vfnd=.false.
            do iv2=1, size(S_VV)
              if ((qvv(iv2,2,ib).eq.v1).and.(qvv(iv2,1,ib).eq.v2)) then
C     vfnd=.true.
                do ik=1, size(S_K)
                  i1=(iv1-1)*size(S_K)+ik
                  i2=(iv2-1)*size(S_K)+ik
                  do i=1, size(S_H)
                    vsum(iv1,i)=vsum(iv1,i)+0.5*(h(i1,i)+h(i2,i))**2
                    vsum(iv2,i)=vsum(iv2,i)+0.5*(h(i1,i)-h(i2,i))**2
                  end do
                end do
                do i=1, size(S_H)
                  if (abs(vsum(iv1,i)-vsum(iv2,i))
     $                 .gt.abs(vvsum(iv1,i)-vvsum(iv2,i))) then
                    vvsum(iv1,i)=vsum(iv1,i)                
                    vvsum(iv2,i)=vsum(iv2,i)                
                  end if
                end do
              end if
            end do
          end if
        end do
      end if
c      write(*,'(A,2I3)') ' vsum: j,gam,',j,gam
c      do i=1, size(S_H)
c        write(*,'(I2,(F10.5,$))') i,
c     $       (vsum(iv,i),iv=1, size(S_VV))
c        write(*,*)
c      end do
c      write(*,'(A,2I3)') ' vvsum: j,gam,',j,gam
c      do i=1, size(S_H)
c        write(*,'(I2,F15.5,(F10.5,$))') i,
c     $       eval(i),(vvsum(iv,i),iv=1, size(S_VV))
c        write(*,*)
c      end do

c        if ((v1.eq.v2).or.(not.vfnd)) then
c          do ik=1, size(S_K)
c            i1=(iv1-1)*size(S_K)+ik
c            vvsum(iv1,i)=vvsum(iv1,i)+(h(i1,i))**2
c          end do
c        end if  
c      end do
c      kof=j+1
c      do iv=1, size(S_H)
c        do i=1, size(S_H)
c          do itop=1, ctlint(C_NTOP)
c            vsum(qvk(iv,Q_V+itop),i,itop)
c     $           = vsum(qvk(iv,Q_V+itop),i,itop)+h(iv,i)**2
c          end do
c          ksum(qvk(iv,Q_K)+kof,i)=ksum(qvk(iv,Q_K)+kof,i)+h(iv,i)**2
c          vvsum(qv(iv),i)       =vvsum(qv(iv),i)+h(iv,i)**2
c        end do  
c      end do
c      kof=j+1

c     if (gam.eq.0) then
c     do i=1, size(S_H)
c       do ik1=1, size(S_K)/2
c         ik2=size(S_K)-ik1+1
c         do iv=1, size(S_VV)
c           i1=(iv-1)*size(S_K)+ik1
c           i2=(iv-1)*size(S_K)+ik2
c           ksum(ik1,i)=ksum(ik1,i)+0.5d0*(h(i1,i)+h(i2,i))**2
c           ksum(ik2,i)=ksum(ik2,i)+0.5d0*(h(i1,i)-h(i2,i))**2
c         end do
c       end do  
c       ik=size(S_K)/2+1
c       do iv=1, size(S_VV)
c         i1=(iv-1)*size(S_K)+ik
c         ksum(ik,i)=ksum(ik,i)+h(i1,i)**2
c       end do
c     end do
c     else
      do i=1, size(S_H)
        do ik=1, size(S_K)
          do iv=1, size(S_VV)
            i1=(iv-1)*size(S_K)+ik
            ksum(ik,i)=ksum(ik,i)+h(i1,i)**2
          end do
        end do  
      end do
c     end if
      call maxof(ksum,DIM2J1,DIMTOT,size(S_K),1,size(S_H),bestk,scndk
     $     ,besth,scndh)
      call maxof(vvsum,DIMVV,DIMTOT,size(S_VV)
     $       ,1,size(S_H),bestvv,scndvv,besth,scndh)
c      do itop=1, ctlint(C_NTOP)
c        call maxof(vsum(1,1,itop),DIMV,DIMTOT,size(S_V+itop)
c     $       ,1,size(S_H),bestv(1,itop),scndv(1,itop),besth,scndh)
c      end do

      qvp(bestvv(1))=1
      do icc=1, size(S_H)-1
        if ((bestvv(icc+1).ne.bestvv(icc)).and.
     $       (qvp(bestvv(icc+1)).eq.0)) then
          qvp(bestvv(icc+1))=qvp(bestvv(icc))+1
        end if
      end do
      if (((ctlint(C_EVAL).gt.0).and.(fistat.eq.0))
     $     .or.(ctlint(C_EVAL).gt.1)) then
        if (ofistat.ne.fistat) then
          write(20,'(A,I3)') ' *********** FISTAT =',fistat
          ofistat=fistat
        end if
        write(20,'(/,4A3,2X,4A3,A15,4A)')
     $         'J','S','B','F','T','t','K','V',' Energy/GHz '
     $        ,'   best K(s)      vec ','  V   vec'
     $        ,' 2nd.K vec',' 2nd.V vec'
      endif 
      do icc=1, size(S_H)
        difk=1.0
        if (size(S_K).gt.1) then
          difk=(dabs(ksum(bestk(icc),icc)-ksum(scndk(icc),icc)))
        end if

        do itop=1, ctlint(C_NTOP)        
          difv(itop)=1.0
        end do
c        do itop=1, ctlint(C_NTOP)
c          if (size(S_V+itop).gt.1) then
c            difv(itop)=(dabs(vsum(bestv(icc,itop),icc,itop)
c     $         - vsum(scndv(icc,itop),icc,itop)))
c          end if
c        end do
        qmk(icc,Q_J)=j
        qmk(icc,Q_S)=gam
        qmk(icc,Q_F)=f
        qmk(icc,Q_B)=ib
        if (difk.lt.degn) then
          kdegen=.true.
          call maxof(h,DIMTOT,DIMTOT,size(S_H),icc,1,bestvk,scndvk,
     $         besth1,scndh1)
          if (abs(qvk(bestvk(1),Q_K)).ne.abs(qvk(bestk(icc),Q_K)))
     $         write(*,*) 'ERROR in assgn: K trouble',j,gam,icc 
          qmk(icc,Q_K)=sign(qvk(bestvk(1),Q_K),
     $         int(besth1(1)*scndh1(1)*10000.0))
          qmk(icc,Q_K2)=qmk(icc,Q_K)
          qmk(icc,Q_GK)=1000
          qtp(bestvv(icc))=qtp(bestvv(icc))+1
          qmk(icc,Q_TJ)=qtp(bestvv(icc))
          qmk(icc,Q_V1)=qvp(bestvv(icc))
          qmk(icc,Q_T )=icc
          sgnchr=' '
          if (qmk(icc,Q_K).gt.0) sgnchr='+'
          if (qmk(icc,Q_K).lt.0) sgnchr='-'
          if (((ctlint(C_EVAL).gt.0).and.(fistat.eq.0))
     $         .or.(ctlint(C_EVAL).gt.1)) then 
            write(20,
     $         '(4i3,2X,2i3,i2,A,i3,F15.8,1X,2I3,2F7.3,I4,F6.3,10X,I4,
     $           F6.3)')
     $         j,gam,ib,f,qmk(icc,Q_T),qmk(icc,Q_TJ),abs(qmk(icc,Q_K))
     $         ,sgnchr,qmk(icc,Q_V1),eval(icc)
     $         ,qvk(bestvk(1),Q_K) ,qvk(scndvk(1),Q_K)
     $         ,besth1(1),scndh1(1) 
     $         ,bestvv(icc),sqrt(vvsum(bestvv(icc),icc))
     $         ,scndvv(icc),sqrt(vvsum(scndvv(icc),icc))
C     $         ,(bestv(icc,itop),itop=1,ctlint(C_NTOP))
          end if
        else
          goodn=int(dabs(1.0-ksum(scndk(icc),icc)/ksum(bestk(icc),icc))
     $         *100.0)
          if (qmk(icc,Q_GK).ne.-200) write(0,'(A,5I3)')
     $         ' Assignment Warning',j,gam,ib,f,icc
          qmk(icc,Q_GK)=goodn
          qmk(icc,Q_K) =qvk(bestk(icc),Q_K) ! this is not correct, but it will
          qmk(icc,Q_K2)=qvk(scndk(icc),Q_K) ! work. 
          qtp(bestvv(icc))=qtp(bestvv(icc))+1
          qmk(icc,Q_TJ)=qtp(bestvv(icc))
          qmk(icc,Q_V1)=qvp(bestvv(icc))
          qmk(icc,Q_T )=icc
          if (((ctlint(C_EVAL).gt.0).and.(fistat.eq.0))
     $         .or.(ctlint(C_EVAL).gt.1)) then
            write(20,
     $       '(4i3,2X,4i3,F15.8,1X,I6,F14.3,I4,F6.3,I4,F6.3,I4,F6.3)')
     $         j,gam,ib,f,qmk(icc,Q_T),qmk(icc,Q_TJ),qmk(icc,Q_K)
     $         ,qmk(icc,Q_V1),eval(icc)
     $         ,qvk(bestk(icc),Q_K),sqrt(ksum(bestk(icc),icc))
     $         ,bestvv(icc),sqrt(vvsum(bestvv(icc),icc))
     $         ,qvk(scndk(icc),Q_K),sqrt(ksum(scndk(icc),icc))
     $         ,scndvv(icc),sqrt(vvsum(scndvv(icc),icc))
C     $         ,(bestv(icc,itop),vsum(bestv(icc,itop),icc,itop)
C     $         ,itop=1,ctlint(C_NTOP))
          end if
        end if         
      end do
      if (gam.gt.0) then
        if (kdegen) gamma(gam,0)=1
      end if
      return

      end

C----------------------------------------------------------------------
       subroutine rotate(d,beta,j,oldj)
       implicit none
       include 'iam.fi'
       real*8  d(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2)
       real*8  beta
       integer j,oldj,ij
       if ((j.eq.0).or.(j.eq.1)) then
         call rotat1(d,beta,j)
         oldj=j
       end if
       if (oldj.eq.j) return
       if (oldj.gt.j) then
         do ij=1, j 
           call rotat1(d,beta,ij)
         end do
       end if
       if (oldj.lt.j) then
         do ij=oldj+1, j 
           call rotat1(d,beta,ij)
         end do
       end if
       oldj=j
       return
       end

C----------------------------------------------------------------------
      subroutine nqcmat_ir(jselect,gam,f,ib,evalv,ovv,rotm,rott
     $ ,tori,atot,qmv,h_2,qvks,qmks,qmvs,qvs,ruse)
      ! The original exact matrix built was removed and instead this routine is used to only fill in semi-rigid rotor matrix elements and internal rotation contribution.
      ! f should be f1 in input.
      implicit none
      include 'iam.fi'
      
      real*8  atot(DIMPAR,DIMVB)

      integer j, gam, f, ib                  
      integer jselect 
      integer sj,ej,cm
      integer oldj(DIMTOP)
      real*8  h(DIMTOT,DIMTOT),evh(DIMTOT)
      real*8  h_2(DIMQ*DIMTOT,DIMQ*DIMTOT)
      real*8  hresort(DIMQ*DIMTOT,DIMQ*DIMTOT) 
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP) 
      real*8  a(DIMPAR)
      real*8  beta_tot
      integer qmv(DIMV)
      integer qmvs(DIMQ,DIMV)
C     quantum numbers
      integer qvk(DIMTOT,Q_K:Q_V+DIMTOP), qmk(DIMTOT,DIMQLP)
      integer qmks(DIMQ,DIMTOT,DIMQLP)
      integer qvks(DIMQ,DIMTOT,Q_K:Q_V+DIMTOP) 
      integer qvs(DIMQ,DIMTOT)
      integer qv(DIMTOT)
C     work
      integer check2
      integer c2o
      integer c2o2
      integer id,ie,i,iv,ik,itop,ivr,ivc,ir,ic,it1,it2,ikr,ikc
      integer ruse(DIMVV,DIMVV,DIMTOP)
      logical complex
      
      complex = .true.

      h_2=0.0
      hresort=0.0

      sj=(2*jselect-(f-ctlint(C_SPIN)))/2!startj offset
      ej=((f+ctlint(C_SPIN))-2*jselect)/2!endj offset
      if ((jselect-sj).lt.0) sj=jselect
      if (f.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs(f+2*(jselect-sj)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       sj=sj-(abs(abs(f-2*(jselect-sj))-ctlint(C_SPIN))/2)
       end if 
      end if
      
      if (f.eq.-1) then
       sj=0
       ej=0
      end if
      !!! for non existing J/F states with J>jselect no matrix has to be built up
      a=atot(:,ib)
      do itop=1,ctlint(C_NTOP) ! resetting oldj to 0.
        oldj(itop)=0           ! resetting oldj to 0.
      end do                   ! resetting oldj to 0.
      

      do j=jselect-sj,jselect+ej
      cm=j-(jselect-sj)!count matrices
      
      if (j.ge.0) then   !START IF HERE 
      do itop=1,ctlint(C_NTOP)
        beta_tot=a(P1_BETA+DIMPIR*(itop-1))
C     $     +log(1.+j*(j+1))*a(P1_BLOGJ+DIMPIR*(itop-1))
     $     +dble((j*(j+1))**1)*a(P1_BETJ1+DIMPIR*(itop-1))
     $     +dble((j*(j+1))**2)*a(P1_BETJ2+DIMPIR*(itop-1))
     $     +dble((j*(j+1))**3)*a(P1_BETJ3+DIMPIR*(itop-1))
     $     +dble((j*(j+1))**4)*a(P1_BETJ4+DIMPIR*(itop-1))
        
          call rotate(rotm(-DIMJ,-DIMJ,1,itop)
     $         ,beta_tot,j,oldj(itop))
      end do

C       if ((abs((j)-df).le.abs(di)).or.(df.eq.-0.5)) then ! For some reason this exception leads to malfunction.
      size(S_K)=2*j+1
      i=0
      do iv=1, size(S_VV)
        do ik=1, size(S_K)
          i=i+1
          qvk(i,Q_K) =ik-j-1
          qv(i)=iv
          do itop=1, ctlint(C_NTOP)
            qvk(i,Q_V+itop)=qvv(iv,itop,ib)
          end do
        end do
      end do
      size(S_H)=i
      
      h(:size(S_H),:size(S_H))=0.0
      
      if (gam.ne.0) then 
        call bld2vjk(j,gam,f
     $         ,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,complex)

      else
        if (size(S_VV).gt.1) stop 'size vv > 1 for rigid rotor!'
      end if


      call addrig(j,gam,f
     $     ,qvk,ruse,h,a,evalv,ovv,rotm,rott,tori,complex)
      if ((abs(2*j-f).le.abs(ctlint(C_SPIN))).or.(f.eq.-1)) then ! J-F > I  
       h_2(1+DIMTOT*cm:1+DIMTOT*cm+size(S_H),       ! same as double do loop above but more compact.
     $    1+DIMTOT*cm:1+DIMTOT*cm+size(S_H))=h(:size(S_H),:size(S_H)) ! same as double do loop above but more compact.
      end if ! J-F > I            


       
       qvks(cm+1,:,:)=qvk
       qmks(cm+1,:,:)=qmk
       qmvs(cm+1,:)=qmv
       qvs(cm+1,:)=qv
        
       end if ! J<0
      end do

        



      check2=0   
      c2o=-1
      c2o2=-1
      do j=jselect-sj,jselect+ej
      
      cm=j-(jselect-sj)!count matrices
      do i=1, 2*j+1                                  
      end do
      hresort(1+check2:2*j+1+check2,1+check2:2*j+1+check2)=
     $ h_2(cm*DIMTOT+1:cm*DIMTOT+2*j+1
     $   ,cm*DIMTOT+1:cm*DIMTOT+2*j+1)  ! Diagonal Blocks dj=0
     
        if (c2o.ne.-1) then
        hresort(1+c2o:2*(j-1)+1+c2o,1+check2:2*j+1+check2)=
     $ h_2((cm-1)*DIMTOT+1:(cm-1)*DIMTOT+2*(j-1)+1
     $   ,cm*DIMTOT+1:cm*DIMTOT+2*j+1)  ! Offdiag 1a dj=1
       hresort(1+check2:2*j+1+check2,1+c2o:2*(j-1)+1+c2o)=
     $ h_2(cm*DIMTOT+1:cm*DIMTOT+2*j+1
     $   ,(cm-1)*DIMTOT+1:(cm-1)*DIMTOT+2*(j-1)+1)  ! Offdiag 1b
        endif
      if (c2o2.ne.-1) then
       hresort(1+c2o2:2*(j-2)+1+c2o2,1+check2:2*j+1+check2)=
     $ h_2((cm-2)*DIMTOT+1:(cm-2)*DIMTOT+2*(j-2)+1
     $   ,cm*DIMTOT+1:cm*DIMTOT+2*j+1)  ! Offdiag 2a dj=2
       hresort(1+check2:2*j+1+check2,1+c2o2:2*(j-2)+1+c2o2)=
     $ h_2(cm*DIMTOT+1:cm*DIMTOT+2*j+1
     $   ,(cm-2)*DIMTOT+1:(cm-2)*DIMTOT+2*(j-2)+1)  ! Offdiag 2b
      endif  
      c2o2=c2o
      c2o=check2
      check2=check2+2*j+1
      end do
      
      h_2=hresort
     
      return 
      end 
C------------------------------------------------------------------  
C------------------------------------------------------------------  
      subroutine nqcmat_NQ2df1eqn(jselect,f1select,f
     $ ,ib,atot,h_3,occupied,n)
C      This subroutine builds the matrix elements for the second nucleus, offdiagonal in f1 by n
C      n=0 means that the elements diagonal in f1 are added.
       implicit none
       include 'iam.fi'
       real*8 h_3(DIMQ2*DIMQ*DIMTOT,DIMQ2*DIMQ*DIMTOT)
       real*8  atot(DIMPAR,DIMVB)
       integer h_sizesr(DIMQ)
       integer sizeprer
       integer sizejr
       integer h_sizesc(DIMQ)
       integer sizeprec
       integer sizejc
       integer totsizef1c
       integer n
       integer f1c, f1r
       integer f1cp1
       integer offsetf1
       integer occupied
       
       integer jselect,f1select,f,ib
       integer jr,jc !j of row, j of column.
       integer jcp1
       integer startjr, endjr, startjc, endjc ! first J and last J in loop
       integer startjcp1, endjcp1
       integer cmr,cmc    !position index in dependence of j (counts matrices)
       integer ir,ic !running indices for rows and columns of a matrix
       integer start_ic, end_ic ! limits the ic loop, since delta K=2 is maximum.
       integer kc,kr ! K quantum numbers for the matrix elements (columns and rows)
       integer t2    !determins the exponent on the (-1) prefactor.
       
       real*8 tj,wsj !output from wigner 3j and wigner 6j soubroutine.
       real*8 tji !the threej symbol involving only I only has to be calculated once and can be kept in memory since start of this subroutine (even earlier actually)
                      !the wsj symbol involving only F,F1,and I can also be calculated once and kept for the elements DIAGONAL in F1, which are calculated in this subroutine.
       real*8 chiq_r(-2:2) ! the various real parts of chiq
       real*8 chiq_i(-2:2) ! the various imaginary parts of chiq
       real*8 q2z, q2d, q2xz, q2yz, q2xy 
       integer q
       
       real*8 totprod    !total value of matrix element
       
       
       f1c=f1select
       f1r=f1select+n*2 !the nth offdiagonal in f1
       f1cp1=f1c+2  !f1+1
      
      
      call threej(ctlint(C_SPIN2), 4, ctlint(C_SPIN2)
     $         ,-ctlint(C_SPIN2), 0 , ctlint(C_SPIN2), tji) !threej for i2 only has to be calculated one times here.
      startjc=(f1c-ctlint(C_SPIN))/2!startj ! for columns
      endjc=(f1c+ctlint(C_SPIN))/2!endj      ! for columns
      startjr=((f1r)-ctlint(C_SPIN))/2!startj ! for rows, offset by two halfs for delf1 eq 1
      endjr=((f1r)+ctlint(C_SPIN))/2!endj      ! for columns
      
      
      startjcp1=(f1cp1-ctlint(C_SPIN))/2 ! these are required to calculate matrix size offset for deltaF1 eq.2
      endjcp1=(f1cp1+ctlint(C_SPIN))/2   ! these are required to calculate matrix size offset for deltaF1 eq.2
      
      !SET this up in dependence of atot!
      q2z  = -atot(P_Q2Z,ib)  ! working hypothesis, sign change required to match previous outputs.
      q2d  = -atot(P_Q2D,ib)  ! working hypothesis, sign change required to match previous outputs.
      q2xz = -atot(P_Q2XZ,ib) ! working hypothesis, sign change required to match previous outputs.
      q2yz = -atot(P_Q2YZ,ib) ! working hypothesis, sign change required to match previous outputs.
      q2xy = atot(P_Q2XY,ib) ! working hypothesis, sign change required to match previous outputs.
                        !q gives the delta K values of the matrix elements.
      chiq_r(0)=q2z
      chiq_i(0)=0.0 !no imaginary part for q=0 (elements diagonal in K must be real)
      
      chiq_r(1)=-(2.0/3.0)**0.5*(q2xz)
      chiq_i(1)=-(2.0/3.0)**0.5*(-q2yz)
      
      chiq_r(-1)=(2.0/3.0)**0.5*(q2xz) ! single sign change in real part
      chiq_i(-1)=(2.0/3.0)**0.5*(q2yz)! double sign change in imaginary part
      
      chiq_r(2) =(1.0/6.0)**0.5*(q2d)
      chiq_i(2) =(1.0/6.0)**0.5*(2*q2xy)
      chiq_r(-2)=chiq_r(2) !no sign change in real part for q=2
      chiq_i(-2)=-chiq_i(2)     
      
      !SET this up in dependence of atot!
      
      if (startjc.lt.0) startjc=0
      if (f1c.lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs(f1c+2*(startjc)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       startjc=(ctlint(C_SPIN)-f1c)/2
       end if 
      end if

      if (startjr.lt.0) startjr=0
      if ((f1r).lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs((f1r)+2*(startjr)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       startjr=(ctlint(C_SPIN)-f1r)/2
       end if 
      end if
      
      if (startjcp1.lt.0) startjcp1=0
      if ((f1cp1).lt.ctlint(C_SPIN))then ! if F is smaller than I, then there are J states that can not exist
       if (abs((f1cp1)+2*(startjcp1)).lt.ctlint(C_SPIN))then !The non existing states are characterized by F+J < I.
       startjcp1=(ctlint(C_SPIN)-f1cp1)/2
       end if 
      end if
      
      
      
      offsetf1=0   !the offset will be used to find the correct position for the diagonal and offdiagonal delta f1 blocks.
      if (n.ge.1) then
       do jc=startjc,endjc
       offsetf1=offsetf1+2*jc+1
       end do
       if (n.ge.2) then
       do jcp1=startjcp1,endjcp1
        offsetf1=offsetf1+2*jcp1+1
       end do
       end if
      end if 

      
      do jr=startjr,endjr
       cmr=jr-startjr!count matrices starting at 0
       sizejr=2*jr+1 
       h_sizesr(cmr+1)=sizejr
       sizeprer=sum(h_sizesr(:cmr+1))-sizejr
       do jc=startjc,endjc
        cmc=jc-startjc!count matrices starting at 0
        sizejc=2*jc+1 
        h_sizesc(cmc+1)=sizejc
        sizeprec=sum(h_sizesc(:cmc+1))-sizejc
C
CCC           All elements delta J in one clause.
CC
        if ((n.ne.0).or.(jr.ge.jc)) then !for elements diagonal in F1, elements with jr lower than jc would fall into the imaginary upper right triangle.
        if (abs(jc-jr).le.2) then
        do kr=-jr,jr
         do kc=-jc,+jc
           do q=-2,2    !the 5 cases for q.   
            if ((-kc-q+kr).eq.0) then !triangle condition for tj symbol.
             
             call sixj(f,ctlint(C_SPIN2),f1c,4
     $                      ,f1r,ctlint(C_SPIN2),wsj)
             totprod=wsj*((f1c+1)*(f1r+1))**0.5 !sj1*((2f1+1)*(2f1'+1))**0.5
 
             call sixj(2*jc,f1c,ctlint(C_SPIN)
     $                               ,f1r,2*jr,4,wsj)
             totprod=totprod*wsj*((2*jc+1)*(2*jr+1))**0.5 !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5
             
             call threej(2*jc, 4, 2*jr, -2*kc, -2*q , 2*kr, tj)
             totprod=totprod*tj !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5*tj1
 
             
             tj=tji
             totprod=totprod/tj
             totprod=totprod/4
             t2=2*kc+2+ctlint(C_SPIN)+f1r
     $                       +ctlint(C_SPIN2)+f1r+f
             if(mod(t2,2).ne.0) then
              write(0,*) "t2*2 is not even, this can not be. dj0nq2", t2
              stop
             end if
             t2=t2/2
             totprod=totprod*(-1)**t2 !(-1)**t2*sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5*tj1/tj2/4 ! only thing missing at this point is the chis.
             
             ir = sizeprer + kr + jr + 1
             ic = sizeprec + kc + jc + 1
C             if ((occupied+offsetf1+ir.eq.10).and.
C     $               (occupied+ic.eq.1).and.(f1c.eq.1)) then
C              write(0,*) f, f1c, jc, kc, f1r, jr, kr, totprod*chiq_r(-q)
C              write(0,*) occupied+offsetf1+ir, occupied+ic
C              write(0,*) startjr, startjc, sizeprec
C             end if
             h_3(occupied+offsetf1+ir,occupied+ic)=totprod*chiq_r(-q) ! real matrix elements, lower triangle and diagonal 
             ! offsetf1 is here used as an offset, since it is equal to the lower state matrix dimension.
             ! occupied jumps the f1select part of the matrix.
             if (q.ne.0) then
              h_3(occupied+ic,occupied+offsetf1+ir)=totprod*chiq_i(-q) !upper triangle, imaginary.
             end if
            end if !check only those where non zero matrix elements are expected (see also the first wigner 3jsymbol for this condition)
           end do !q
          end do !ic
         end do   !ir 
        end if ! Jpm2 condition
        end if !delta F1=0 condition for imaginary triangle.
        
       end do!!! jc loop ends.
      end do!!! jr loop ends.
      return 
      end 
      
      
C------------------------------------------------------------------  
      subroutine nqcmat_NQ1(jselect,f1select,ib,h_2,atot,two_I)
C      Same as nqcmat_NQ2, but for testing if the old implementation can be changed with this new one for a single nucleus (general expression vs explicit matrix elelemts)
C      arguments change a bit compared to NQ2, since the wigner symbols also change, ie f is not used anymore, also i give spin and quadrupole coupling tensor components as 
C      arguments here, so I can swap for testing purposes
       implicit none
       include 'iam.fi'
       real*8 h_2(DIMQ*DIMTOT,DIMQ*DIMTOT)
       real*8  atot(DIMPAR,DIMVB)
       integer h_sizes(DIMQ)
       integer sizepre
       integer sizej
       integer jselect,f1select,ib
       integer j
       integer sj,ej !offset of Jmin and Jmax relative to Jselect
       integer cm    !position index in dependence of j (counts matrices)
       integer ir,ic !running indices for rows and columns of a matrix
       integer start_ic, end_ic ! limits the ic loop, since delta K=2 is maximum.
       integer kc,kr ! K quantum numbers for the matrix elements (columns and rows)
       integer t2    !determins the exponent on the (-1) prefactor.
       integer two_I !replaces ctlint(C_SPIN) in this subroutine
       
       real*8 di,dii1,df,dff1 ! first nucleus
       real*8 tj,wsj !output from wigner 3j and wigner 6j soubroutine.
       real*8 tji!the threej symbol involving only I only has to be calculated once and can be kept in memory since start of this subroutine (even earlier actually)
                      !the wsj symbol involving only F,F1,and I can also be calculated once and kept for the elements DIAGONAL in F1, which are calculated in this subroutine.
       real*8 chiq_r(-2:2) ! the various real parts of chiq
       real*8 chiq_i(-2:2) ! the various imaginary parts of chiq
       real*8 qz, qd, qxz, qyz, qxy 
       integer q
       
       real*8 totprod    !total value of matrix element
       
      h_2=0.0
      
      call threej(two_I, 4, two_I
     $         ,-two_I, 0 , two_I, tji) !threej for i2 only has to be calculated one times here.
      
      sj=(2*jselect-(f1select-two_I))/2!startj offset
      ej=((f1select+two_I)-2*jselect)/2!endj offset
      
      
      !SET this up in dependence of atot!
      qz  = -atot(P_QZ,ib)  ! working hypothesis, sign change required to match previous outputs.
      qd  = -atot(P_QD,ib)  ! working hypothesis, sign change required to match previous outputs.
      qxz = -atot(P_QXZ,ib) ! working hypothesis, sign change required to match previous outputs.
      qyz = -atot(P_QYZ,ib) ! working hypothesis, sign change required to match previous outputs.
      qxy = atot(P_QXY,ib) ! working hypothesis, sign change required to match previous outputs.
                        !q gives the delta K values of the matrix elements.
      chiq_r(0)=qz
      chiq_i(0)=0.0 !no imaginary part for q=0 (elements diagonal in K must be real)
      
      chiq_r(1)=-(2.0/3.0)**0.5*(qxz)
      chiq_i(1)=-(2.0/3.0)**0.5*(-qyz)
      
      chiq_r(-1)=(2.0/3.0)**0.5*(qxz) ! single sign change in real part
      chiq_i(-1)=(2.0/3.0)**0.5*(qyz)! double sign change in imaginary part
      
      chiq_r(2) =(1.0/6.0)**0.5*(qd)
      chiq_i(2) =(1.0/6.0)**0.5*(2*qxy)
      chiq_r(-2)=chiq_r(2) !no sign change in real part for q=2
      chiq_i(-2)=-chiq_i(2)     
      

      
      !SET this up in dependence of atot!
      
      if ((jselect-sj).lt.0) sj=jselect
      if (f1select.lt.two_I)then ! if F is smaller than I, then there are J states that can not exist
       if (abs(f1select+2*(jselect-sj)).lt.two_I)then !The non existing states are characterized by F+J < I.
       sj=sj-(abs(abs(f1select-2*(jselect-sj))-two_I)/2)
       end if 
      end if
C      if (f1select.eq.3) then
C      write(0,*) jselect, f1select, jselect-sj, jselect+ej
C      end if
C    First adding the elements diagonal in j 

      do j=jselect-sj,jselect+ej
       cm=j-(jselect-sj)!count matrices
       sizej=2*j+1 
       h_sizes(cm+1)=sizej
       sizepre=sum(h_sizes(:cm+1))-sizej
      



CC           Diagonal Elements delta J=0
C
       do kr=-j,j
        do kc=-j,+j
          do q=-2,2    !the 5 cases for q.   
           if ((-kc-q+kr).eq.0) then !triangle condition for tj symbol.
            
            call sixj(f1select,two_I,2*j,4,2*j,two_I,wsj)
            totprod=wsj*((2*j+1)*(2*j+1))**0.5 !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5

            call threej(2*j, 4, 2*j, -2*kc, -2*q , 2*kr, tj)
            totprod=totprod*tj !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5*tj1
            
            tj=tji
            totprod=totprod/tj
            totprod=totprod/4
            t2=2*kc+two_I+f1select+2
            if(mod(t2,2).ne.0) then
             write(0,*) "t2*2 is not even, this can not be. dj0nq1", t2
             stop
            end if
            t2=t2/2
            totprod=totprod*(-1)**t2 !(-1)**t2*sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5*tj1/tj2/4 ! only thing missing at this point is the chis.
            
            ir = sizepre + kr + j + 1
            ic = sizepre + kc + j + 1
            h_2(ir,ic)=totprod*chiq_r(-q) ! real matrix elements, lower triangle and diagonal
            if (q.ne.0) then
             h_2(ic,ir)=totprod*chiq_i(-q) !upper triangle, imaginary.
            end if
           end if !check only those where non zero matrix elements are expected (see also the first wigner 3jsymbol for this condition)
          end do !q
         end do !ic
        end do   !ir 
        
        
        
C           Off-Diagonal Elements delta J=1

       if ((j).le.(jselect+ej-1)) then !always between j and j+1, requires to go one less j in the loop
       do kr=-j-1,j+1
        do kc=-j,+j
          do q=-2,2    !the 5 cases for q.
           if ((-kc-q+kr).eq.0) then !triangle condition for tj symbol.
            call sixj(f1select,two_I,2*j,4,2*(j+1),two_I,wsj)
            totprod=wsj*((2*j+1)*(2*(j+1)+1))**0.5 !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5
            call threej(2*j, 4, 2*(j+1), -2*kc, -2*q , 2*kr, tj) !also here. 
            totprod=totprod*tj
            tj=tji
            totprod=totprod/tj
            totprod=totprod/4
            t2=2*kc+two_I+f1select+2+2*j+2*(j+1) 
            if(mod(t2,2).ne.0) then
             write(0,*) "t2*2 is not even, this can not be. dj2nq1", t2
             stop
            end if
            t2=t2/2
            totprod=totprod*(-1)**t2
            
            ir=sizepre+(2*j+1)+(kr+j+2)
            ic=sizepre+(kc+j+1)
            h_2(ir,ic)=totprod*chiq_r(-q) ! real matrix elements, lower triangle and diagonal !row index shifts by 2j+1 +1 to match with the new K values
            if (q.ne.0) then
             h_2(ic,ir)= totprod*chiq_i(-q) !upper triangle, imaginary.
            end if
            
           end if !check only those where non zero matrix elements are expected (see also the first wigner 3jsymbol for this condition)
          end do !q
         end do !ic
        end do   !ir 
       end if ! condition for delta j=1
       
       
C           Off-Diagonal Elements delta J=2
       
       if ((j).le.(jselect+ej-2)) then !always between j and j+2, requires to go two less j in the loop
       do kr=-j-2,j+2
        do kc=-j,+j
          do q=-2,2    !the 5 cases for q.
           if ((-kc-q+kr).eq.0) then !triangle condition for tj symbol.
            call sixj(f1select,two_I,2*j,4,2*(j+2),two_I,wsj)
            
            totprod=wsj*((2*j+1)*(2*(j+2)+1))**0.5 !sj1*((2f1+1)*(2f1'+1))**0.5*sj2*((2j+1)*(2j+2))**0.5
            
            call threej(2*j, 4, 2*(j+2), -2*kc, -2*q , 2*kr, tj) !also here. 
            

            
            
            totprod=totprod*tj
            tj=tji
            totprod=totprod/tj
            totprod=totprod/4

            t2=2*kc+two_I+f1select+2+2*j+2*(j+2)  !stays same as for diagonal elements.
            if(mod(t2,2).ne.0) then
             write(0,*) "t2*2 is not even, this can not be. dj2nq1", t2
             stop
            end if

            t2=t2/2
            totprod=totprod*(-1)**t2
            
            ir=sizepre+(2*j+1)+(2*(j+1)+1)+(kr+j+3)
            ic=sizepre+(kc+j+1)

            h_2(ir,ic)=totprod*chiq_r(-q) ! real matrix elements, lower triangle and diagonal !row index shifts by 2j+1 +1 to match with the new K values
C            if ((f1select.eq.3).and.(q.eq.2)) then
C             write(0,*) j,f1select,two_I,q,kr,kc
C             write(0,*) totprod*chiq_r(-q)
C             stop
C            end if
            if (q.ne.0) then
             h_2(ic,ir)= totprod*chiq_i(-q) !upper triangle, imaginary.
            end if
            
           end if !check only those where non zero matrix elements are expected (see also the first wigner 3jsymbol for this condition)
          end do !q
         end do !ic
        end do   !ir 
       end if ! condition for delta j=2
       
       
      end do !j
     
      !!!
      
      return 
      end 
C ---------------------------------------------------------------------
       subroutine rotat1(d,beta,j)
       implicit none
       include 'iam.fi'
       real*8  d(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2)
       real*8  da(-DIMJ-1:DIMJ+1,-DIMJ-1:DIMJ+1,1:2)
       real*8 d1(-2:2,-2:2, 1:2)
       real*8 cg(-1:1,-DIMJ-1:DIMJ+1)
       real*8 beta,djj1
       real*8 dnew,ddnew,cb,sb,dcb,dsb
       integer k,j,k1,k2,m1,m2

       if (j.eq.0) then
         d(0,0,1)=1.0d0
         d(0,0,2)=0.0d0
         return
       endif
C
C       D1(beta) Matrix initialisieren
C
       if(beta.ne.0.0d0) then
         cb=0.5d0 * dcos(beta)
         sb=0.5d0 * dsqrt(2.0d0) * dsin(beta)
         d1(-1,-1,1)=0.5d0 + cb
         d1(-1, 1,1)=0.5d0 - cb
         d1(-1, 0,1)=sb
         d1( 0, 0,1)=2.d0 * cb
         d1( 1, 1,1)= d1(-1,-1,1)
         d1( 0,-1,1)=-d1(-1, 0,1)
         d1( 1,-1,1)= d1(-1, 1,1)
         d1( 0, 1,1)= d1(-1, 0,1)
         d1( 1, 0,1)= d1( 0,-1,1)         
       else
         d1(-1,-1,1)=1.0d0
         d1( 0, 0,1)=1.0d0
         d1( 1, 1,1)=1.0d0
         d1(-1, 1,1)=0.0d0
         d1( 1,-1,1)=0.0d0
         d1( 0,-1,1)=0.0d0
         d1(-1, 0,1)=0.0d0
         d1( 0, 1,1)=0.0d0
         d1( 1, 0,1)=0.0d0
       end if
       dcb=-0.5d0 * dsin(beta)
       dsb=0.5d0 * dsqrt(2.0d0) * dcos(beta)
       d1(-1,-1,2)=dcb
       d1(-1, 1,2)=-dcb
       d1(-1, 0,2)=dsb
       d1( 0, 0,2)=2.d0 * dcb
       d1( 1, 1,2)= d1(-1,-1,2)
       d1( 0,-1,2)=-d1(-1, 0,2)
       d1( 1,-1,2)= d1(-1, 1,2)
       d1( 0, 1,2)= d1(-1, 0,2)
       d1( 1, 0,2)= d1( 0,-1,2)
C
       if ((j.eq.1)) then
         do k1=-1,1
           do k2=-1,1
             d(k1,k2,1)=d1(k1,k2,1)
             d(k1,k2,2)=d1(k1,k2,2)
           end do
         end do
       end if
       if (j.eq.1) return

       do k1=-j,j
         do k2=-j,j
           da(k1,k2,1)=d(k1,k2,1)
           da(k1,k2,2)=d(k1,k2,2)
         end do
       end do
       do k1=1,2
         do k2=-j-1,j+1
           da( j+1,k2,k1)=0.0d0
           da(-j-1,k2,k1)=0.0d0
           da(k2, j+1,k1)=0.0d0
           da(k2,-j-1,k1)=0.0d0
         end do
       end do

C
C   Clebsch Gordan Koeffizienten
C
       do k=-j,j
         djj1=dble(J*(2*J-1))
C       <J-1  M   1  0|JM>
         cg(0,k)=dsqrt(dble((j-k) * (j+k))/djj1)
C       <J-1  M+1 1 -1|JM>
         cg(-1,k)=dsqrt(dble((j-k-1)*(j-k))/(2.0d0*djj1))
C       <J-1  M-1 1 +1|JM>
         cg( 1,k)=dsqrt(dble((j+k-1)*(j+k))/(2.0d0*djj1))
C         write(*,1003) j,k,cg(-1,k), cg(0,k), cg(1,k), djj1
 1003  format(' J  K ',2I4,'   m+1  m  m-1',4F10.4)
       end do
C
C       do k1=0,j
C         do k2=k1,j
       do k1=-j,j
         do k2=-j,j
           dnew=0.0d0
           ddnew=0.0d0
           do m1=-1,1
             do m2=-1,1
             dnew=dnew
     $            + cg(m1,k1)*d1(m1,m2,1)*da(k1-m1,k2-m2,1)*cg(m2,k2)
             ddnew=ddnew
     $            + cg(m1,k1)*d1(m1,m2,1)*da(k1-m1,k2-m2,2)*cg(m2,k2)
     $            + cg(m1,k1)*d1(m1,m2,2)*da(k1-m1,k2-m2,1)*cg(m2,k2)
             end do
           end do
           d(k1,k2,1)=dnew
C          da(-k2,-k1)=dnew
           d(k1,k2,2)=ddnew
C          dda(-k2,-k1)=ddnew
C          if (myand((k2-k1),1).eq.1) then
C            da(k2,k1)= dnew
C            dda(k2,k1)= ddnew
C          else
C            da(k2,k1)=-dnew
C            dda(k2,k1)=-ddnew
C          endif
C          da(-k1,-k2) = da(k2,k1)
C          dda(-k1,-k2) = dda(k2,k1)
         end do
       end do
C
       return
       end
C---------------------
       subroutine argsort(Array,Indices,DIM_A,D) !herbers2024
       implicit none
C      Simple routine that takes array, sorts it decending, and returns the flipped indixes positions, starting from 1 to D
       integer DIM_A
       real*8 Array(DIM_A)
       real*8 T,TI
       integer Indices(DIM_A)
       integer D
       integer i,j
       Indices=0
       DO i=1,D
        Indices(i)=i
       end do
       DO i=1,D-1
          DO j=i+1,D
             IF (Array(i) .LT. Array(j)) THEN
               T = Array(i)  
               TI = Indices(i)
               Array(i)=Array(j) 
               Indices(i) = Indices(j)
               Array(j) = T  
               Indices(j) = TI
             END IF
          END DO
       END DO
 
       return
       end
       
C---------------------
       subroutine argsortAC(Array,Indices,DIM_A,D) !herbers2024
       implicit none
C      Simple routine that takes array, sorts it Acending, and returns the flipped indixes positions, starting from 1 to D
       integer DIM_A
       real*8 Array(DIM_A)
       real*8 T,TI
       integer Indices(DIM_A)
       integer D
       integer i,j
       Indices=0
       DO i=1,D
        Indices(i)=i
       end do
       DO i=1,D-1
          DO j=i+1,D
             IF (Array(i) .GT. Array(j)) THEN
               T = Array(i)  
               TI = Indices(i)
               Array(i)=Array(j) 
               Indices(i) = Indices(j)
               Array(j) = T  
               Indices(j) = TI
             END IF
          END DO
       END DO
 
       return
       end
       
C---------------------
       subroutine twosort(Array,Indices,DIM_A,D) !herbers2024
       implicit none
C      Simple routine that takes array, sorts it Array And Indicies, order of Array values Decending.
       integer DIM_A
       real*8 Array(DIM_A)
       real*8 T,TI
       integer Indices(DIM_A)
       integer D
       integer i,j
       DO i=1,D-1
          DO j=i+1,D
             IF (Array(i) .LT. Array(j)) THEN
               T = Array(i)  
               TI = Indices(i)
               Array(i)=Array(j) 
               Indices(i) = Indices(j)
               Array(j) = T  
               Indices(j) = TI
             END IF
          END DO
       END DO
C
       return
       end
       