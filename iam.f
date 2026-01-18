      program xiam
C     Fitting of Rotational Constants and Centrifugal Distortion Constants 
C     and Internal-rotation parameters of a 3 fold top
C     Syntax of Input-File:
C
C     Phase Convention: <JK|Px|JK+1>=0.5*dsqrt(J*(J+1)-K*(K+1)) 
C  
      implicit none
   
      
      
      include 'iam.fi'

C     arrays for fit-proc
C     DIMFIT: Dimension for paramaters to be fitted simultaneouly
      real*8 covar(DIMFIT,DIMFIT), alpha(DIMFIT,DIMFIT),
     $     evec(DIMFIT,DIMFIT), w(DIMFIT),
     $     beta(DIMFIT), freed(DIMFIT)
C     One fit paramater is a linearcombination of the parameters a 
C     PArameter Linear Combinations as fit parameters 
C     DIMPLC maximum number parameters in one fitparameters
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
C     dfit controls the differences quotient calculation
      integer dfit(DIMFIT)

C     a is the parameter array for the Hamiltonian
      real*8  a(DIMPAR,DIMVB), anew(DIMPAR,DIMVB), da(DIMPAR,DIMVB)
C     ifit controls if analytic derivatives are to be calculated
      integer ifit(DIMPAR,DIMVB)
C     real*8  ab(DIMPAR)

      real*8  chi2(DIMVB,DIMVV,2),wght(DIMVB,DIMVV,2)
      integer ndat(DIMVB,DIMVV,2)
      real*8  sigsq,sigsqold,stepw,cmax
      integer npar,i,izyk,j,iconv,iq,nsvfit,fstat,ift,nfit,itop,deriv
      integer myints,oib,maxf,maxb,maxg,ifunpr,ming
      integer maxf1 !h2024
      integer ij,ib,iv,is,jc,ic,i2f,isf,i2f1
      integer endij
      integer itwof,itwof1
      real*8  maxchn
      character*15 astr
      character*12 dstr
      character*30 binfname
      logical sortend
      real*8  pi,indeg,inkj,inkc,incm,fold,condno
      real*8 tj
      integer two_j1,two_j2,two_j3
      integer two_m1,two_m2,two_m3
      integer myand,myor
      external myand,myor
      integer sig_stat
      common/sig_com/sig_stat

      include 'iamdata.fi'
      call mysignal()
      
      pi=dacos(-1.0d0)
      inkj=3.9903132D-04
      inkc=3.9903132D-04/4.186903D0
      indeg=180.0d0/pi
      incm=1.0d0/29.9792458d0! 
      write(*,'(A,A)')
     $     ' Rotational, Centrifugal Distortion,'
     $     ,' Internal Rotation Calculation (V2.5e)' 
      write(*,'(23X,A)')
     $     'Holger Hartwig 08-Nov-96 (hartwig@phc.uni-kiel.de)'
      write(*,'(/,2A,/)') ' Please cite:',
     $     ' H.Hartwig and H.Dreizler, Z.Naturforsch, 51a (1996) 923.'
C      write(*,'(A,$)') ' Calculation date and time: '
      call mydate()
      write(*,*)
      write(*,'(A,A)') 'Modified Version: XIAM-2NQ v0.32 -' 
     $                  ,'By Sven Herbers 18-January-2026' 
      write(*,*) 'sven_herbers@web.de'
      write(*,*) 'Cite: J. Chem. Phys., 2025, '
     $         ,'162, 234304, DOI: 10.1063/5.0267651 ' 
      write(*,*)
      write(*,'(A)') ' Type help now for the list of parameters : '
      call parinp(a,palc,pali,ifit,dfit,npar,nfit)
      call mysignal()
      
C     the array todo is the list of (good) quantum no.s which
C     have to be calculated to assign all transitions
C     the good quantum no.s are j, f, b and sigma, no symmetry adapted k used ! herbers 2024, added f1 to the list
      size(S_MAXK)=0
      ncalc=0
      maxf=0
      maxf1=0 !h2024
      maxg=0
      maxb=0
      ming=100

C     qlin is the list of transitions
      do ic=1, DIMTDO
        todo(ic,Q_STAT)=0
      end do
          
C     copy qlin into todo
      do i=1, ctlint(C_NDATA)
        do iq=Q_UP, Q_LO
          if (maxg.lt.qlin(i,Q_S,iq)) maxg=qlin(i,Q_S,iq)
          if (ming.gt.qlin(i,Q_S,iq)) ming=qlin(i,Q_S,iq)
          if (maxb.lt.qlin(i,Q_B,iq)) maxb=qlin(i,Q_B,iq)
          if (maxf.lt.qlin(i,Q_F,iq)) maxf=qlin(i,Q_F,iq)
          if (maxf1.lt.qlin(i,Q_F1,iq)) maxf1=qlin(i,Q_F1,iq) !h2024
c          vb(qlin(i,Q_B,iq))=qlin(i,Q_V1,iq)
          if (size(S_MAXK).lt.qlin(i,Q_J,iq))
     $         size(S_MAXK)=qlin(i,Q_J,iq) 
          do ic=1, ncalc
            if ((qlin(i,Q_J,iq).eq.todo(ic,Q_J))
     $         .and.(qlin(i,Q_S,iq).eq.todo(ic,Q_S))
     $         .and.(qlin(i,Q_B,iq).eq.todo(ic,Q_B))
     $         .and.(qlin(i,Q_F,iq).eq.todo(ic,Q_F)) 
     $         .and.(qlin(i,Q_F1,iq).eq.todo(ic,Q_F1))) then !h2024
              goto 50
            end if
          end do
          ncalc=ncalc+1
          if (ncalc.gt.DIMTDO) stop 'ERROR: todo >DIMTDO'
          todo(ncalc,Q_J)=qlin(i,Q_J,iq)
          todo(ncalc,Q_S)=qlin(i,Q_S,iq)
          todo(ncalc,Q_B)=qlin(i,Q_B,iq)
          todo(ncalc,Q_F)=qlin(i,Q_F,iq)
          todo(ncalc,Q_F1)=qlin(i,Q_F1,iq)
          ic=ncalc
 50       continue
          if ((dln(i,LN_ERR).ne.NOFIT).or.(ctlint(C_DFRQ).gt.0)) 
     $         todo(ic,Q_STAT)=myor(todo(ic,Q_STAT),2) 
        end do
      end do
      if (((ctlint(C_SPIN).ne.0).or.(ctlint(C_SPIN2).ne.0))
     $  .and.(ctlint(C_DW).eq.3))    !Herbers2024
     $              size(S_MAXK)=size(S_MAXK)+ctlint(C_SPIN)!Herbers2024
     $                           +ctlint(C_SPIN2)!h2024
C     write(0,*) size(S_MAXK),ctlint(C_SPIN),ctlint(C_SPIN2)
C     for intensities all j,b,f,gam are needed
      if (ctlint(C_NTOP).le.0) then
      if (ctlint(C_INTS).eq.2) ctlint(C_INTS)=3
      end if
      if (ctlint(C_INTS).ge.2) then
        write(*,'(/,A)') ' Intensity Calculation Mode '  
        if (ctlint(C_NZYK).ne.1) then
          ctlint(C_NZYK)=1
          write(*,'(A)') ' Number of iteration cylcles reset to one !'  
        end if
c        if (ctlint(C_NTOP).gt.0) ming=min(1,ming)
        do ib=1, size(S_NB) ! maxb
          endij=size(S_MAXK)
          if (((ctlint(C_SPIN).ne.0).or.(ctlint(C_SPIN2).ne.0))
     $     .and.(ctlint(C_DW).eq.3)) ! maxK changes its meaning a bit throughout here
     $     endij=endij-ctlint(C_SPIN)-ctlint(C_SPIN2)                        ! 
          do ij=0,endij ! used to end at maxK, but now I have to make a DWctrl=3 exception
C     because nq-hfs intensities are not implemented: ignore f  ! Herbers 2024: now they are implemented, but the maxK to be used changes then.
C            do i2f=
            i2f=-1
            i2f1=-1
            isf=0
            if (ctlint(C_INTS).eq.2) isf=1
              do is=isf, size(S_G) 
                do ic=1, ncalc
                  if (      (ij.eq.todo(ic,Q_J))
     $                 .and.(is.eq.todo(ic,Q_S))
     $                 .and.(ib.eq.todo(ic,Q_B))
     $                 .and.(i2f.eq.todo(ic,Q_F)) 
     $                 .and.(i2f1.eq.todo(ic,Q_F1)))goto 51
                end do
                ncalc=ncalc+1
                if (ncalc.gt.DIMTDO) stop 'ERROR: todo >DIMTDO'
                todo(ncalc,Q_J)=ij
                todo(ncalc,Q_S)=is
                todo(ncalc,Q_B)=ib
                todo(ncalc,Q_F)=i2f
                todo(ncalc,Q_F1)=i2f1
 51             continue
              end do
C            end do
          end do
        end do
      end if
C     sort b values in first place, followed by ascending J values
      do iq=1, ncalc
        sortend=.true.
        do i=2, ncalc
          if (todo(i-1,Q_B).gt.todo(i,Q_B)) then
            call swptdo(i,i-1)
            sortend=.false.
          end if
        end do
        if (sortend) goto 13
      end do
      
 13   continue
      do iq=1, ncalc
        sortend=.true.
        do i=2, ncalc
          if ((todo(i-1,Q_J).gt.todo(i,Q_J))
     $         .and.(todo(i-1,Q_B).eq.todo(i,Q_B))) then
            call swptdo(i,i-1)
            sortend=.false.
          end if
        end do
        if (sortend) goto 12
      end do
12    continue
      do iq=1, ncalc
        sortend=.true.
        do i=2, ncalc
          if ((todo(i-1,Q_S).gt.todo(i,Q_S))
     $         .and.(todo(i-1,Q_J).eq.todo(i,Q_J))
     $         .and.(todo(i-1,Q_B).eq.todo(i,Q_B))) then
            call swptdo(i,i-1)
            sortend=.false.
          end if
        end do
        if (sortend) goto 11
      end do
 11   continue
      if (((ctlint(C_SPIN).ne.0).or.(ctlint(C_SPIN2).ne.0))
     $    .and.(ctlint(C_DW).ne.0)) then    !Herbers2024
      write(*,'(/,2X,A,I4)') '\\  Maximal K = J ='               !Herbers2024
     $       ,size(S_MAXK)-ctlint(C_SPIN)-ctlint(C_SPIN2)        !Herbers2024 h2024
      write(*,'(/,2X,A,I4)') '\\  Used maxJ in NQC matrix ='    !Herbers2024
     $       ,size(S_MAXK)                                       !Herbers2024
      else                                                       !Herbers2024
      write(*,'(/,2X,A,I4)') '\\  Maximal K = J =', size(S_MAXK) !Herbers2024
      end if                                                     !Herbers2024
      if (myand(ctlint(C_PRI),AP_ST).ne.0) then
        write(*,'(/,A)') ' ToDo: list of quantum no.'
        write(*,'(3X,A)') '   J  Sym B  F  Stat'
        do i=1,ncalc
          write(*,'(3X,5I4)') todo(i,Q_J),todo(i,Q_S)
     $         ,todo(i,Q_B),todo(i,Q_F),todo(i,Q_STAT)
        end do
      end if
C --------------------------------------------------------------
      oib =-1
      do ic=1, ncalc
        ib =todo(ic,Q_B)
        if (ib.ne.oib) then
          oib=ib 
          ctlint(C_WOODS)=ctlnb(CB_WDS,ib)
          ctlint(C_ADJF) =ctlnb(CB_ADJ,ib)
          write(*,'(2X,2(A,I3))') '\\ B=',ib,'    adj=',ctlint(C_ADJF)
          if (ctlint(C_NTOP).le.0) goto 33
          if (myand(ctlint(C_WOODS),1).ne.0) write(*,'(2X,A)')
     $         '\\ (1) calculate torsional integrals ' 
          if (myand(ctlint(C_WOODS),2).ne.0) write(*,'(2X,A)')
     $         '\\ (2) use scaled torsional integrals ' 
          if (myand(ctlint(C_WOODS),4).ne.0) write(*,'(2X,A)')
     $         '\\ (4) use torsional integrals in the rotation matrix '
          if (myand(ctlint(C_WOODS),8).ne.0) write(*,'(2X,A)')
     $         '\\ (8) use torsional integrals in H_ir for other top(s)' 
          if (myand(ctlint(C_WOODS),16).ne.0) write(*,'(2X,A)')
     $         '\\ (16)use torsional integrals in twotop terms'
          if (myand(ctlint(C_WOODS),32).ne.0) write(*,'(2X,A)')
     $         '\\ (32)use torsional integrals in rigid rotor H_rr ' 
          if (myand(ctlint(C_WOODS),64).ne.0) write(*,'(2X,A)')
     $         '\\ (64)multiply torsional integrals like Demaison'
 33       continue  
        end if
      end do

      do i=1, DIMPAR
        do ib=1, size(S_NB)
          anew(i,ib)=a(i,ib)
        end do
      end do
      sigsqold=1.0d30
      iconv=1
      deriv=1
      fstat=0
      stepw=1.0d0
      myints=ctlint(C_INTS)
      ctlint(C_INTS)=0
      do izyk=1, ctlint(C_NZYK)
        ifunpr=0
        if (ctlint(C_EVAL).ne.0) write(20,'(/,A,i5)') ' Iteration',izyk
        if (izyk.ne.1) write(*,'(/,A,/,A,i3)')
     $       ' -------------------------------------------------------',
     $       ' Iteration :', izyk-1
        xcy=0
        if (myand(ctlint(C_XPR),XP_EC).ne.0) xcy=1
        if ((izyk.eq.1).and.(myand(ctlint(C_XPR),XP_FI).ne.0)) xcy=1
        if ((iconv.ge.2).and.(myand(ctlint(C_XPR),XP_LA).ne.0)) xcy=1

        if ((myand(ctlint(C_PRI),AP_PC).ne.0).and.(xcy.ge.1))
     $       call prpar(a)
                
C --  calc. the spectrum
        if (ctlint(C_NZYK).eq.1) ctlint(C_INTS)=myints
        if (iconv.eq.2) ctlint(C_INTS)=myints
        call calc1(anew,ifit,dfit,npar,nfit,palc,pali,0)
        if (ctlint(C_NZYK).eq.1) goto 40

C --  get sigsq
        fstat=0
        call lmfit(ctlint(C_NDATA),npar,size(S_NB)
     $       ,DIMFIT,DIMPAR,DIMVB,DIMPLC
     $       ,ctlint(C_PRI),nsvfit,nfit,ifit,dfit
     $       ,alpha,covar,evec,beta,w,a,anew,da,freed,palc,pali
     $       ,ctlpar(C_ROFIT),sigsq,sigsqold,ctlpar(C_EPS)
     $       ,stepw,fstat,ctlint(C_ORGER),0
     $       ,ctlpar(C_LMBDA),ctlint(C_FITSC),ctlint(C_SVDER))

C --  test if sigsq was improved    
        if ((sigsq/sigsqold).le.ctlpar(C_CNVG)) iconv=1
        if (((sigsq/sigsqold).gt.ctlpar(C_CNVG)) 
     $       .and.((sigsq/sigsqold).le.(2.0d0-ctlpar(C_CNVG))))
     $       iconv=iconv+1

C --  switch to better derivatives if the fit is not converging
        if ((sigsq/sigsqold).gt.(2.0d0-ctlpar(C_CNVG))) then
          if (iconv.gt.0) iconv=0
          if (iconv.le.0) iconv=iconv-1
        end if
C--   set deriv to 3 at the the first time; deriv is set to 2 later on.
        if ((iconv.eq.-3).and.(deriv.eq.1)) then 
          deriv=3
          write(*,*)'Switching to better derivatives (takes more time)'
        end if 

        if ((sigsq.lt.sigsqold).and.(myand(ctlint(C_XPR),XP_CC).ne.0))
     $       xcy=1

        write(*,'(A,D12.6,A,F8.6,A,I2)')
     $       ' Sigma:',sqrt(sigsq)
     $       ,'  Sigma/OldSigma:',sqrt(sigsq/sigsqold)
     $       ,'  conv:',iconv

C --  print the parameters if anew is better than before
        if ((izyk.ne.1).and.(xcy.ge.1)) then
          write(*,'(/,15X,A,5X,A)') 'Parameters','Change'
          call pra(anew,da,ifit,myand(ctlint(C_PRI),AP_PL),1,0)
        end if

C --  print the parameter with the max change
        if (izyk.ne.1) then
          maxchn=0.0
          ic=0
          do i=1, npar
            do ib=1, size(S_NB)
              if ((anew(i,ib).ne.0.0).and.(ifit(i,ib).ne.0)) then
                if (abs(da(i,ib)/anew(i,ib)).gt.maxchn) then
                  maxchn=abs(da(i,ib)/anew(i,ib))
                  ic=i
                  ij=ib
                end if
              end if
            end do
          end do
          if (ic.gt.0) then
            call pra_f(anew(ic,ij),da(ic,ij),astr,dstr)
            write(*,'(1X,2A,I1,4A,F8.3,A)') parstr(ic),'(',ij,')'
     $           ,astr,'  ',dstr,maxchn*100.0,'% Max. Change'
          end if
        end if

C --  print the transitions
        if (((myand(ctlint(C_PRI),AP_TL).ne.0).and.(xcy.ge.1)).or.
     $       ((myand(ctlint(C_PRI),AP_TF).ne.0).and.(izyk.eq.1))) then
          write(*,*)
          ifunpr=1
          call funpr(ctlpar(C_ROFIT),ndat,chi2,wght)
        end if 
        
C --  if sigsq is better then oldsigsq calc the dqu derivatives
        if ((sigsq.lt.sigsqold).or.(ctlint(C_DFRQ).gt.1)
     $       .or.(deriv.eq.3)) then
          call calc1(anew,ifit,dfit,npar,nfit,palc,pali,deriv)
          if (ctlint(C_DFRQ).ne.0) then
            write(21,'(/,A,i5)') ' Iteration',izyk
            call prderiv(dfit,nfit)
          end if
        end if
        if (deriv.eq.3) deriv=2
        
C --  if convergence is achived: fstat=-1 to calc. the errors
        fstat=1
        if (iconv.gt.2) fstat=-1  
        if (izyk.eq.ctlint(C_NZYK)) fstat=-1001
        if (sig_stat.eq.1) fstat=-1002
        if (ctlpar(C_LMBDA).gt.1.0d3) fstat=-1003
        call lmfit(ctlint(C_NDATA),npar,size(S_NB)
     $       ,DIMFIT,DIMPAR,DIMVB,DIMPLC
     $       ,ctlint(C_PRI)
     $       ,nsvfit,nfit,ifit,dfit
     $       ,alpha,covar,evec,beta,w,a,anew,da,freed,palc,pali
     $       ,ctlpar(C_ROFIT),sigsq,sigsqold,ctlpar(C_EPS)
     $       ,stepw,fstat,ctlint(C_ORGER),0
     $       ,ctlpar(C_LMBDA),ctlint(C_FITSC),ctlint(C_SVDER))

        if (fstat.ge.0) then
          do ib=1, size(S_NB)
            call adjusta(anew(1,ib),npar,ctlnb(CB_ADJ,ib))
          end do
        end if
        condno=0.0d0
        do i = 1,nfit
          if (w(i).gt.0.0d0) then
            condno=w(nfit)/w(i)
            goto 60
          end if
        end do
 60     continue
        write(*,'(/,A,I3,A,F6.4,A,D9.3,A,D9.3)')
     $       ' indep.par:',nsvfit,
     $       ' stepw:',stepw,
     $       ' lambda:',ctlpar(C_LMBDA),
     $       ' cond.no:',condno

        if (((myand(ctlint(C_PRI),AP_SV).ne.0).and.(xcy.ge.1))
     $       .and.(fstat.gt.-1)) then
          write(*,'(/,A,A,/)') ' Eigenvalues and '
     $         ,'Eigenvector Matrix of SVD-FIT '
          do i = 1, nfit
            write(*,'(1D12.6,2X,30F6.3)') w(i),
     $           (evec(j,i),j=1,nfit)
          end do
        end if        
        if (sigsq.lt.sigsqold) sigsqold=sigsq
C --  stop here if convergence        
        if (fstat.lt.0) goto 20
      end do
 30   continue
      write(*,'(/,A)') ' No Convergence, fit aborted '
 20   continue
      
      if ((iconv.eq.0).or.(ctlint(C_INTS).ne.myints)
     $     .or.(sigsq.gt.sigsqold)) then
        ctlint(C_INTS)=myints
        ifunpr=0
        write(*,'(/,A)') ' Recalculation of the spectrum'
        call calc2(a,ifit,dfit,npar,nfit,palc,pali,0)
      end if
      izyk=izyk-1
      
      write(*,'(A,/,A,I3,/)')
     $     '##########################################################',
     $     '  End at Cycle ',izyk 

 40   continue
      if (ifunpr.eq.0) then
        call funpr(ctlpar(C_ROFIT),ndat,chi2,wght)
      end if
      if (nfit.gt.0) then
        write(*,'(/,A)')' RMS deviations (MHz), B and V sorted'
        do j=1, 2
         if (j.eq.1) write(*,'(2A3,A4,A)') 'B','V','n',' splittings MHz'
         if (j.eq.2) write(*,'(2A3,A4,A)') 'B','V','n',' abs. freq. MHz'
          do ib=1, DIMVB
            do iv=1, DIMVV
              if (ndat(ib,iv,j).ne.0) then
                write(*,'(2I3,I4,2F18.6)') ib,iv,ndat(ib,iv,j)
     $               ,sqrt(chi2(ib,iv,j)/wght(ib,iv,j))*1.D3
     $               ,sqrt(chi2(ib,iv,j)/wght(ib,iv,j))*1.D3
     $               *(ndat(ib,iv,j)+nfit)/ndat(ib,iv,j)
              end if
            end do
          end do
        end do
      end if
      if (ctlint(C_PRINT).eq.6) then
      if (nfit.gt.0) then
        write(*,'(/,A)')' RMS deviations (cm-1), B and V sorted'
        do j=1, 2
        if (j.eq.1) write(*,'(2A3,A4,A)') 'B','V','n',' splittings cm-1'
        if (j.eq.2) write(*,'(2A3,A4,A)') 'B','V','n',' abs. freq. cm-1'
          do ib=1, DIMVB
            do iv=1, DIMVV
              if (ndat(ib,iv,j).ne.0) then
                write(*,'(2I3,I4,2F18.6)') ib,iv,ndat(ib,iv,j)
     $               ,sqrt(chi2(ib,iv,j)/wght(ib,iv,j))*incm
     $               ,sqrt(chi2(ib,iv,j)/wght(ib,iv,j))*incm
     $               *(ndat(ib,iv,j)+nfit)/ndat(ib,iv,j)
              end if
            end do
          end do
        end do
      end if
	  end if      
      if (ctlint(C_NZYK).eq.1) goto 100
      if (nsvfit.ne.nfit) then
        write(*,'(/,A,I3,A,I3,A)')' **** Warning **** : only',nsvfit
     $       ,' of',nfit,' Parameters are linear independent!'
        if (ctlint(C_SVDER).ne.0) write(*,'(A)')
     $       ' Errors may be wrong!  set ''svderr'' to 0 recommended'
      end if
      write(*,'(/,12X,A)') 'Parameters  and  Errors'
      call pra(a,da,ifit,min(myand(ctlint(C_PRI),AP_PL)+1,5),1,0)

      write(*,'(/,A,F15.6,A,/)') ' Standard Deviation'
     $     ,dsqrt(sigsqold)*1000.0d0,' MHz'

      do ib=1, size(S_NB)
        write(*,'(A,I3)')'------------------------------------- B =',ib
        call adjusta(a(1,ib),npar,ctlnb(CB_ADJ,ib))
        call prrrp(a(1,ib),da(1,ib),covar,palc,pali,nfit,ib)
        if (ctlint(C_NTOP).gt.0) then
          call prirp(a(1,ib),da(1,ib),indeg)
          call prpot(a(1,ib),da(1,ib),inkj,inkc,incm)
          if (ctlint(C_NTOP).gt.1) then
            fold=a(P1_F,ib)
            call adjusta(a(1,ib),npar,myor(ctlnb(CB_ADJ,ib),1))
            if (a(P1_F,ib).ne.fold) write(*,'(A,/,20X,A)')
     $           '**** Warning **** : F values not consistent '
     $           ,'Use adjf = 1 to keep F values right.'
            do itop=1, ctlint(C_NTOP)
              ift=DIMPIR*(itop-1)
              write(*,'(A,1F18.9)') 'F(calc)  ',a(P1_F+ift,ib)
            end do
            if (ctlint(C_NTOP).gt.2) write(*,'(A,1F18.9)')
     $           'F12(calc)',a(P_FF,ib)
          end if
        end if
      end do

      write(*,'(/,A)') ' Errors of fitted linear combinations' 
      write(*,'(5F15.9)')
     $       (anew(mod(i-1,DIMPAR)+1,int((i-1)/DIMPAR)+1),i=1, nfit)
      write(*,'(/,A)') 
     $     ' Correlation Matrix of fitted linear combinations ' 
      cmax=0.0
      do i=1, nfit
        do j=1, i-1
          if (abs(covar(j,i)).gt.cmax) then
            ic=i
            jc=j
            cmax=abs(covar(j,i))
          end if
        end do
        write(*,'(1X,A6,$)') parstr(pali(i,1,1))
        write(*,'(20F7.3)') (covar(j,i), j=1,i-1),1.0d0
      end do
      write(*,'(A,I3,A,I3,A,F7.4,A)') ' strongest correlation between'
     $     ,ic,' and',jc,' (',covar(jc,ic),')'
      write(*,'(/,A)') ' Freedom Cofreedom Matrix of linear comb.' 
      cmax=2.0
      do i=1, nfit
        do j=1, i-1
          if (covar(i,j).lt.cmax) then
            ic=i
            jc=j
            cmax=covar(i,j)
          end if
        end do
        write(*,'(1X,A6,$)') parstr(pali(i,1,1))
        write(*,'(20F7.3)') (covar(i,j), j=1,i)
      end do
      write(*,'(A,I3,A,I3,A,F7.4,A)') ' minimum cofreedom between'    
     $     ,ic,' and',jc,' (',cmax,')'

      write(*,'(/,A,A,/)') ' Eigenvalues and '
     $     ,'Eigenvector Matrix of SVD-FIT '
      do i = 1, nfit
        write(*,'(1D12.6,2X,30F6.3)') w(i),
     $       (evec(j,i),j=1,nfit)
      end do

 100  continue
      write(*,*)
C      if (myand(ctlint(C_PRI),AP_LT).ne.0) call ltxpr(ctlpar(C_ROFIT)) ! Herbers2024 - removed the latex print.
      if (myand(ctlint(C_INTS),16).eq.0) then
         do i=1, ncalc
            call binnam(todo(i,Q_J),todo(i,Q_S),todo(i,Q_F),todo(i,Q_B)
     $           ,todo(i,Q_F1),binfname)
c            write(*,*) binfname
            open(99,file=binfname)
            close(99,status='delete')
         end do
      end if
      
      if (myand(ctlint(C_INTS),16).eq.0) then
      if ((ctlint(C_DW).eq.3).and.
     $   ((ctlint(C_SPIN).ne.0).or.(ctlint(C_SPIN2).ne.0))) then ! when the number of temporary files also includes J+-2, then these have to be also delted
      do j=0, size(S_MAXK) ! But I do not want to figure out the names of these, so I just go over all temporary files again.
      do i=0, int(3**ctlint(C_NTOP)/2+1)
      do itwof1=2*j-ctlint(C_SPIN),2*j+ctlint(C_SPIN),2 
      do itwof=itwof1-ctlint(C_SPIN2),itwof1+ctlint(C_SPIN2),2 
      do ib=1, DIMVB
         call binnam(j,i,itwof,ib
     $           ,itwof1,binfname)
            open(99,file=binfname)
C           write(0,*) binfname
            close(99,status='delete')
      end do
      end do
      end do
      end do
      end do
      end if
      end if
      open(99,file='tori.b1')
      close(99,status='delete')
      stop 
      end

C----------------------------------------------------------------------
      subroutine sig_func(sig_no)
C     called with signal(2)=SIGINT=control C
C     see procedure mysignal() in iamsys.f
      implicit none
      integer sig_stat,sig_no
      common/sig_com/sig_stat
      sig_stat=1
      write(0,*)'Signal no.',sig_no
      write(0,*)'Control-C: premature termiation of xiam, finishing...'
      write(*,*)'Control-C: premature termiation of xiam, finishing...'
      return
      end

C----------------------------------------------------------------------
      subroutine funcs(ix,df,dfda,a,sig,nfit,ifit,dfit,idfrq)
C     interface subroutine between LM Fit and the dnv matrix
      implicit none
      include 'iam.fi'
      integer ix, nfit, idfrq
      integer ifit(DIMPAR,DIMVB),dfit(DIMFIT)
      real*8  df, sig, dfda(DIMFIT),a(DIMPAR,DIMVB)
C     local..
      real*8  dfdao(DIMFIT)
      real*8  dfobs,dfcal,dupda(DIMFIT),dloda(DIMFIT),dfsig,dfor,dfcr
      integer avg,i,ref
      logical firstt
      save firstt
      external myand,myor
      integer myand,myor
      data firstt/.true./

      dfor=0.0d0
      dfcr=0.0d0
      do  i=1, nfit
        dfdao(i)=0.0d0
        dfda(i)=0.0d0
      end do

      dfobs=dln(ix,LN_FREQ)
      dfsig=dln(ix,LN_ERR)

C     average of the data, normaly avg=ix: this has no effect 
      avg=qlin(ix,Q_AVG,Q_UP)
      if (avg.eq.0) avg = ix
      dfcal=((dnv(avg,NV_ENG,Q_UP)-dnv(avg,NV_ENG,Q_LO))
     $     +(dnv(ix,NV_ENG,Q_UP)-dnv(ix,NV_ENG,Q_LO)))/2.0d0
      do i=1, nfit
        dupda(i)=(dnv(avg,i+NV_DEF,Q_UP)+dnv(ix,i+NV_DEF,Q_UP))/2.0d0
        dloda(i)=(dnv(avg,i+NV_DEF,Q_LO)+dnv(ix,i+NV_DEF,Q_LO))/2.0d0
      end do
     
C     check if a difference is fitted 
      ref=qlin(ix,Q_REF,Q_UP)
      if ((ref.ne.0).and.(ref.ne.ix)) then
        dfor=dln(ref,LN_FREQ)
        dfcr=dnv(ref,NV_ENG,Q_UP)-dnv(ref,NV_ENG,Q_LO)
        do i=1, nfit
          dfdao(i)=-(dnv(ref,i+NV_DEF,Q_UP)-dnv(ref,i+NV_DEF,Q_LO))
        end do
      end if  
      
      sig=dfsig
      if ((myand(qlin(ix,Q_STAT,Q_LO),1).ne.1).or.
     $     (myand(qlin(ix,Q_STAT,Q_UP),1).ne.1)) sig=NOFIT

      df=dfobs-dfor-dfcal+dfcr
      if ((ref.ne.0).and.(ref.ne.ix).and.(sig.lt.1d3*ctlpar(C_DEFER)))
     $     then
        if (abs(dfor-dfobs-dfcal+dfcr).lt.abs(df)) then 
          if (firstt) write(*,'(/,A,A,/)') '   *Warning*:'
     $      ,' possible assignment error in line(s) marked with ''x'''
          if (firstt) firstt=.false.
          qlin(ix,Q_STAT,Q_UP)=myor(qlin(ix,Q_STAT,Q_UP),8)
        end if
      end if
      do i=1, nfit
        dfda(i)=dfdao(i)+dupda(i)-dloda(i)
      end do
      return
      end
      
C----------------------------------------------------------------------
      subroutine calc1(a,ifit,dfit,npar,nfit,palc,pali,istat)
C     simple calculation of the spectrum if istat .le. 0
C     calculation of derivatives if istat.gt.0
C     if istat.gt.1 better derivatives are used

      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR,DIMVB)
      integer ifit(DIMPAR,DIMVB),dfit(DIMFIT),npar,nfit,istat
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
C     local ..
      integer ifitmp(DIMPAR,DIMVB)
      real*8  ad(DIMPAR,DIMVB)
      real*8  dnvtmp(DIMLIN,Q_UP:Q_LO),dnvsav(DIMLIN,Q_UP:Q_LO)
      integer i,j,pri,ib
      real*8  diffup,difflo,delta,dsum,devar(DIMPAR)
      integer myand
      real*8  myrand
      external myand,myrand
      save devar, dnvsav
      data devar /DIMPAR*1.0/
      if (istat.le.0) then
C     clear the old eigenvalues und derivativs
        do i=1, ctlint(C_NDATA)
          do j=1,DIMDNV
            dnv(i,j,Q_UP)=0.0d0
            dnv(i,j,Q_LO)=0.0d0
          end do
        end do
        xde=0
        if (xcy.ge.1) xde=1                   ! For output control
C     calculate with original values, obtain derivatives if dfit(i).gt.0
        call calc2(a,ifit,dfit,npar,nfit,palc,pali,0)
        
C     save the original frequencies
        do i=1, ctlint(C_NDATA)
          dnvsav(i,Q_UP)=dnv(i,NV_ENG,Q_UP)
          dnvsav(i,Q_LO)=dnv(i,NV_ENG,Q_LO)
        end do

      else
C     calculate derivatives with difference quotient if dfit(i).lt.0
C     but do not calc. analtic deriv. here
        do ib=1, DIMVB
          do i=1, DIMPAR
            ifitmp(i,ib)=0
          end do
        end do
        xde=0
        if ((myand(ctlint(C_PRI),XP_DE).ne.0).and.(xcy.ge.1)) xde=1
        pri=ctlint(C_PRI)
        ctlint(C_PRI)=myand(pri,not(AP_EH))
        ctlint(C_PRI)=myand(pri,not(AP_MH))

C     calculate the ad = parameters + delta
        do i=1, nfit
          if (dfit(i).lt.0) then
            devar(i)=-devar(i)
            dsum=0.0d0
            do j=1, pali(i,0,1)
              dsum=dsum+dabs(a(pali(i,j,1),pali(i,j,2))*palc(i,j))
            end do
            delta=(palc(i,0)/100.0d0)*dsum*devar(i)/dble(pali(i,0,1))
            do ib=1, size(S_NB)
              do j=1, npar
                ad(j,ib)=a(j,ib)
              end do
            end do
            do j=1, pali(i,0,1)
              ad(pali(i,j,1),pali(i,j,2))
     $             =ad(pali(i,j,1),pali(i,j,2))+delta*palc(i,j)
            end do
            call calc2(ad,ifitmp,dfit,npar,nfit,palc,pali,i)

C     do an opposite variation for precise derivatives 
            if ((dfit(i).le.-2).or.(istat.gt.1)) then
              do ib=1, size(S_NB)
                do j=1, npar
                  ad(j,ib)=a(j,ib)
                end do
              end do
              do j=1, pali(i,0,1)
                ad(pali(i,j,1),pali(i,j,2))
     $               =ad(pali(i,j,1),pali(i,j,2))-delta*palc(i,j)
              end do
              do j=1, ctlint(C_NDATA)
                dnvtmp(j,Q_UP)=dnv(j,NV_ENG,Q_UP)
                dnvtmp(j,Q_LO)=dnv(j,NV_ENG,Q_LO)
              end do
              call calc2(ad,ifitmp,dfit,npar,nfit,palc,pali,i)

C     calculate differential quotient          
              do j=1,ctlint(C_NDATA)
                diffup=dnvtmp(j,Q_UP)-dnv(j,NV_ENG,Q_UP)
                difflo=dnvtmp(j,Q_LO)-dnv(j,NV_ENG,Q_LO)
                dnv(j,i+NV_DEF,Q_UP)=diffup/(2.0d0*delta)
                dnv(j,i+NV_DEF,Q_LO)=difflo/(2.0d0*delta)
              end do
            else
C     calculate differential quotient          
              do j=1,ctlint(C_NDATA)
                diffup=dnv(j,NV_ENG,Q_UP)-dnvsav(j,Q_UP)
                difflo=dnv(j,NV_ENG,Q_LO)-dnvsav(j,Q_LO)
                dnv(j,i+NV_DEF,Q_UP)=diffup/delta
                dnv(j,i+NV_DEF,Q_LO)=difflo/delta
              end do
            end if
          end if
        end do
        ctlint(C_PRI)=pri 
C     restore the original frequencies
        do i=1, ctlint(C_NDATA)
          dnv(i,NV_ENG,Q_UP)=dnvsav(i,Q_UP)
          dnv(i,NV_ENG,Q_LO)=dnvsav(i,Q_LO)
        end do
      end if
      return
      end

C----------------------------------------------------------------------
      subroutine calc2(a,ifit,dfit,npar,nfit,palc,pali,fistat)
C     calculation of the eigenvalues 
C     the evalues are put in the field of dnv(1..ndata,NV_ENG,Q_UP/LO)
C     the deviations DE/DPi in dnv(1..ndata,2-DIMPAR,Q_UP/LO(i))

      implicit none
      include 'iam.fi'
      real*8  a(DIMPAR,DIMVB)
      integer ifit(DIMPAR,DIMVB),dfit(DIMFIT),npar,nfit,fistat
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
      integer fstatus(-1:DIMJ+(DIMQ+1)/2+(DIMQ2+1)/2,0:1,
     $          0:DIMGAM,DIMVB)!
      integer hf,hf1
      integer caseff1
      
C      real*8  hsdq(DIMQ2,DIMQ,DIMTOT,DIMTOT)               
      real*8  evhsdq(DIMQ2,DIMQ+DIMQ2,DIMTOT)                      
      integer jsaved
      integer jcheck
      integer runj
      integer newf
      integer stepf
      integer newf1
      integer stepf1
      
C     work
      real*8  hsdw(DIMDW,DIMTOT,DIMTOT) 
      real*8  evhsdw(DIMDW,DIMTOT)!one extra to keep double well working...
      real*8  h(DIMTOT,DIMTOT)
      real*8  evh   (DIMTOT)!one extra to keep double well working...
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      integer qmv(DIMV),oldj(DIMTOP)
      real*8  ints
      character*30 fmtstr2
      integer j,gam,qf,ib,oib,ic,i,itop,maxi,mini,iv,ntop,is      
      integer qf1
      integer myand,not
      external myand
      save    evalv,ovv,rotm,rott,tori,qmv,oib
      fstatus=-1


      do i=1,ctlint(C_NDATA)
        qlin(i,Q_STAT,Q_UP)=myand(qlin(i,Q_STAT,Q_UP),(not(1)))
        qlin(i,Q_STAT,Q_LO)=myand(qlin(i,Q_STAT,Q_LO),(not(1)))
        qlin(i,Q_STAT,Q_UP)=myand(qlin(i,Q_STAT,Q_UP),(not(8)))
        qlin(i,Q_STAT,Q_LO)=myand(qlin(i,Q_STAT,Q_LO),(not(8)))
        qlin(i,Q_GK,Q_UP)=-1000
        qlin(i,Q_GK,Q_LO)=-1000
      end do
      do i=1, size(S_G)
        gamma(i,0)=0
      end do

      oib =-1
      do ic=1, ncalc
        j  =todo(ic,Q_J)
        gam=todo(ic,Q_S)
        qf =todo(ic,Q_F)
        qf1 =todo(ic,Q_F1) !h2024
        ib =todo(ic,Q_B)
        if ((myand(todo(ic,Q_STAT),2).eq.0).and.
     $       (fistat.ne.0)) goto 99
        ctlint(C_WOODS)=ctlnb(CB_WDS,ib)
        ctlint(C_ADJF) =ctlnb(CB_ADJ,ib)
        size(S_VV)=1
        do iv=1, DIMVV
          if (qvv(iv,1,ib).eq.-1) goto 11
          size(S_VV)=iv
        end do
 11     continue
        if (size(S_VV).gt.DIMVV) stop ' max VV > DIM VV'
        do itop=1, ctlint(C_NTOP)
          size(S_V+itop)=0
        end do
        do itop=1, ctlint(C_NTOP)
          mini=99
          maxi=-1
          do iv=1, size(S_VV)
            if (qvv(iv,itop,ib).lt.mini) mini=qvv(iv,itop,ib) 
            if (qvv(iv,itop,ib).gt.maxi) maxi=qvv(iv,itop,ib) 
          end do
          size(S_V+itop)=maxi-mini+1
          size(S_MINV+itop)=mini
          if (size(S_V+itop).lt.0) stop 'error in calc2'
        end do

        if (ib.ne.oib) then
          oib=ib 
          do itop=1,ctlint(C_NTOP)
            oldj(itop)=0
          end do

          call adjusta(a(1,ib),npar,ctlint(C_ADJF))
          if ((myand(ctlint(C_PRI),AP_PC).ne.0)
     $         .and.(xde.ne.0)) then
            write(*,'(2(A,I3))') 'Parameter for B=',ib
     $           ,'   deriv.=',fistat
            write(*,'(A,3F18.9)') 'BJ  BK    B-  '
     $           ,a(P_BJ,ib),a(P_BK,ib),a(P_BD,ib)
            do itop=1, ctlint(C_NTOP)
              write(*,'(A,3F18.9)') 'rho gamma beta'
     $             ,a(P1_RHO +DIMPIR*(itop-1),ib)
     $             ,a(P1_GAMA+DIMPIR*(itop-1),ib)
     $             ,a(P1_BETA+DIMPIR*(itop-1),ib)
            end do
          end if
C     calc the |m> and |K> part in the rho-system
          call calmk(ib,h,evalv,ovv,rotm,rott,tori
     $         ,a(1,ib),qmv,ifit(1,ib),npar,fistat,0)
C     and save the torsional integrals for the intensities
          if ((ctlint(C_INTS).gt.1).and.(fistat.eq.0))
     $         call wrtori(tori,ib)          
        end if
C     set up the rotation matrix 
        do itop=1,ctlint(C_NTOP)
          call rotate(rotm(-DIMJ,-DIMJ,1,itop)
     $         ,a(P1_BETA+DIMPIR*(itop-1),ib),j,oldj(itop))
        end do
C  !!!!!!!!!!!!!!!!!!  ntop !!!!!!!!!!!!!!!!!!!
        ntop=ctlint(C_NTOP)
        if (gam.eq.0) ctlint(C_NTOP)=0
      if (ctlint(C_DW).eq.1) then ! IF-DW
        call calvjk_d(j,gam,qf,ib,qf1,evalv,ovv,rotm,rott,tori
     $    ,a,qmv,ifit,dfit,palc,pali,npar,fistat
     $    ,hsdw,evhsdw)      
      else ! ELSE-DW
       if ((ctlint(C_DW).eq.3).and.((ctlint(C_SPIN).ne.0))) then !IF-Q
          if (mod(qf,2).eq.1) then  !IF-ODD-F
            hf=(qf+1)/2
          else !ELSE-ODD-F
            hf=(qf)/2
          endif !END-IF-ODD-F
          
          if (mod(qf1,2).eq.1) then  !IF-ODD-F1
            hf1=(qf1+1)/2
          else !ELSE-ODD-F1
            hf1=(qf1)/2
          endif !END-IF-ODD-F1
          caseff1=0
          if (qf1.ge.0) then !IF-QF1
          if (qf.eq.-1) then !IF-QF-NOTUSED
          hf=hf1 ! use f1 in case F is not used (single nucleus.
          caseff1=0
          else !ELSE-QF-NOTUSED
          caseff1=1
          end if !END-IF-QF-NOTUSED
            if (fstatus(hf,caseff1,gam,ib).ne.-1)then !IF-FSTATUS (checks if F matrix was already calculated
C              DO NOTHING!
C              DO NOTHING!
C              DO NOTHING!
            else !ELSE-FSTATUS
                 call calvjk_qdj2f12(j,gam,qf,ib,qf1,evalv,ovv,rotm
     $        ,rott,tori,a,qmv,ifit,dfit,palc,pali,npar,fistat
     $        ,evhsdq)
              fstatus(hf,caseff1,gam,ib)=j
            end if !END-IF-FSTATUS
          else ! ELSE-QF1 !no exception yet implemented. But this is the case for intensities. I will put standard xiam routine
C                 call calvjk_qdj2f12(j,gam,qf,ib,qf1,evalv,ovv,rotm
C     $        ,rott,tori,a,qmv,ifit,dfit,palc,pali,npar,fistat
C     $        ,evhsdq)
          call calvjk(j,gam,qf,ib,qf1,h,evalv,ovv,rotm,rott,tori
     $    ,a(1,ib),qmv,ifit(1,ib),dfit,palc,pali,npar,fistat,evh)
          if ((fistat.eq.0).and.(ctlint(C_INTS).ge.1)) then   !IF-INT                 !Intensity f loop for energy calculation and saving
                                                                               !Intensity f loop
          do stepf=0,ctlint(C_SPIN2)+ctlint(C_SPIN)                            !only need f loop, since f1 and j are produced along the way.                                     
          newf=2*j-ctlint(C_SPIN2)-ctlint(C_SPIN)+2*stepf   
          if (newf.lt.0) then
           newf=-1*newf
          end if
          if (newf.le.2*j) then                              !IF-leJ
                 newf1=newf+ctlint(C_SPIN2)!to not trigger stop conditions I implemented in imav.f
          else                                               !Else-leJ
                 newf1=newf-ctlint(C_SPIN2)!to not trigger stop conditions I implemented in imav.f
          end if                                             !END-IF-leJ
          if (newf1.lt.0) then
           newf1=-1*newf1
          end if
           
           if (newf.ge.0) then !IF-F only calculate F matrices with F bigger equal to 0
                 call calvjk_qdj2f12(j,gam,newf,ib,newf1,evalv,ovv,rotm          !Intensity f loop     ! some of these are redundant.
     $        ,rott,tori,a,qmv,ifit,dfit,palc,pali,npar,fistat                 !Intensity f loop       ! some of these are redundant.
     $        ,evhsdq)                                                         !Intensity f loop
           end if !END-IF-F
          end do                                                               !Intensity f loop
          end if                                            !END-INT
          end if !END-QF1 (end of case that qf1<0, a marker for intensity calculation, but could also be used in linelists to deactivate hyperfine effects.
          
        else !ELSE-Q!Standard XIAM, with no offdiagonal elements in J
          call calvjk(j,gam,qf,ib,qf1,h,evalv,ovv,rotm,rott,tori
     $    ,a(1,ib),qmv,ifit(1,ib),dfit,palc,pali,npar,fistat,evh)
        endif !END-IF-Q
      endif !END-IF-DW

C     calc the intens
C     calc the intensities for this ib of all eigenvalues are done
        if (((ic.eq.ncalc).or.(ib.ne.todo(ic+1,Q_B)))
     $       .and.(fistat.eq.0).and.(ctlint(C_INTS).ge.1)) then
          do i=1,ctlint(C_NDATA)
            if ((qlin(i,Q_B,Q_UP).eq.ib).and.(qlin(i,Q_B,Q_LO).eq.ib))
     $           then
              call intens(i,ints,a(P_MUX,ib),a(P_MUY,ib),a(P_MUZ,ib)
     $             ,tori)
              dln(i,LN_INT)=ints
            end if
          end do
        end if
        ctlint(C_NTOP)=ntop
C  !!!!!!!!!!!!!!!!!!  ntop !!!!!!!!!!!!!!!!!!!
 99     continue
      end do

C     test if all eigenvalues are calculated and assigned
      do i=1,ctlint(C_NDATA)
        if ((myand(qlin(i,Q_STAT,Q_UP),16).eq.16)
     $       .or.(myand(qlin(i,Q_STAT,Q_LO),16).eq.16)
     $       .or.(fistat.eq.0)) then
          if ((myand(qlin(i,Q_STAT,Q_UP),1).ne.1)
     $         .or.(myand(qlin(i,Q_STAT,Q_LO),1).ne.1))
     $         write(*,'(A,I4,A,2I4)')
     $         ' Assign error at line',i,' reason',
     $         qlin(i,Q_STAT,Q_UP),qlin(i,Q_STAT,Q_LO)
        end if
      end do

C     call the routines for intensity calculation
      if ((ctlint(C_INTS).gt.1).and.(fistat.eq.0)) then 
        if (ctlint(C_INTS).eq.2)
     $       call intall(a(P_MUX,ib),a(P_MUY,ib),a(P_MUZ,ib)
     $       ,ctlpar(C_TEMP))
        if (ctlint(C_INTS).eq.3)
     $       call intal2(a(P_MUX,ib),a(P_MUY,ib),a(P_MUZ,ib)
     $       ,ctlpar(C_TEMP))
      end if
      return
      end

C----------------------------------------------------------------------
      subroutine calmk(ib,h,evalv,ovv,rotm,rott,tori
     $     ,a,qmv,ifit,npar,fistat,imaxm)
C     calculation of the eigenvalues and matrixelements 
C     of the internal rotation part
      implicit none
      include 'iam.fi'
      integer ib,imaxm
      integer ifit(DIMPAR),npar,fistat
      integer qmv(DIMV)
      real*8  a(DIMPAR)
      real*8  h(DIMTOT,DIMTOT)
      real*8  evalv(DIMV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  ovv(DIMV,DIMV,DIMOVV,-DIMSIG:DIMSIG,-DIMJ:DIMJ,DIMTOP)
      real*8  rotm(-DIMJ:DIMJ,-DIMJ:DIMJ,1:2,DIMTOP)
      real*8  rott(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,DIMTOP)
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)

C     work
      real*8  am(DIMPM),ai(DIMPIR)
      real*8  mvec(DIMM,DIMV,-DIMJ:DIMJ)
      integer imfit(DIMOVV),sigma,maxm
      integer i,itop,k,ip,ivc,ivr,isig,im
      integer sdone(-20:20,DIMTOP)
      integer myand
      external myand

      do itop=1,ctlint(C_NTOP)
        do isig=-20,20
          sdone(isig,itop)=0
        end do
      end do
      do isig=1,size(S_G)
        do itop=1,ctlint(C_NTOP) 
          if (myand(ctlint(C_WOODS),4).eq.0) then
            size(S_FIRV+itop)=size(S_MINV+itop)
            size(S_MAXV+itop)=size(S_V+itop)
          else
            size(S_FIRV+itop)=1
            if (size(S_MAXV+itop).eq.0)
     $          size(S_MAXV+itop)=size(S_V+itop)
          end if
          if (size(S_FIRV+itop).le.0) then
            write(*,'(A,I3)') 'FIRST V: spezify V lines for B=',ib
            stop 'FIRST V < 0 error'
          end if
          sigma=gamma(isig,itop)
          if (sdone(sigma,itop).eq.0) then
            sdone(sigma,itop)=1
            if (imaxm.eq.0) then
              maxm=size(S_MAXM+itop)
            end if
            if (imaxm.lt.0) then
              maxm=size(S_MINV+itop)-1+size(S_V+itop)
            end if
            if (imaxm.gt.0) then
              maxm=imaxm
            end if
            do i=1, DIMPM 
              am(i)=a(DIMPRR+(itop-1)*DIMPIR+i)
              imfit(i)=ifit(DIMPRR+(itop-1)*DIMPIR+i)
            end do
            do i=1, DIMPIR
                ai(i)=a(DIMPRR+(itop-1)*DIMPIR+i)
            end do
            imfit(PM_PI)=ifit(P_FF)
            if (a(P_FF).ne.0.0) imfit(PM_PI)=1 
            imfit(PM_COS)=ifit(P_VCC)
            if (a(P_VCC).ne.0.0) imfit(PM_COS)=1 
            imfit(PM_SIN)=ifit(P_VSS)
            if (a(P_VSS).ne.0.0) imfit(PM_SIN)=1
            do k=-size(S_MAXK),size(S_MAXK)
              call calcm(sigma,h
     $             ,evalv(1,sigma,k,itop)
     $             ,ovv(1,1,1,sigma,k,itop) ! this sends ovv(:,:,:,sigma,k,itop) to calcm, which then uses ovv(DIMV,DIMV,DIMOVV) to calculate for the first three entries.
     $             ,mvec(1,1,k)
     $             ,am
     $             ,qmv 
     $             ,imfit
     $             ,k
     $             ,maxm
     $             ,size(S_FIRV+itop)
     $             ,size(S_MAXV+itop),ai)
            end do
            if ((myand(ctlint(C_PRI),AP_MO).ne.0).and.(xde.ge.1)) then
              write(*,'(A,5I3)') ' itop,sigma,k,B,fistat'
     $             ,itop,sigma,k,ib,fistat
              do ip=1,DIMOVV
                do ivr= 1, size(S_MAXV+itop)
                  write(*,'(I3,A,30F11.6)') ip,' ovv',
     $                 ((ovv(ivr,ivc,ip,sigma,k,itop)
     $                 ,k=-size(S_MAXK),size(S_MAXK))
     $                 ,ivc=1, size(S_MAXV+itop))
                end do
              end do
              write(*,*) 'mvec'
              do im=1,2*size(S_MAXM+itop)+1 
                write(*,'(40F10.6)') ((mvec(im,ivc,k)
     $               ,k=-size(S_MAXK),size(S_MAXK))
     $               ,ivc=1, size(S_MAXV+itop))
              end do
            end if
            call caltori(sigma,itop,mvec,tori,fistat
     $           ,maxm,size(S_FIRV+itop), size(S_MAXV+itop))
          end if
        end do ! itop
      end do ! isig

      if ((myand(ctlint(C_PRI),AP_EO).ne.0).and.(xde.ge.1)) then
        do itop=1, ctlint(C_NTOP)
          write(*,'(/,A,I3)')' Eigenvalues of one Top.  B=',ib
          do k=-size(S_MAXK),size(S_MAXK)
            do i=1, size(S_MAXV+itop)
              do isig=1, size(S_G)
                if (gamma(isig,itop).eq.-NaQN) goto 10
                write(*,'(F15.6,$)')
     $               evalv(i,gamma(isig,itop),k,itop)
              end do
 10           continue
            end do
            write(*,*)
          end do
        end do
      end if
      return
      end

C----------------------------------------------------------------------
      subroutine caltori(sigma,itop,mvec,tori,fistat,maxm,minv,sizv)
C     calculation of the eigenvalues of one matrix with specified j,k,m
C     the evalues are put in the field of dnv(1..ndata,NV_ENG,Q_UP/LO)
C     the deviations DE/DPi in dnv(1..ndata,2-DIMPAR,Q_UP/LO(i))
      
      implicit none
      include 'iam.fi'
      integer sigma,itop,fistat,maxm,minv,sizv
      real*8  tori(-DIMJ:DIMJ,-DIMJ:DIMJ,DIMV,DIMV,
     $     -DIMSIG:DIMSIG,DIMTOP)
      real*8  mvec(DIMM,DIMV,-DIMJ:DIMJ)
      
      integer kr,kc,vr,vc,im
      real*8  sum1,sum2, scale
      integer myand
      external myand

      if (myand(ctlint(C_WOODS),1).ne.0) then
C     calc the torsional integrals
        do kr=-size(S_MAXK), size(S_MAXK)
          do kc=-size(S_MAXK), size(S_MAXK)
            do vr=1, sizv
              do vc=1, sizv
                tori(kr,kc,vr,vc,sigma,itop)=0.0d0
                do im=1, 2*maxm+1
                  tori(kr,kc,vr,vc,sigma,itop) =
     $                 tori(kr,kc,vr,vc,sigma,itop)
     $                 +mvec(im,vr,kr)*mvec(im,vc,kc)
                end do
              end do
            end do
          end do
        end do

        do kr=-size(S_MAXK), size(S_MAXK)
          do vr=1,sizv 
            do vc=vr+1,sizv 
              tori(kr,kr,vr,vc,sigma,itop)=0.0d0
              tori(kr,kr,vc,vr,sigma,itop)=0.0d0
            end do
          end do
          do vr=1,sizv 
            tori(kr,kr,vr,vr,sigma,itop)=1.0d0
          end do
        end do

        if ((myand(ctlint(C_PRI),AP_TI).ne.0).and.(xde.ge.1)) then
          write(*,'(/,2(A,I3,3X))') ' torsional integrals  sigma=',sigma
     $         ,'deriv.=',fistat
          do vr=1,sizv 
            do kr=-size(S_MAXK), size(S_MAXK)
              write(*,'(20F8.5)')
     $             ((tori(kr,kc,vr,vc,sigma,itop)
     $             ,kc=-size(S_MAXK), size(S_MAXK))
     $             ,vc=1,sizv)
            end do
          end do
        end if
      end if

      if (myand(ctlint(C_WOODS),1).eq.0) then     
        do kr=-size(S_MAXK), size(S_MAXK)
          do kc=-size(S_MAXK), size(S_MAXK)
            do vr=1,sizv 
              do vc=1,sizv
                tori(kr,kc,vr,vc,sigma,itop)=0.0d0
              end do
            end do
          end do
        end do
        do kr=-size(S_MAXK), size(S_MAXK)
          do kc=-size(S_MAXK), size(S_MAXK)
            do vr=1,sizv
              tori(kr,kc,vr,vr,sigma,itop)=1.0d0
            end do
          end do
        end do
      end if

      if (myand(ctlint(C_WOODS),2).ne.0) then
        do kr=-size(S_MAXK), size(S_MAXK)
          do kc= -size(S_MAXK), size(S_MAXK)
            do vr=1,sizv
              sum1=0.0d0
              sum2=0.0d0
              do vc=1,vr-1
                sum1=sum1+tori(kr,kc,vr,vc,sigma,itop)**2
              end do
              do vc=vr,sizv
                sum2=sum2+tori(kr,kc,vr,vc,sigma,itop)**2
              end do
              do vc=vr,sizv
                scale=(dsqrt(1.0d0-sum1))/dsqrt(sum2)
                tori(kr,kc,vr,vc,sigma,itop)=
     $               tori(kr,kc,vr,vc,sigma,itop)*scale
                if (vr.ne.vc) 
     $               tori(kc,kr,vc,vr,sigma,itop)=
     $               tori(kr,kc,vr,vc,sigma,itop)
              end do
            end do
          end do
        end do
        if ((myand(ctlint(C_PRI),AP_TI).ne.0).and.(xde.ge.1)) then
          write(*,'(/,A)') '----- scaled torsional integrals -----'
          do vr=1,sizv
            do kr=-size(S_MAXK), size(S_MAXK)
              write(*,'(20F8.5)')
     $             ((tori(kr,kc,vr,vc,sigma,itop)
     $             ,kc=-size(S_MAXK), size(S_MAXK))
     $             ,vc=1,sizv)
            end do
          end do
        end if
      end if
      return
      end

C ---------------------------------------------------------------------
      subroutine maxof(mat, dimrow, dimcol, nrow, fcol, ncol,
     $     besti, scndi, bestv, scndv)
C     find the biggest and second abs values in column=fcol..fcol+ncol-1
      implicit none
      integer dimrow,dimcol,nrow,fcol,ncol
      integer besti(ncol), scndi(ncol)
      real*8  mat(dimrow,dimcol), bestv(ncol), scndv(ncol)
C     work
      integer ic,ir,i

      if ((ncol+fcol-1).gt.dimcol)  stop ' ERROR: DIMCOL  exceeded'
      if (nrow.gt.dimrow)  stop ' ERROR: DIMROW  exceeded'
      do ic=1, ncol
        besti(ic)=1
        scndi(ic)=0
        bestv(ic)=0.0d0
        scndv(ic)=0.0d0
      end do

      do ir=1, nrow
        do i=1, ncol
          ic=i+fcol-1
          if (scndv(i).le.dabs(mat(ir,ic))) then
            if (bestv(i).le.dabs(mat(ir,ic))) then
              scndi(i)=besti(i)
              scndv(i)=bestv(i)
              besti(i)=ir
              bestv(i)=dabs(mat(ir,ic))
            else
              scndi(i)=ir
              scndv(i)=dabs(mat(ir,ic))
            end if
          end if          
        end do
      end do
      do i=1, ncol
        ic=i+fcol-1
        bestv(i)=mat(besti(i),ic)
        scndv(i)=mat(scndi(i),ic)
      end do
      return
      end

C----------------------------------------------------------------------
      subroutine savepsi(i,psi)
      implicit none
      real*8  psi
      integer i
      include 'iam.fi'
      dln(i,LN_PSI)=psi
      return
      end

C--------------------------------------------------------------
      subroutine swptdo(i,j)
      implicit none
      include 'iam.fi'
      integer i,j,itmp
      
      itmp=todo(i,Q_J)
      todo(i,Q_J)=todo(j,Q_J)
      todo(j,Q_J)=itmp
      itmp=todo(i,Q_S)
      todo(i,Q_S)=todo(j,Q_S)
      todo(j,Q_S)=itmp
      itmp=todo(i,Q_B)
      todo(i,Q_B)=todo(j,Q_B)
      todo(j,Q_B)=itmp
      itmp=todo(i,Q_F)
      todo(i,Q_F)=todo(j,Q_F)
      todo(j,Q_F)=itmp
      itmp=todo(i,Q_F1)
      todo(i,Q_F1)=todo(j,Q_F1) ! h2024, adding F1 to swapping funciton
      todo(j,Q_F1)=itmp         ! h2024, adding F1 to swapping funciton
      itmp=todo(i,Q_STAT)
      todo(i,Q_STAT)=todo(j,Q_STAT)
      todo(j,Q_STAT)=itmp
      return
      end

C=====================================================================
      block data bdata
C ---------------------------------------------------------------------
      include 'iam.fi'
      include 'iamdata.fi'
      integer idat
 	
      data gamstr(1)/'G'/
      data gamstr(2)/'S'/
      data gamstr(3)/'G'/
      data gamstr(4)/'S'/
      data gamstr(5)/'G'/
      data gamstr(6)/'S'/
      data gamstr(7)/'G'/!Herbers2024
      data gamstr(8)/'S'/!Herbers2024

C      data gamstr /DIMTOP*'G'/

      data parstr(P_BJ    ) /'BJ      '/, parfit(P_BJ    ) /0/
      data parstr(P_BK    ) /'BK      '/, parfit(P_BK    ) /0/
      data parstr(P_BD    ) /'B-      '/, parfit(P_BD    ) /0/
      data parstr(P_DZX   ) /'Dzx     '/, parfit(P_DZX   ) /0/
      data parstr(P_DZXJ  ) /'DzxJ    '/, parfit(P_DZXJ  ) /0/
      data parstr(P_DZXK  ) /'DzxK    '/, parfit(P_DZXK  ) /0/
      data parstr(P_DJ    ) /'DJ      '/, parfit(P_DJ    ) /0/
      data parstr(P_DJK   ) /'DJK     '/, parfit(P_DJK   ) /0/
      data parstr(P_DK    ) /'DK      '/, parfit(P_DK    ) /0/
      data parstr(P_DJD   ) /'dj      '/, parfit(P_DJD   ) /0/
      data parstr(P_DKD   ) /'dk      '/, parfit(P_DKD   ) /0/
C     data parstr(P_R6    ) /'R6      '/, parfit(P_R6    ) /0/
      data parstr(P_HJ    ) /'H_J     '/, parfit(P_HJ    ) /0/
      data parstr(P_HJK   ) /'HJK     '/, parfit(P_HJK   ) /0/
      data parstr(P_HKJ   ) /'HKJ     '/, parfit(P_HKJ   ) /0/
      data parstr(P_HK    ) /'H_K     '/, parfit(P_HK    ) /0/
      data parstr(P_HJD   ) /'h_j     '/, parfit(P_HJD   ) /0/
      data parstr(P_HJKD  ) /'hjk     '/, parfit(P_HJKD  ) /0/
      data parstr(P_HKD   ) /'h_k     '/, parfit(P_HKD   ) /0/ !/1/ means analytic gradient not available
      data parstr(P_LJ    ) /'L_J     '/, parfit(P_LJ    ) /0/
      data parstr(P_LJJK  ) /'LJJK    '/, parfit(P_LJJK  ) /0/
      data parstr(P_LJK   ) /'LJK     '/, parfit(P_LJK   ) /0/
      data parstr(P_LKKJ  ) /'LKKJ    '/, parfit(P_LKKJ  ) /0/
      data parstr(P_LK    ) /'L_K     '/, parfit(P_LK    ) /0/
      data parstr(P_LJD   ) /'l_j     '/, parfit(P_LJD   ) /0/
      data parstr(P_LJKD  ) /'ljk     '/, parfit(P_LJKD  ) /0/
      data parstr(P_LKJD  ) /'lkj     '/, parfit(P_LKJD  ) /0/
      data parstr(P_LKD   ) /'l_k     '/, parfit(P_LKD   ) /0/
      data parstr(P_E     ) /'E       '/, parfit(P_E     ) /1/ !Herbers2023
      data parstr(P_GX12  ) /'Gx12    '/, parfit(P_GX12  ) /1/ !Herbers2024 Wilsons coriolis coupling   The parameterts were tested against SPFIT SPCAT, while all combinations of Gx,y,z OR Fxy,xz,yz reproduced the 
      data parstr(P_GY12  ) /'Gy12    '/, parfit(P_GY12  ) /1/ !Herbers2024 Wilsons coriolis coupling   simulated spectra, contradictions where found when Mixing the Gs and the Fs. I could not find the error, so I 
      data parstr(P_GZ12  ) /'Gz12    '/, parfit(P_GZ12  ) /1/ !Herbers2024 Wilsons coriolis coupling   do not recommend mixing Gs and Fs at the moment.
      data parstr(P_FXY1  ) /'Fxy12   '/, parfit(P_FXY1  ) /1/ !Herbers2024 Picketts coriolis coupling 
      data parstr(P_FYZ1  ) /'Fyz12   '/, parfit(P_FYZ1  ) /1/ !Herbers2024 Picketts coriolis coupling 
      data parstr(P_FXZ1  ) /'Fxz12   '/, parfit(P_FXZ1  ) /1/ !Herbers2024 Picketts coriolis coupling 
      data parstr(P_GX34  ) /'Gx34    '/, parfit(P_GX34  ) /1/ !Herbers2024 Wilsons coriolis coupling  
      data parstr(P_GY34  ) /'Gy34    '/, parfit(P_GY34  ) /1/ !Herbers2024 Wilsons coriolis coupling  
      data parstr(P_GZ34  ) /'Gz34    '/, parfit(P_GZ34  ) /1/ !Herbers2024 Wilsons coriolis coupling  
      data parstr(P_FXY3  ) /'Fxy34   '/, parfit(P_FXY3  ) /1/ !Herbers2024 Picketts coriolis coupling 
      data parstr(P_FYZ3  ) /'Fyz34   '/, parfit(P_FYZ3  ) /1/ !Herbers2024 Picketts coriolis coupling 
      data parstr(P_FXZ3  ) /'Fxz34   '/, parfit(P_FXZ3  ) /1/ !Herbers2024 Picketts coriolis coupling       
      data parstr(P_DX1   ) /'DxS1    '/, parfit(P_DX1   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DX2   ) /'DxS2    '/, parfit(P_DX2   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DX3   ) /'DxS3    '/, parfit(P_DX3   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DX4   ) /'DxS4    '/, parfit(P_DX4   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DX5   ) /'DxS5    '/, parfit(P_DX5   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5      
      data parstr(P_DX6   ) /'DxS6    '/, parfit(P_DX6   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DX7   ) /'DxS7    '/, parfit(P_DX7   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DX8   ) /'DxS8    '/, parfit(P_DX8   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DX9   ) /'DxS9    '/, parfit(P_DX9   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DX10  ) /'DxS10   '/, parfit(P_DX10  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DX11  ) /'DxS11   '/, parfit(P_DX11  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DY1   ) /'DyS1    '/, parfit(P_DY1   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DY2   ) /'DyS2    '/, parfit(P_DY2   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DY3   ) /'DyS3    '/, parfit(P_DY3   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DY4   ) /'DyS4    '/, parfit(P_DY4   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DY5   ) /'DyS5    '/, parfit(P_DY5   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5   
      data parstr(P_DY6   ) /'DyS6    '/, parfit(P_DY6   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DY7   ) /'DyS7    '/, parfit(P_DY7   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DY8   ) /'DyS8    '/, parfit(P_DY8   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DY9   ) /'DyS9    '/, parfit(P_DY9   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DY10  ) /'DyS10   '/, parfit(P_DY10  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DY11  ) /'DyS11   '/, parfit(P_DY11  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ1   ) /'DzS1    '/, parfit(P_DZ1   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ2   ) /'DzS2    '/, parfit(P_DZ2   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DZ3   ) /'DzS3    '/, parfit(P_DZ3   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DZ4   ) /'DzS4    '/, parfit(P_DZ4   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5
      data parstr(P_DZ5   ) /'DzS5    '/, parfit(P_DZ5   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ6   ) /'DzS6    '/, parfit(P_DZ6   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ7   ) /'DzS7    '/, parfit(P_DZ7   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ8   ) /'DzS8    '/, parfit(P_DZ8   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DZ9   ) /'DzS9    '/, parfit(P_DZ9   ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5  
      data parstr(P_DZ10  ) /'DzS10   '/, parfit(P_DZ10  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5 
      data parstr(P_DZ11  ) /'DzS11   '/, parfit(P_DZ11  ) /1/ !Same as P_x P_y P_z by hartwig but separate for S1,2,3,4,5       
      data parstr(P_DUMP  ) /'dummy   '/, parfit(P_DUMP  ) /1/ !Herbers2024 A dummy parameter not used for anything.                                             
      data parstr(P_QZ    ) /'chi_z   '/, parfit(P_QZ    ) /0/
      data parstr(P_QD    ) /'chi_-   '/, parfit(P_QD    ) /0/
      data parstr(P_QXY   ) /'chi_xy  '/, parfit(P_QXY   ) /0/
      data parstr(P_QXZ   ) /'chi_xz  '/, parfit(P_QXZ   ) /0/
      data parstr(P_QYZ   ) /'chi_yz  '/, parfit(P_QYZ   ) /0/
      data parstr(P_Q2Z   ) /'chi2_z  '/, parfit(P_Q2Z   ) /0/ !h2024, parameters of 2nd nucleus.
      data parstr(P_Q2D   ) /'chi2_-  '/, parfit(P_Q2D   ) /0/ !h2024, parameters of 2nd nucleus.
      data parstr(P_Q2XY  ) /'chi2_xy '/, parfit(P_Q2XY  ) /0/ !h2024, parameters of 2nd nucleus.
      data parstr(P_Q2XZ  ) /'chi2_xz '/, parfit(P_Q2XZ  ) /0/ !h2024, parameters of 2nd nucleus.
      data parstr(P_Q2YZ  ) /'chi2_yz '/, parfit(P_Q2YZ  ) /0/ !h2024, parameters of 2nd nucleus.   
      data parstr(P_FF    ) /'F12     '/, parfit(P_FF    ) /1/
      data parstr(P_VSS   ) /'Vss     '/, parfit(P_VSS   ) /1/
      data parstr(P_VCC   ) /'Vcc     '/, parfit(P_VCC   ) /1/
                                             
      data parstr(P_CP    ) /'C+      '/, parfit(P_CP    ) /0/
      data parstr(P_CZ    ) /'C_z     '/, parfit(P_CZ    ) /0/
      data parstr(P_CD    ) /'C-      '/, parfit(P_CD    ) /0/
                                      
      data parstr(P_MUX   ) /'mu_x    '/, parfit(P_MUX   ) /9/
      data parstr(P_MUY   ) /'mu_y    '/, parfit(P_MUY   ) /9/
      data parstr(P_MUZ   ) /'mu_z    '/, parfit(P_MUZ   ) /9/
                                             
      data parstr(P_PX    ) /'P_x     '/, parfit(P_PX    ) /0/
      data parstr(P_PY    ) /'P_y     '/, parfit(P_PY    ) /0/
      data parstr(P_PZ    ) /'P_z     '/, parfit(P_PZ    ) /0/
                                            
      data parstr(P1_VN1  ) /'V1n_1   '/, parfit(P1_VN1  ) /1/
      data parstr(P1_VN2  ) /'V2n_1   '/, parfit(P1_VN2  ) /1/
      data parstr(P1_F    ) /'F_1     '/, parfit(P1_F    ) /1/
      data parstr(P1_RHO  ) /'rho_1   '/, parfit(P1_RHO  ) /1/
      data parstr(P1_BETA ) /'beta_1  '/, parfit(P1_BETA ) /1/
      data parstr(P1_GAMA ) /'gamma_1 '/, parfit(P1_GAMA ) /1/
      data parstr(P1_DPIJ ) /'Dpi2J_1 '/, parfit(P1_DPIJ ) /1/
      data parstr(P1_DPIK ) /'Dpi2K_1 '/, parfit(P1_DPIK ) /1/
      data parstr(P1_DPID ) /'Dpi2-_1 '/, parfit(P1_DPID ) /1/
      data parstr(P1_DP4J ) /'Dpi4J_1 '/, parfit(P1_DP4J ) /1/!Herbers2026
      data parstr(P1_DP4K ) /'Dpi4K_1 '/, parfit(P1_DP4K ) /1/!Herbers2026
      data parstr(P1_DP4D ) /'Dpi4-_1 '/, parfit(P1_DP4D ) /1/!Herbers2026
      data parstr(P1_DC3K ) /'Dc3K_1  '/, parfit(P1_DC3K ) /1/ !Herbers2018
      data parstr(P1_DC3D ) /'Dc3-_1  '/, parfit(P1_DC3D ) /1/ !Herbers2018
      data parstr(P1_DC3J ) /'Dc3J_1  '/, parfit(P1_DC3J ) /1/
      data parstr(P1_D3K2 ) /'D3KK_1  '/, parfit(P1_D3K2 ) /1/ !Herbers2024
      data parstr(P1_D3J2 ) /'D3JJ_1  '/, parfit(P1_D3J2 ) /1/ !Herbers2026
      data parstr(P1_D3_J ) /'D3-j_1  '/, parfit(P1_D3_J ) /1/ !Herbers2026
      data parstr(P1_D3_K ) /'D3-k_1  '/, parfit(P1_D3_K ) /1/ !Herbers2026
      data parstr(P1_D3JK ) /'D3JK_1  '/, parfit(P1_D3JK ) /1/ !Herbers2026
      data parstr(P1_DPK2 ) /'Dp2KK_1 '/, parfit(P1_DPK2 ) /1/ !Herbers2024
      data parstr(P1_DPJ2 ) /'Dp2JJ_1 '/, parfit(P1_DPJ2 ) /1/ !Herbers2026
      data parstr(P1_DP_J ) /'Dp2-j_1 '/, parfit(P1_DP_J ) /1/ !Herbers2026
      data parstr(P1_DPJK ) /'Dp2JK_1 '/, parfit(P1_DPJK ) /1/ !Herbers2026
      data parstr(P1_DP_K ) /'Dp2-k_1 '/, parfit(P1_DP_K ) /1/ !Herbers2026
      data parstr(P1_DP4K2 ) /'Dp4KK_1 '/, parfit(P1_DP4K2 ) /1/ !Herbers2024
      data parstr(P1_DP4J2 ) /'Dp4JJ_1 '/, parfit(P1_DP4J2 ) /1/ !Herbers2026
      data parstr(P1_DP4_J ) /'Dp4-j_1 '/, parfit(P1_DP4_J ) /1/ !Herbers2026
      data parstr(P1_DP4JK ) /'Dp4JK_1 '/, parfit(P1_DP4JK ) /1/ !Herbers2026
      data parstr(P1_DP4_K ) /'Dp4-k_1 '/, parfit(P1_DP4_K ) /1/ !Herbers2026
      data parstr(P1_DFM2 ) /'DFm2_1  '/, parfit(P1_DFM2 ) /1/ !Herbers2026
      data parstr(P1_FMK  ) /'Fmk_1   '/, parfit(P1_FMK  ) /1/ !Herbers2026
      data parstr(P1_FMKJ ) /'FmkJ_1  '/, parfit(P1_FMKJ ) /1/ !Herbers2026
      data parstr(P1_FMKK ) /'FmkK_1  '/, parfit(P1_FMKK ) /1/ !Herbers2026
      data parstr(P1_FMKD ) /'Fmk-_1  '/, parfit(P1_FMKD ) /1/ !Herbers2026
      data parstr(P1_FMKJ2) /'FmkJJ_1 '/, parfit(P1_FMKJ2) /1/ !Herbers2026
      data parstr(P1_FMKK2) /'FmkKK_1 '/, parfit(P1_FMKK2) /1/ !Herbers2026
      data parstr(P1_FMKJK) /'FmkJK_1 '/, parfit(P1_FMKJK) /1/ !Herbers2026
      data parstr(P1_FMK_J) /'Fmk-j_1 '/, parfit(P1_FMK_J) /1/ !Herbers2026
      data parstr(P1_FMK_K) /'Fmk-k_1 '/, parfit(P1_FMK_K) /1/ !Herbers2026
      data parstr(P1_MK3J ) /'mk3J_1  '/, parfit(P1_MK3J ) /1/ !Herbers2026
      data parstr(P1_MK3K ) /'mk3K_1  '/, parfit(P1_MK3K ) /1/ !Herbers2026
      data parstr(P1_MK3D ) /'mk3-_1  '/, parfit(P1_MK3D ) /1/ !Herbers2026
      data parstr(P1_M3KJ ) /'m3kJ_1  '/, parfit(P1_M3KJ ) /1/ !Herbers2026
      data parstr(P1_M3KK ) /'m3kK_1  '/, parfit(P1_M3KK ) /1/ !Herbers2026
      data parstr(P1_M3KD ) /'m3k-_1  '/, parfit(P1_M3KD ) /1/ !Herbers2026
      data parstr(P1_MK3  ) /'mk3_1   '/, parfit(P1_MK3  ) /1/ !Herbers2026
      data parstr(P1_M3K  ) /'m3k_1   '/, parfit(P1_M3K  ) /1/ !Herbers2026
      data parstr(P1_DPI4 ) /'Dpi4_1  '/, parfit(P1_DPI4 ) /1/ !Herbers2026
      data parstr(P1_MK   ) /'mk_1    '/, parfit(P1_MK   ) /1/ !Herbers2026
      data parstr(P1_MKJ  ) /'mkJ_1   '/, parfit(P1_MKJ  ) /1/ !Herbers2026
      data parstr(P1_MKK  ) /'mkK_1   '/, parfit(P1_MKK  ) /1/ !Herbers2026
      data parstr(P1_MKD  ) /'mkD_1   '/, parfit(P1_MKD  ) /1/ !Herbers2026
      data parstr(P1_M2K2 ) /'m2k2_1  '/, parfit(P1_M2K2 ) /1/ !Herbers2026
      data parstr(P1_RHOJ ) /'rhoJ_1  '/, parfit(P1_RHOJ ) /1/ !Herbers2026
      data parstr(P1_RHOK ) /'rhoK_1  '/, parfit(P1_RHOK ) /1/ !Herbers2026
      data parstr(P1_RHOD ) /'rho-_1  '/, parfit(P1_RHOD ) /1/ !Herbers2026
      data parstr(P1_RHOJ2) /'rhoJJ_1 '/, parfit(P1_RHOJ2) /1/ !Herbers2026
      data parstr(P1_RHOK2) /'rhoKK_1 '/, parfit(P1_RHOK2) /1/ !Herbers2026
      data parstr(P1_RHOJK) /'rhoJK_1 '/, parfit(P1_RHOJK) /1/ !Herbers2026
      data parstr(P1_RHO_J) /'rho-j_1 '/, parfit(P1_RHO_J) /1/ !Herbers2026
      data parstr(P1_RHO_K) /'rho-k_1 '/, parfit(P1_RHO_K) /1/ !Herbers2026
      data parstr(P1_F0   ) /'F0_1    '/, parfit(P1_F0   ) /1/
      data parstr(P1_ANGX ) /'epsil_1 '/, parfit(P1_ANGX ) /1/
      data parstr(P1_ANGZ ) /'delta_1 '/, parfit(P1_ANGZ ) /1/
      data parstr(P1_FM2K2) /'Fm2k2_1 '/, parfit(P1_FM2K2) /1/ !Herbers2026
      data parstr(P1_FMK3 )  /'Fmk3_1  '/,  parfit(P1_FMK3 ) /1/ !Herbers2026
      data parstr(P1_FM3K )  /'Fm3k_1  '/,  parfit(P1_FM3K ) /1/ !Herbers2026
      data parstr(P1_FK2  )   /'Fk2_1   '/, parfit(P1_FK2  ) /1/ !Herbers2026
      data parstr(P1_FM3K3) /'Fm3k3_1 '/,   parfit(P1_FM3K3) /1/ !Herbers2026
                                            
      data parstr(P2_VN1  ) /'V1n_2   '/, parfit(P2_VN1  ) /1/
      data parstr(P2_VN2  ) /'V2n_2   '/, parfit(P2_VN2  ) /1/
      data parstr(P2_F    ) /'F_2     '/, parfit(P2_F    ) /1/
      data parstr(P2_RHO  ) /'rho_2   '/, parfit(P2_RHO  ) /1/
      data parstr(P2_BETA ) /'beta_2  '/, parfit(P2_BETA ) /1/
      data parstr(P2_GAMA ) /'gamma_2 '/, parfit(P2_GAMA ) /1/
      data parstr(P2_DPIJ ) /'Dpi2J_2 '/, parfit(P2_DPIJ ) /1/
      data parstr(P2_DPIK ) /'Dpi2K_2 '/, parfit(P2_DPIK ) /1/
      data parstr(P2_DPID ) /'Dpi2-_2 '/, parfit(P2_DPID ) /1/
      data parstr(P2_DP4J ) /'Dpi4J_2 '/, parfit(P2_DP4J ) /1/!Herbers2026
      data parstr(P2_DP4K ) /'Dpi4K_2 '/, parfit(P2_DP4K ) /1/!Herbers2026
      data parstr(P2_DP4D ) /'Dpi4-_2 '/, parfit(P2_DP4D ) /1/!Herbers2026
      data parstr(P2_DC3K ) /'Dc3K_2  '/, parfit(P2_DC3K ) /1/ !Herbers2018
      data parstr(P2_DC3D ) /'Dc3-_2  '/, parfit(P2_DC3D ) /1/ !Herbers2018
      data parstr(P2_DC3J ) /'Dc3J_2  '/, parfit(P2_DC3J ) /1/
      data parstr(P2_D3K2 ) /'D3KK_2  '/, parfit(P2_D3K2 ) /1/ !Herbers2024
      data parstr(P2_D3J2 ) /'D3JJ_2  '/, parfit(P2_D3J2 ) /1/ !Herbers2026
      data parstr(P2_D3_J ) /'D3-j_2  '/, parfit(P2_D3_J ) /1/ !Herbers2026
      data parstr(P2_D3_K ) /'D3-k_2  '/, parfit(P2_D3_K ) /1/ !Herbers2026
      data parstr(P2_D3JK ) /'D3JK_2  '/, parfit(P2_D3JK ) /1/ !Herbers2026
      data parstr(P2_DPK2 ) /'Dp2KK_2 '/, parfit(P2_DPK2 ) /1/ !Herbers2024
      data parstr(P2_DPJ2 ) /'Dp2JJ_2 '/, parfit(P2_DPJ2 ) /1/ !Herbers2026
      data parstr(P2_DP_J ) /'Dp2-j_2 '/, parfit(P2_DP_J ) /1/ !Herbers2026
      data parstr(P2_DPJK ) /'Dp2JK_2 '/, parfit(P2_DPJK ) /1/ !Herbers2026
      data parstr(P2_DP_K ) /'Dp2-k_2 '/, parfit(P2_DP_K ) /1/ !Herbers2026
      data parstr(P2_DP4K2 ) /'Dp4KK_2 '/, parfit(P2_DP4K2 ) /1/ !Herbers2024
      data parstr(P2_DP4J2 ) /'Dp4JJ_2 '/, parfit(P2_DP4J2 ) /1/ !Herbers2026
      data parstr(P2_DP4_J ) /'Dp4-j_2 '/, parfit(P2_DP4_J ) /1/ !Herbers2026
      data parstr(P2_DP4JK ) /'Dp4JK_2 '/, parfit(P2_DP4JK ) /1/ !Herbers2026
      data parstr(P2_DP4_K ) /'Dp4-k_2 '/, parfit(P2_DP4_K ) /1/ !Herbers2026
      data parstr(P2_DFM2 ) /'DFm2_2   '/, parfit(P2_DFM2 ) /1/ !Herbers2026
      data parstr(P2_FMK ) /'Fmk_2   '/, parfit(P2_FMK ) /1/ !Herbers2026
      data parstr(P2_FMKJ ) /'FmkJ_2  '/, parfit(P2_FMKJ ) /1/ !Herbers2026
      data parstr(P2_FMKK ) /'FmkK_2  '/, parfit(P2_FMKK ) /1/ !Herbers2026
      data parstr(P2_FMKD ) /'Fmk-_2  '/, parfit(P2_FMKD ) /1/ !Herbers2026
      data parstr(P2_FMKJ2) /'FmkJJ_2 '/, parfit(P2_FMKJ2) /1/ !Herbers2026
      data parstr(P2_FMKK2) /'FmkKK_2 '/, parfit(P2_FMKK2) /1/ !Herbers2026
      data parstr(P2_FMKJK) /'FmkJK_2 '/, parfit(P2_FMKJK) /1/ !Herbers2026
      data parstr(P2_FMK_J) /'Fmk-j_2 '/, parfit(P2_FMK_J) /1/ !Herbers2026
      data parstr(P2_FMK_K) /'Fmk-k_2 '/, parfit(P2_FMK_K) /1/ !Herbers2026
      data parstr(P2_MK3J ) /'mk3J_2  '/, parfit(P2_MK3J ) /1/ !Herbers2026
      data parstr(P2_MK3K ) /'mk3K_2  '/, parfit(P2_MK3K ) /1/ !Herbers2026
      data parstr(P2_MK3D ) /'mk3-_2  '/, parfit(P2_MK3D ) /1/ !Herbers2026
      data parstr(P2_M3KJ ) /'m3kJ_2  '/, parfit(P2_M3KJ ) /1/ !Herbers2026
      data parstr(P2_M3KK ) /'m3kK_2  '/, parfit(P2_M3KK ) /1/ !Herbers2026
      data parstr(P2_M3KD ) /'m3k-_2  '/, parfit(P2_M3KD ) /1/ !Herbers2026
      data parstr(P2_MK3  ) /'mk3_2   '/, parfit(P2_MK3 ) /1/ !Herbers2026
      data parstr(P2_M3K  ) /'m3k_2   '/, parfit(P2_M3K ) /1/ !Herbers2026
      data parstr(P2_DPI4 ) /'Dpi4_2  '/, parfit(P2_DPI4 ) /1/ !Herbers2026
      data parstr(P2_MK   ) /'mk_2    '/, parfit(P2_MK   ) /1/ !Herbers2026
      data parstr(P2_MKJ  ) /'mkJ_2   '/, parfit(P2_MKJ  ) /1/ !Herbers2026
      data parstr(P2_MKK  ) /'mkK_2   '/, parfit(P2_MKK  ) /1/ !Herbers2026
      data parstr(P2_MKD  ) /'mkD_2   '/, parfit(P2_MKD  ) /1/ !Herbers2026
      data parstr(P2_M2K2 ) /'m2k2_2  '/, parfit(P2_M2K2 ) /1/ !Herbers2026
      data parstr(P2_F0   ) /'F0_2    '/, parfit(P2_F0   ) /1/
      data parstr(P2_RHOJ ) /'rhoJ_2  '/, parfit(P2_RHOJ ) /1/ !Herbers2026
      data parstr(P2_RHOK ) /'rhoK_2  '/, parfit(P2_RHOK ) /1/ !Herbers2026
      data parstr(P2_RHOD ) /'rho-_2  '/, parfit(P2_RHOD ) /1/ !Herbers2026
      data parstr(P2_RHOJ2) /'rhoJJ_2 '/, parfit(P2_RHOJ2) /1/ !Herbers2026
      data parstr(P2_RHOK2) /'rhoKK_2 '/, parfit(P2_RHOK2) /1/ !Herbers2026
      data parstr(P2_RHOJK) /'rhoJK_2 '/, parfit(P2_RHOJK) /1/ !Herbers2026
      data parstr(P2_RHO_J) /'rho-j_2 '/, parfit(P2_RHO_J) /1/ !Herbers2026
      data parstr(P2_RHO_K) /'rho-k_2 '/, parfit(P2_RHO_K) /1/ !Herbers2026
      data parstr(P2_ANGX ) /'epsil_2 '/, parfit(P2_ANGX ) /1/
      data parstr(P2_ANGZ ) /'delta_2 '/, parfit(P2_ANGZ ) /1/
      data parstr(P2_FM2K2) /'Fm2k2_2 '/,   parfit(P2_FM2K2) /1/ !Herbers2026
      data parstr(P2_FMK3 )  /'Fmk3_2  '/,  parfit(P2_FMK3 ) /1/ !Herbers2026
      data parstr(P2_FM3K )  /'Fm3k_2  '/,  parfit(P2_FM3K ) /1/ !Herbers2026
      data parstr(P2_FK2  )   /'Fk2_2   '/, parfit(P2_FK2  ) /1/ !Herbers2026
      data parstr(P2_FM3K3) /'Fm3k3_2 '/,   parfit(P2_FM3K3) /1/ !Herbers2026

      data parstr(P3_VN1  ) /'V1n_3   '/, parfit(P3_VN1  ) /1/
      data parstr(P3_VN2  ) /'V2n_3   '/, parfit(P3_VN2  ) /1/
      data parstr(P3_F    ) /'F_3     '/, parfit(P3_F    ) /1/
      data parstr(P3_RHO  ) /'rho_3   '/, parfit(P3_RHO  ) /1/
      data parstr(P3_BETA ) /'beta_3  '/, parfit(P3_BETA ) /1/
      data parstr(P3_GAMA ) /'gamma_3 '/, parfit(P3_GAMA ) /1/
      data parstr(P3_DPIJ ) /'Dpi2J_3 '/, parfit(P3_DPIJ ) /1/
      data parstr(P3_DPIK ) /'Dpi2K_3 '/, parfit(P3_DPIK ) /1/
      data parstr(P3_DPID ) /'Dpi2-_3 '/, parfit(P3_DPID ) /1/
      data parstr(P3_DP4J ) /'Dpi4J_3 '/, parfit(P3_DP4J ) /1/!Herbers2026
      data parstr(P3_DP4K ) /'Dpi4K_3 '/, parfit(P3_DP4K ) /1/!Herbers2026
      data parstr(P3_DP4D ) /'Dpi4-_3 '/, parfit(P3_DP4D ) /1/!Herbers2026
      data parstr(P3_DC3K ) /'Dc3K_3  '/, parfit(P3_DC3K ) /1/ !Herbers2018  
      data parstr(P3_DC3D ) /'Dc3-_3  '/, parfit(P3_DC3D ) /1/ !Herbers2018      
      data parstr(P3_DC3J ) /'Dc3J_3  '/, parfit(P3_DC3J ) /1/
      data parstr(P3_D3K2 ) /'D3KK_3  '/, parfit(P3_D3K2 ) /1/ !Herbers2024
      data parstr(P3_D3J2 ) /'D3JJ_3  '/, parfit(P3_D3J2 ) /1/ !Herbers2026
      data parstr(P3_D3_J ) /'D3-j_3  '/, parfit(P3_D3_J ) /1/ !Herbers2026
      data parstr(P3_D3_K ) /'D3-k_3  '/, parfit(P3_D3_K ) /1/ !Herbers2026
      data parstr(P3_D3JK ) /'D3JK_3  '/, parfit(P3_D3JK ) /1/ !Herbers2026
      data parstr(P3_DPK2 ) /'Dp2KK_3 '/, parfit(P3_DPK2 ) /1/ !Herbers2024
      data parstr(P3_DPJ2 ) /'Dp2JJ_3 '/, parfit(P3_DPJ2 ) /1/ !Herbers2026
      data parstr(P3_DP_J ) /'Dp2-j_3 '/, parfit(P3_DP_J ) /1/ !Herbers2026
      data parstr(P3_DPJK ) /'Dp2JK_3 '/, parfit(P3_DPJK ) /1/ !Herbers2026
      data parstr(P3_DP4_K ) /'Dp2-k_3 '/, parfit(P3_DP4_K ) /1/ !Herbers2026
      data parstr(P3_DP4K2 ) /'Dp4KK_3 '/, parfit(P3_DP4K2 ) /1/ !Herbers2024
      data parstr(P3_DP4J2 ) /'Dp4JJ_3 '/, parfit(P3_DP4J2 ) /1/ !Herbers2026
      data parstr(P3_DP4_J ) /'Dp4-j_3 '/, parfit(P3_DP4_J ) /1/ !Herbers2026
      data parstr(P3_DP4JK ) /'Dp4JK_3 '/, parfit(P3_DP4JK ) /1/ !Herbers2026
      data parstr(P3_DP4_K ) /'Dp4-k_3 '/, parfit(P3_DP4_K ) /1/ !Herbers2026
      data parstr(P3_DFM2 ) /'DFm2_3   '/, parfit(P3_DFM2 ) /1/ !Herbers2026
      data parstr(P3_FMK ) /'Fmk_3   '/, parfit(P3_FMK ) /1/ !Herbers2026
      data parstr(P3_FMKJ ) /'FmkJ_3  '/, parfit(P3_FMKJ ) /1/ !Herbers2026
      data parstr(P3_FMKK ) /'FmkK_3  '/, parfit(P3_FMKK ) /1/ !Herbers2026
      data parstr(P3_FMKD ) /'Fmk-_3  '/, parfit(P3_FMKD ) /1/ !Herbers2026
      data parstr(P3_FMKJ2) /'FmkJJ_3 '/, parfit(P3_FMKJ2) /1/ !Herbers2026
      data parstr(P3_FMKK2) /'FmkKK_3 '/, parfit(P3_FMKK2) /1/ !Herbers2026
      data parstr(P3_FMKJK) /'FmkJK_3 '/, parfit(P3_FMKJK) /1/ !Herbers2026
      data parstr(P3_FMK_J) /'Fmk-j_3 '/, parfit(P3_FMK_J) /1/ !Herbers2026
      data parstr(P3_FMK_K) /'Fmk-k_3 '/, parfit(P3_FMK_K) /1/ !Herbers2026
      data parstr(P3_MK3J ) /'mk3J_3  '/, parfit(P3_MK3J ) /1/ !Herbers2026
      data parstr(P3_MK3K ) /'mk3K_3  '/, parfit(P3_MK3K ) /1/ !Herbers2026
      data parstr(P3_MK3D ) /'mk3-_3  '/, parfit(P3_MK3D ) /1/ !Herbers2026
      data parstr(P3_M3KJ ) /'m3kJ_3  '/, parfit(P3_M3KJ ) /1/ !Herbers2026
      data parstr(P3_M3KK ) /'m3kK_3  '/, parfit(P3_M3KK ) /1/ !Herbers2026
      data parstr(P3_M3KD ) /'m3k-_3  '/, parfit(P3_M3KD ) /1/ !Herbers2026
      data parstr(P3_MK3  ) /'mk3_3   '/, parfit(P3_MK3 ) /1/ !Herbers2026
      data parstr(P3_M3K  ) /'m3k_3   '/, parfit(P3_M3K ) /1/ !Herbers2026
      data parstr(P3_DPI4 ) /'Dpi4_3  '/, parfit(P3_DPI4 ) /1/ !Herbers2026
      data parstr(P3_MK   ) /'mk_3    '/, parfit(P3_MK   ) /1/ !Herbers2026
      data parstr(P3_MKJ  ) /'mkJ_3   '/, parfit(P3_MKJ  ) /1/ !Herbers2026
      data parstr(P3_MKK  ) /'mkK_3   '/, parfit(P3_MKK  ) /1/ !Herbers2026
      data parstr(P3_MKD  ) /'mkD_3   '/, parfit(P3_MKD  ) /1/ !Herbers2026
      data parstr(P3_M2K2 ) /'m2k2_3  '/, parfit(P3_M2K2 ) /1/ !Herbers2026
      data parstr(P3_RHOJ ) /'rhoJ_3  '/, parfit(P3_RHOJ ) /1/ !Herbers2026
      data parstr(P3_RHOK ) /'rhoK_3  '/, parfit(P3_RHOK ) /1/ !Herbers2026
      data parstr(P3_RHOD ) /'rho-_3  '/, parfit(P3_RHOD ) /1/ !Herbers2026
      data parstr(P3_RHOJ2) /'rhoJJ_3 '/, parfit(P3_RHOJ2) /1/ !Herbers2026
      data parstr(P3_RHOK2) /'rhoKK_3 '/, parfit(P3_RHOK2) /1/ !Herbers2026
      data parstr(P3_RHOJK) /'rhoJK_3 '/, parfit(P3_RHOJK) /1/ !Herbers2026
      data parstr(P3_RHO_J) /'rho-j_3 '/, parfit(P3_RHO_J) /1/ !Herbers2026
      data parstr(P3_RHO_K) /'rho-k_3 '/, parfit(P3_RHO_K) /1/ !Herbers2026
      data parstr(P3_F0   ) /'F0_3    '/, parfit(P3_F0   ) /1/
      data parstr(P3_ANGX ) /'epsil_3 '/, parfit(P3_ANGX ) /1/
      data parstr(P3_ANGZ ) /'delta_3 '/, parfit(P3_ANGZ ) /1/
      data parstr(P3_FM2K2) /'Fm2k2_3 '/,   parfit(P3_FM2K2) /1/ !Herbers2026
      data parstr(P3_FMK3 )  /'Fmk3_3  '/,  parfit(P3_FMK3 ) /1/ !Herbers2026
      data parstr(P3_FM3K )  /'Fm3k_3  '/,  parfit(P3_FM3K ) /1/ !Herbers2026
      data parstr(P3_FK2  )   /'Fk2_3   '/, parfit(P3_FK2  ) /1/ !Herbers2026
      data parstr(P3_FM3K3) /'Fm3k3_3 '/,   parfit(P3_FM3K3) /1/ !Herbers2026
      
      data parstr(P4_VN1  ) /'V1n_4   '/, parfit(P4_VN1  ) /1/ !Herbers2024
      data parstr(P4_VN2  ) /'V2n_4   '/, parfit(P4_VN2  ) /1/ !Herbers2024
      data parstr(P4_F    ) /'F_4     '/, parfit(P4_F    ) /1/ !Herbers2024
      data parstr(P4_RHO  ) /'rho_4   '/, parfit(P4_RHO  ) /1/ !Herbers2024
      data parstr(P4_BETA ) /'beta_4  '/, parfit(P4_BETA ) /1/ !Herbers2024
      data parstr(P4_GAMA ) /'gamma_4 '/, parfit(P4_GAMA ) /1/ !Herbers2024
      data parstr(P4_DPIJ ) /'Dpi2J_4 '/, parfit(P4_DPIJ ) /1/ !Herbers2024
      data parstr(P4_DPIK ) /'Dpi2K_4 '/, parfit(P4_DPIK ) /1/ !Herbers2024
      data parstr(P4_DPID ) /'Dpi2-_4 '/, parfit(P4_DPID ) /1/ !Herbers2024
      data parstr(P4_DP4J ) /'Dpi4J_4 '/, parfit(P4_DP4J ) /1/ !Herbers2026
      data parstr(P4_DP4K ) /'Dpi4K_4 '/, parfit(P4_DP4K ) /1/ !Herbers2026
      data parstr(P4_DP4D ) /'Dpi4-_4 '/, parfit(P4_DP4D ) /1/ !Herbers2026
      data parstr(P4_DC3K ) /'Dc3K_4  '/, parfit(P4_DC3K ) /1/ !Herbers2024 
      data parstr(P4_DC3D ) /'Dc3-_4  '/, parfit(P4_DC3D ) /1/ !Herbers2024     
      data parstr(P4_DC3J ) /'Dc3J_4  '/, parfit(P4_DC3J ) /1/ !Herbers2024
      data parstr(P4_D3K2 ) /'D3KK_4  '/, parfit(P4_D3K2 ) /1/ !Herbers2024
      data parstr(P4_D3J2 ) /'D3JJ_4  '/, parfit(P4_D3J2 ) /1/ !Herbers2026
      data parstr(P4_D3_J ) /'D3-j_4  '/, parfit(P4_D3_J ) /1/ !Herbers2026
      data parstr(P4_D3_K ) /'D3-k_4  '/, parfit(P4_D3_K ) /1/ !Herbers2026
      data parstr(P4_D3JK ) /'D3JK_4  '/, parfit(P4_D3JK ) /1/ !Herbers2026
      data parstr(P4_DPK2 ) /'Dp2KK_4 '/, parfit(P4_DPK2 ) /1/ !Herbers2024
      data parstr(P4_DPJ2 ) /'Dp2JJ_4 '/, parfit(P4_DPJ2 ) /1/ !Herbers2026
      data parstr(P4_DP_J ) /'Dp2-j_4 '/, parfit(P4_DP_J ) /1/ !Herbers2026
      data parstr(P4_DPJK ) /'Dp2JK_4 '/, parfit(P4_DPJK ) /1/ !Herbers2026
      data parstr(P4_DP_K ) /'Dp2-k_4 '/, parfit(P4_DP_K ) /1/ !Herbers2026
      data parstr(P4_DP4K2 ) /'Dp4KK_4 '/, parfit(P4_DP4K2 ) /1/ !Herbers2024
      data parstr(P4_DP4J2 ) /'Dp4JJ_4 '/, parfit(P4_DP4J2 ) /1/ !Herbers2026
      data parstr(P4_DP4_J ) /'Dp4-j_4 '/, parfit(P4_DP4_J ) /1/ !Herbers2026
      data parstr(P4_DP4JK ) /'Dp4JK_4 '/, parfit(P4_DP4JK ) /1/ !Herbers2026
      data parstr(P4_DP4_K ) /'Dp4-k_4 '/, parfit(P4_DP4_K ) /1/ !Herbers2026
      data parstr(P4_DFM2 ) /'DFm2_4   '/, parfit(P4_DFM2 ) /1/ !Herbers2026
      data parstr(P4_FMK ) /'Fmk_4   '/, parfit(P4_FMK ) /1/ !Herbers2026
      data parstr(P4_FMKJ ) /'FmkJ_4  '/, parfit(P4_FMKJ ) /1/ !Herbers2026
      data parstr(P4_FMKK ) /'FmkK_4  '/, parfit(P4_FMKK ) /1/ !Herbers2026
      data parstr(P4_FMKD ) /'Fmk-_4  '/, parfit(P4_FMKD ) /1/ !Herbers2026
      data parstr(P4_FMKJ2) /'FmkJJ_4 '/, parfit(P4_FMKJ2) /1/ !Herbers2026
      data parstr(P4_FMKK2) /'FmkKK_4 '/, parfit(P4_FMKK2) /1/ !Herbers2026
      data parstr(P4_FMKJK) /'FmkJK_4 '/, parfit(P4_FMKJK) /1/ !Herbers2026
      data parstr(P4_FMK_J) /'Fmk-j_4 '/, parfit(P4_FMK_J) /1/ !Herbers2026
      data parstr(P4_FMK_K) /'Fmk-k_4 '/, parfit(P4_FMK_K) /1/ !Herbers2026
      data parstr(P4_MK3J ) /'mk3J_4  '/, parfit(P4_MK3J ) /1/ !Herbers2026
      data parstr(P4_MK3K ) /'mk3K_4  '/, parfit(P4_MK3K ) /1/ !Herbers2026
      data parstr(P4_MK3D ) /'mk3-_4  '/, parfit(P4_MK3D ) /1/ !Herbers2026
      data parstr(P4_M3KJ ) /'m3kJ_4  '/, parfit(P4_M3KJ ) /1/ !Herbers2026
      data parstr(P4_M3KK ) /'m3kK_4  '/, parfit(P4_M3KK ) /1/ !Herbers2026
      data parstr(P4_M3KD ) /'m3k-_4  '/, parfit(P4_M3KD ) /1/ !Herbers2026
      data parstr(P4_MK3  ) /'mk3_4   '/, parfit(P4_MK3 ) /1/ !Herbers2026
      data parstr(P4_M3K  ) /'m3k_4   '/, parfit(P4_M3K ) /1/ !Herbers2026
      data parstr(P4_DPI4 ) /'Dpi4_4  '/, parfit(P4_DPI4 ) /1/ !Herbers2026
      data parstr(P4_MK   ) /'mk_4    '/, parfit(P4_MK   ) /1/ !Herbers2026
      data parstr(P4_MKJ  ) /'mkJ_4   '/, parfit(P4_MKJ  ) /1/ !Herbers2026
      data parstr(P4_MKK  ) /'mkK_4   '/, parfit(P4_MKK  ) /1/ !Herbers2026
      data parstr(P4_MKD  ) /'mkD_4   '/, parfit(P4_MKD  ) /1/ !Herbers2026
      data parstr(P4_M2K2 ) /'m2k2_4  '/, parfit(P4_M2K2 ) /1/ !Herbers2026
      data parstr(P4_RHOJ ) /'rhoJ_4  '/, parfit(P4_RHOJ ) /1/ !Herbers2026
      data parstr(P4_RHOK ) /'rhoK_4  '/, parfit(P4_RHOK ) /1/ !Herbers2026
      data parstr(P4_RHOD ) /'rho-_4  '/, parfit(P4_RHOD ) /1/ !Herbers2026
      data parstr(P4_RHOJ2) /'rhoJJ_4 '/, parfit(P4_RHOJ2) /1/ !Herbers2026
      data parstr(P4_RHOK2) /'rhoKK_4 '/, parfit(P4_RHOK2) /1/ !Herbers2026
      data parstr(P4_RHOJK) /'rhoJK_4 '/, parfit(P4_RHOJK) /1/ !Herbers2026
      data parstr(P4_RHO_J) /'rho-j_4 '/, parfit(P4_RHO_J) /1/ !Herbers2026
      data parstr(P4_RHO_K) /'rho-k_4 '/, parfit(P4_RHO_K) /1/ !Herbers2026
      data parstr(P4_F0   ) /'F0_4    '/, parfit(P4_F0   ) /1/ !Herbers2024
      data parstr(P4_ANGX ) /'epsil_4 '/, parfit(P4_ANGX ) /1/ !Herbers2024
      data parstr(P4_ANGZ ) /'delta_4 '/, parfit(P4_ANGZ ) /1/ !Herbers2024
      data parstr(P4_FM2K2) /'Fm2k2_4 '/,   parfit(P4_FM2K2) /1/ !Herbers2026
      data parstr(P4_FMK3 )  /'Fmk3_4  '/,  parfit(P4_FMK3 ) /1/ !Herbers2026
      data parstr(P4_FM3K )  /'Fm3k_4  '/,  parfit(P4_FM3K ) /1/ !Herbers2026
      data parstr(P4_FK2  )   /'Fk2_4   '/, parfit(P4_FK2  ) /1/ !Herbers2026
      data parstr(P4_FM3K3) /'Fm3k3_4 '/,   parfit(P4_FM3K3) /1/ !Herbers2026

      data ctlstr(C_NZYK ) /'nzyk  '/
      data ctlstr(C_NCYCL) /'ncycl '/
      data ctlstr(C_PRI)   /'aprint'/
      data ctlstr(C_PRINT) /'print '/
      data ctlstr(C_XPR)   /'xprint'/
      data ctlstr(C_INTS ) /'ints  '/
      data ctlstr(C_ORGER) /'orger '/
      data ctlstr(C_EVAL ) /'eval  '/
      data ctlstr(C_DFRQ ) /'dfreq '/
      data ctlstr(C_MAXM ) /'maxm  '/
c      data ctlstr(C_MAXM1) /'maxm1 '/
c      data ctlstr(C_MAXM2) /'maxm2 '/
c      data ctlstr(C_MAXM3) /'maxm3 '/

      data ctlstr(C_MAXV)  /'maxvm '/
c      data ctlstr(C_MAXV1) /'maxvm1'/
c      data ctlstr(C_MAXV2) /'maxvm2'/
c      data ctlstr(C_MAXV3) /'maxvm3'/
      data ctlstr(C_WOODS) /'woods '/
c      data ctlstr(C_WOOD1) /'woods1'/
c      data ctlstr(C_WOOD2) /'woods2'/
c      data ctlstr(C_WOOD3) /'woods3'/
c      data ctlstr(C_WOOD4) /'woods4'/
      data ctlstr(C_NDATA) /'ndata '/
      data ctlstr(C_NFOLD) /'nfold '/
      data ctlstr(C_SPIN ) /'spin  '/
      data ctlstr(C_NTOP ) /'ntop  '/
      data ctlstr(C_ADJF ) /'adjf  '/
c      data ctlstr(C_ADJ1 ) /'adjf1 '/
c      data ctlstr(C_ADJ2 ) /'adjf2 '/
c      data ctlstr(C_ADJ3 ) /'adjf3 '/
c      data ctlstr(C_ADJ4 ) /'adjf4 '/
      data ctlstr(C_ROFIT) /'rofit '/
      data ctlstr(C_DEFER) /'defer '/
      data ctlstr(C_EPS  ) /'eps   '/
      data ctlstr(C_WEIGF) /'weigf '/
      data ctlstr(C_CNVG ) /'convg '/
      data ctlstr(C_LMBDA) /'lambda'/
      data ctlstr(C_FITSC) /'fitscl'/
      data ctlstr(C_FRQLO) /'freq_l'/
      data ctlstr(C_FRQUP) /'freq_h'/
      data ctlstr(C_INTLM) /'limit '/
      data ctlstr(C_SVDER) /'svderr'/
      data ctlstr(C_TEMP ) /'temp  '/
      data ctlstr(C_RED  ) /'reduct'/
      data ctlstr(C_QSUM ) /'qsum  '/
      data ctlstr(C_DW   ) /'ctrl'/
C      data ctlstr(C_NQI ) /'NQInt '/
      data ctlstr(C_SPIN2) /'spin2'/
      data ctlstr(C_SORT)  /'fsort '/

      data qnostr( 1     ) /'Jup  '/ 
      data qnostr( 2     ) /'K-up '/
      data qnostr( 3     ) /'K+up '/
      data qnostr( 4     ) /'Jlo  '/
      data qnostr( 5     ) /'K-lo '/
      data qnostr( 6     ) /'K+lo '/
      data qnostr( 7     ) /'=    '/
      data qnostr( 8     ) /'Err  '/
      data qnostr( 9     ) /'Sup  '/
      data qnostr(10     ) /'Slo  '/ 
      data qnostr(11     ) /'V1up '/
      data qnostr(12     ) /'V1lo '/
      data qnostr(13     ) /'V2up '/
      data qnostr(14     ) /'V2lo '/
      data qnostr(15     ) /'Bup '/
      data qnostr(16     ) /'Blo '/
      data qnostr(17     ) /'Fup  '/
      data qnostr(18     ) /'Flo  '/
      data qnostr(19     ) /'Tup  '/
      data qnostr(20     ) /'Tlo  '/
      data qnostr(21     ) /'tup  '/
      data qnostr(22     ) /'tlo  '/
      data qnostr(23     ) /'#    '/
      data qnostr(24     ) /'&    '/
      data qnostr(25     ) /'Kup  '/
      data qnostr(26     ) /'Klo  '/
      data qnostr(27     ) /'diff '/ 
      data qnostr(28     ) /'GHz  '/
      data qnostr(29     ) /'MHz  '/
      data qnostr(30     ) /'cm-1 '/
      data qnostr(31     ) /'F1up'/ !Herbers 2024
      data qnostr(32     ) /'F1lo'/ !Herbers 2024

      data (qnostr(MAXQC+idat),idat=1,DIMGAM) /DIMGAM*'     '/

      data (q_q(idat),idat=1,MAXQC) !translates the Q_X integers into positions of quantum numbers saved in qlin.
     $     /Q_J ,0   ,0   ,Q_J ,0   ,0 
     $     ,0   ,0   ,Q_S, Q_S, Q_V1,Q_V1,Q_V2,Q_V2
     $     ,Q_B, Q_B, Q_F, Q_F, Q_T ,Q_T ,Q_TJ,Q_TJ
     $     ,Q_REF,Q_AVG,  Q_K, Q_K
     $     ,0,0,0,0,Q_F1,Q_F1/ !Herbers 2024 added Q_F1
      data (q_q(MAXQC+idat),idat=1,DIMGAM) /DIMGAM*0/

      data (qul(idat),idat=1,MAXQC) /Q_UP,0   ,0   ,Q_LO,0   ,0   
     $     ,0   ,0   ,Q_UP,Q_LO,Q_UP,Q_LO,Q_UP,Q_LO
     $     ,Q_UP,Q_LO,Q_UP,Q_LO,Q_UP,Q_LO,Q_UP,Q_LO
     $     ,0    ,0    ,  Q_UP,Q_LO
     $     ,0,0,0,0,Q_UP,Q_LO/!added two zeros Herbers 2024
      data (qul(MAXQC+idat),idat=1,DIMGAM) /DIMGAM*0/

c     data  q_q /Q_J, Q_K, Q_J, Q_K
c    $          ,Q_V1,Q_V2,Q_B, Q_S, Q_V1,Q_V2,Q_B, Q_S
c    $          ,Q_F, Q_F, Q_T, Q_T, Q_TJ, Q_TJ, Q_REF, Q_AVG 
c    $          ,0,0,0,0,0,0,0,0,0,0,0/
c     data  qul /Q_UP,Q_UP,Q_LO,Q_LO
c    $          ,Q_UP,Q_UP,Q_UP, Q_UP,Q_LO,Q_LO,Q_LO, Q_LO
c    $          ,Q_UP,Q_LO,Q_UP,Q_LO,Q_UP, Q_LO, 0,    0
c    $          ,0,0,0,0,0,0,0,0,0,0,0 /

C      data rpstr(1) /'Ir (a,b,c <-> z,x,y)   '/
C      data rpstr(2) /'IIr (a,b,c <-> y,z,x)  '/
C      data rpstr(3) /'IIIr (a,b,c <-> x,y,z) '/					 
      end

