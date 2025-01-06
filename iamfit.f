      subroutine lmfit(ndata,npar,nb,DIMFIT,DIMPAR,DIMVB,DIMPLC
     $     ,iprint,nsvfit,nfit,ifit,dfit
     $     ,alpha,covar,evec,beta,w,a,anew,da,freed,palc,pali
     $     ,rofit,sigsq,sigsqold,tol,stepw,fstat,usex,const,alamda
     $     ,fitscl,svderr)
     
C     fit routine Levenberg Marquardt Method
C     Parameters:
C     ndata:  number of data points (transitions)
C     npar:   total number of parameters (fitted and constrained)
C     nb:     number of b's
C     DIMPAR: physical dimension of a,da, ...
C     DIMVB : physical dimension of a,da, ...
C     DIMFIT: physical dimension of alpha,beta, ...
C     DIMPLC: physical dimension of pali, palc
C     iprint: output level 0-5
C     nsvfit: number of indepentent variables
C     nfit:   number of parameters to fit
C     ifit:   array(1:DIMPAR) 1=fit 0=fixed 
C     dfit:   array(1:DIMPAR) 1=fit 0=fixed for each linear comb.
C     alpha:  matrix
C     covar:  freedom \ covariance matrix
C     evec:   eigenvectors of alpha
C     beta:   vector
C     w:      eigenvalues of alpha
C     a:      parameters (old, but good)
C     anew:   parameters (new, just a try)
C     da:     error of paramaters
C     freed:  freedom of parameter
C     palc:   array of linear combinations (no.of.lin.comb, coefficient)  
C     pali:   array of linear combinations (no.of.lin.comb, no.of.parameter or B)  
C     rofit:  robust fitting parameter
C     sigsq:  sum_i (dy_i^2 * weight_i)
C     sigsqold: 
C     tol:    tolerance value to neglect eigenvalues
C     stepw:  stepw: 
C     fstat:  = 0: calc. only chisq  
C             > 0: calc. new parameters (using derivatives dyda)
C             < 0: calc. errors                  
C     usex:   use experimental sigma
C     const:  no. of additional datapoints 
C     alamda: lambda parameter 
C     fitscl: = 0: normalize alpha-Matrix  = 1: don't change alpha
C     svderr: treatment of error-caclulation of nsvfit < nfit
 
      implicit none
C     ..
C     .. Scalar Arguments ..
      real*8 rofit,sigsq,sigsqold,tol,stepw,alamda
      integer ndata,npar,dimfit,nb,dimpar,dimvb,dimplc
      integer iprint,nsvfit,nfit,fstat,usex,const
      integer svderr,fitscl
C     ..
C     .. Array Arguments ..  
      integer ifit(DIMPAR,DIMVB),  dfit(DIMFIT)
      real*8  alpha(DIMFIT,DIMFIT),covar(DIMFIT,DIMFIT),
     $     evec(DIMFIT,DIMFIT)
      real*8  w(DIMFIT),beta(DIMFIT),freed(DIMFIT) 
      real*8  a(DIMPAR,DIMVB),anew(DIMPAR,DIMVB),da(DIMPAR,DIMVB)
      real*8  palc(DIMFIT,-1:DIMPLC)
      integer pali(DIMFIT, 0:DIMPLC,2)
      
C     .. local Scalars
      real*8  chisqex,sig2isum,dy,wt,wmax,wmin,sig2i,thresh,sigex
      real*8  dy2av,psi,psisum,psiav,EPS
      integer i,j,k,l,imin,imax,DIMLMF,b,ic,ir
      parameter (DIMLMF=75,EPS=1.0d-12)!Herbers2024

C     .. local Arrays
      real*8  dyda(DIMLMF),tmp(DIMLMF),dpar(DIMLMF)
      save dy2av

      if (DIMFIT.gt.DIMLMF) stop 'DIMENSION ERROR in iamfit1'
      if (nfit.gt.DIMLMF) stop 'DIMENSION ERROR in iamfit2'
      if (const.ne.0) stop 'FITNEW: const .ne. 0'

      if (fstat.eq.0) then
C     calc sigma**2 
        chisqex=0.0d0
        psisum=0.0
        sig2isum=0.0d0
        do j= 1, ndata+const
          if (j.le.ndata)
     $         call funcs(j,dy,dyda,anew,sigex,nfit,ifit,dfit,1)
C         if (j.gt.ndata)
C     $       call fconst(j-ndata,dy,dyda,anew,sigex,DIMPAR,1,ifit,dfit)
          psi=1.0/(1.0+0.5*((rofit*dy)/sigex)**2)
          call savepsi(j,psi)
          psisum  =psisum   +psi
          chisqex =chisqex  +(dy**2/sigex**2)*psi
          sig2isum=sig2isum +1.0d0/sigex**2
        end do
        if ((ndata-nfit-1).le.0) then
          write(0,*)' FIT: Warning: ndata-nfit = 0: no fit possible!'
          sigsq=100000.0
          dy2av=100000.0
          return
        end if
        psiav=psisum/dble(ndata)
        sigsq=(chisqex/(sig2isum*psiav))
     $       *(dble(ndata)/(dble(ndata-nfit-1)))
        dy2av=chisqex/(dble(ndata-nfit-1))
        
        return
      end if
C ----------------------------------------------------------------------
      if (fstat.gt.0) then
C     Check if new sigma square is better than old sigma square
        if (sigsq.lt.sigsqold) then
C     succsess: calc. new alpha matrix and beta vector
          alamda=alamda*0.4
          do b = 1, nb
            do j = 1, npar
              a(j,b) = anew(j,b)
            end do
          end do
          do j = 1,nfit
            do k = 1,j
              alpha(j,k) = 0.d0
            end do
            beta(j) = 0.d0
          end do
          do i = 1,ndata+const
            do j = 1, nfit
              dyda(j)=0.0d0
            end do
            if (i.le.ndata)
     $           call funcs(i,dy,dyda,anew,sigex,nfit,ifit,dfit,2)
C     if (i.gt.ndata)
C     $         call fconst(i-ndata,dy,dyda,anew,sigex,npar,2,ifit,dfit)
            if (usex.eq.1) then 
              sig2i=1.0d0/(sigex**2)
            else
              sig2i=1.0d0/(sigex**2*dy2av)
            end if
            psi=1.0/(1.0+0.5*((rofit*dy)/sigex)**2)
            sig2i=sig2i*psi
C     build alpha and beta
            do j = 1,nfit
              wt = dyda(j)*sig2i
              do k = 1,j
                alpha(j,k) = alpha(j,k) + wt*dyda(k)
              end do
              beta(j) = beta(j) + dy*wt
            end do
          end do
C     scale beta
          if (fitscl.eq.0) then
            do j = 1,nfit
              beta(j)=beta(j)/dsqrt(alpha(j,j))
            end do
          end if
        else
          alamda=alamda*10.0
        end if

C     put scaled alpha in evec with lambda
        do j = 1,nfit
          if (fitscl.eq.0) then
            do k = 1,j
              evec(j,k) = alpha(j,k)/dsqrt(alpha(j,j)*alpha(k,k))
            end do
          else
            do k = 1,j
              evec(j,k) = alpha(j,k)
            end do
          end if
        end do
        do j = 1,nfit
          evec(j,j) = evec(j,j)*(1.0+alamda) !+alamda
        end do
        call SVDSYDC(evec,w,tmp,nfit,DIMFIT)
        call minmax(w,nfit,wmin,wmax,imin,imax,nsvfit) ! get wmax
        thresh = tol*wmax
        do j= 1,nfit
          if (w(j).lt.thresh) w(j) = 0.D0
        end do  
        call minmax(w,nfit,wmin,wmax,imin,imax,nsvfit) ! get nsvfit
        call svdsybk(evec,w,nfit,DIMFIT,beta,dpar,tmp)

C     rescale
        if (fitscl.eq.0) then
          do j = 1, nfit
            dpar(j)=dpar(j)/dsqrt(alpha(j,j))
          end do
        end if

C     multiply the scale factor palc(i,-1) 
        do j = 1, nfit
          if (palc(j,-1).ne.0.0) dpar(j)=dpar(j)*palc(j,-1)
        end do

C     put the variations of a into da
        do b= 1, nb
          do j = 1, npar
            da(j,b)=0.0d0
          end do
        end do
        do j = 1, nfit
          do i = 1, pali(j,0,1)
            k=pali(j,i,1)
            b=pali(j,i,2)
            da(k,b)=da(k,b)+dpar(j)*palc(j,i)
          end do
        end do
        do b= 1, nb
          do j= 1, npar
            anew(j,b)=a(j,b)+da(j,b)*stepw
          end do
        end do
      end if

C ----------------------------- calc. of errors ----------------
C    don't use lambda in last cylce 
      if (fstat.lt.0) then
        do j = 1,nfit
          if (fitscl.eq.0) then
            do k = 1,j
              evec(j,k) = alpha(j,k)/dsqrt(alpha(j,j)*alpha(k,k))
            end do
          else
            do k = 1,j
              evec(j,k) = alpha(j,k)
            end do
          end if
        end do

C     calculate the covar = inverse of alpha
        call SVDSYDC(evec,w,tmp,nfit,DIMFIT)
        call minmax(w,nfit,wmin,wmax,imin,imax,nsvfit) ! get wmax
        if (svderr.eq.0   ) then
          thresh = EPS*wmax
        else
          thresh = tol*wmax
        end if
        do j= 1,nfit
          if (w(j).lt.thresh) w(j) = 0.D0
        end do
        call minmax(w,nfit,wmin,wmax,imin,imax,nsvfit) ! get nsvfit
        if (svderr.le.0   ) then
          do j= 1,nfit
            if (w(j).eq.0.0d0) w(j) = thresh 
          end do
        end if
        call svdsyinv(evec,w,nfit,DIMFIT,covar)

C     rescale covar
        if (fitscl.eq.0) then
          do j = 1, nfit
            do k = 1, nfit
              covar(j,k)=covar(j,k)/dsqrt(alpha(j,j)*alpha(k,k))
            end do
          end do
        end if

C     da contains the errors of the input parameters
C     dyda is used as a temporary vector        
        do b = 1, nb
          do k = 1, npar
            if (ifit(k,b).le.0) goto 10
            do i = 1, nfit
              dyda(i) = 0.0d0
            end do
            do i = 1, nfit
              do j = 1, pali(i,0,1)  
                if ((pali(i,j,1).eq.k).and.(pali(i,j,2).eq.b)) then
                  dyda(i)=palc(i,j)
                end if
              end do
            end do
            do ir = 1, nfit
              tmp(ir) = 0.0d0
            end do
            do ir = 1, nfit
              do ic = 1, nfit
                tmp(ir)=tmp(ir)+covar(ir,ic)*dyda(ic)
              end do
            end do
            da(k,b)=0.0d0
            do ir = 1, nfit
              da(k,b)=da(k,b)+dyda(ir)*tmp(ir)
            end do
            da(k,b)=dsqrt(da(k,b))
 10         continue
          end do
        end do

C     anew contains now the errors of the fitted linear comb. 
        do j = 1, nfit
          ir=mod(j-1,DIMPAR)+1
          ic=int((j-1)/DIMPAR)+1
          anew(ir,ic)=dsqrt(covar(j,j))
        end do
        do j=1, nfit
            freed(j) = -99.0
        end do

        do i = 1,nfit
C     calc cofreedom matrix on lower left side of covar
          do j = 1,i - 1
            covar(i,j) =
     $           dsqrt(dsqrt((1.0d0-alpha(i,j)**2
     $           / (alpha(i,i)*alpha(j,j)))*
     $           (1.0d0-covar(i,j)**2
     $           / (covar(i,i)*covar(j,j)))))
          end do
C     calc correlation matrix on upper right side of covar
          do j = i + 1, nfit
            covar(i,j) = covar(i,j)/dsqrt(covar(i,i)*covar(j,j))
          end do
        end do
C     calc freedom parameters on diagonal of covar
        do j = 1, nfit
          covar(j,j) = 1.0d0/dsqrt(alpha(j,j)*covar(j,j))
        end do
      end if
      return
      end

C ---------------------------------------------------------------------
      subroutine minmax(w,n,min,max,imin,imax,nw)
      implicit none
      integer n,imin,imax,nw
      real*8  w(n),min,max
      integer i
      min= 1.0d100
      max=-1.0d100
      imin=0
      imax=0
      nw=0
      do i=1, n
        if (w(i).gt.max) then
          max=w(i)
          imax=i
        end if
        if (w(i).ne.0.0d0) then
          nw=nw+1
          if (w(i).lt.min) then
            min=w(i)
            imin=i
          end if
        end if
      end do      
      return
      end

C ---------------------------------------------------------------------
      subroutine SVDSYDC(A,W,TMP,N,DIMA)
C     SVD decomposition of a reell symmetric matrix
      integer dima,n,ierr
      real*8 a(DIMA,DIMA),w(DIMA),tmp(DIMA)
      if (n.gt.dima) stop 'Dimension Error in SVDSYDC'
      call hdiag(dima,n,a,w,tmp,ierr)
      if (ierr.ne.0) then
        write (*,'(a,a,i5)') 'SVDSYDC: Error in hdiag (sydsycd).', 
     $  'Please try to use dqx instead of fit on all parameters',ierr
        write (*,'(60F10.4)') (a(i,i),i=1,n)
        stop
      endif
      call eigsrt(w,a,n,dima)
      return
      end

C ---------------------------------------------------------------------
      subroutine svdsybk(a,w,n,dima,b,x,tmp)
      implicit none
      integer n, dima
      real*8 a(dima,dima), w(n), b(n), x(n), tmp(n)
      integer i,j
C     tmp=A*b
      do i=1, n
        tmp(i)=0.0
        if (w(i).ne.0.0d0) then
          do j=1, n
            tmp(i)=tmp(i)+a(j,i)*b(j)
          end do
          tmp(i)=tmp(i)/w(i)
        end if
      end do
C
      do i=1,n
        x(i)=0.0d0
        do j=1,n
          x(i)=x(i)+a(i,j)*tmp(j)
        end do       
      end do
      return
      end

C ---------------------------------------------------------------------
      subroutine svdsyinv(evec,w,n,dima,ai)
      implicit none
      integer n,dima
      real*8 evec(dima,dima), w(n), ai(dima,dima)
      integer i,j,k
      do i=1,n
        do j=1,n
          ai(i,j)=0.0d0
        end do
        do k=1,n
          if (w(k).ne.0.0d0) then
            do j=1,n
              ai(i,j)=ai(i,j)+evec(i,k)*evec(j,k)/w(k)
            end do
          end if
        end do
      end do
      return
      end

