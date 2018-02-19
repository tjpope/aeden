!---------------------------------------------------------------------!
!                                                                     !
!     Module contains general mathematical functions:                 !
!                                                                     !
!     - Square matrix inverter                                        !
!     - Unitary tranformation matrix finder                           !
!     - Gauss-Jacobi routes and wieghts finder                        !
!     - Log[Gamma[x]] calculator                                      ! 
!                                                                     !
!---------------------------------------------------------------------!
      module functions
!---------------------------------------------------------------------!
      contains
!---------------------------------------------------------------------!
!                                                                     !
!     Function to invert a square matrix. Calls out to LAPACK:        !
!                                                                     !
!   - Edward Anderson, Zhaojun Bai, Christian Bischof, Susan          !
!     Blackford, James Demmel, Jack Dongarra, Jeremy Du Croz, Anne    !
!     Greenbaum, S Hammerling, Alan McKenney, et al. LAPACK Users'    !
!     guide, volume 9. Siam, 1999.                                    !
!                                                                     !
!---------------------------------------------------------------------!
      function invert(a,n) result(b)
      implicit none
      integer::n,info        
      integer,dimension(n)::ipiv
      double precision,dimension(n,n)::a,b
      double precision,dimension(n*n)::work_inv
      intent(in)a,n
      b=a
      call dgetrf(n,n,b,n,ipiv,info)
      call dgetri(n,b,n,ipiv,work_inv,n*n,info)
      end function invert
!---------------------------------------------------------------------!
!                                                                     !
!     Function to find the unitary transformation matrix that         !
!      diagonalizes the input matrix. Calls out to LAPACK:            !
!                                                                     !
!   - Edward Anderson, Zhaojun Bai, Christian Bischof, Susan          !
!     Blackford, James Demmel, Jack Dongarra, Jeremy Du Croz, Anne    !
!     Greenbaum, S Hammerling, Alan McKenney, et al. LAPACK Users'    !
!     guide, volume 9. Siam, 1999.                                    !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine transform(na,M,nb)
      implicit none 
      integer::na,nb,info,i
      double precision,dimension(na,na)::M
      double precision,dimension(na)::e
      double precision,dimension(na*3)::wrk
      intent(in)na
      info=0; call dsyev('V','U',na,M,na,e,wrk,na*3,info)
      nb=0; do i=1,na; if(e(i)>1e-10) then
       nb=nb+1; M(:,nb)=M(:,i)/sqrt(e(i))
      endif; enddo
      end subroutine transform
!---------------------------------------------------------------------!
!                                                                     !
!     Routine to calculate the routes and wieghts for a Gauss-Jacobi  !
!      quadrature integration - copied verbatim from:                 !
!                                                                     !
!   - William H Press, Saul A Teukolsky, William T Vetterling, and    !
!     Brian P Flannery. Numerical recipes in Fortran 77: the art of   !
!     scientific computing, 1992.                                     !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine gauss_jacobi(alf,bet)
      use rundata, only: ngj,xgj,wgj
      implicit none
      integer::m,i,its,j
      double precision::alf,bet,tol,alfbet,an,bn,r1,r2,r3,a,b,c,
     .                                           p1,p2,p3,pp,temp,z,z1
      parameter(tol=3e-14,m=10)
      intent(in)alf,bet
      do i=1,ngj
       if(i.eq.1)then
        an=alf/ngj; bn=bet/ngj
        r1=(1+alf)*(2.78/(4+ngj**2)+.768*an/ngj)
        r2=1+1.48*an+.96*bn+.452*an**2+.83*an*bn
        z=1-r1/r2
       elseif(i.eq.2)then
        r1=(4.1+alf)/((1+alf)*(1+.156*alf))
        r2=1.+.06*(ngj-8.)*(1.+.12*alf)/ngj
        r3=1.+.012*bet*(1.+.25*abs(alf))/ngj
        z=z-(1.-z)*r1*r2*r3
       elseif(i.eq.3)then
        r1=(1.67+.28*alf)/(1.+.37*alf)
        r2=1.+.22*(ngj-8.)/ngj
        r3=1.+8.*bet/((6.28+bet)*ngj*ngj)
        z=z-(xgj(1)-z)*r1*r2*r3
       elseif(i.eq.ngj-1)then
        r1=(1.+.235*bet)/(.766+.119*bet)
        r2=1./(1.+.639*(ngj-4.)/(1.+.71*(ngj-4.)))
        r3=1./(1.+20.*alf/((7.5+alf)*ngj*ngj))
        z=z+(z-xgj(ngj-3))*r1*r2*r3
       elseif(i.eq.ngj)then
        r1=(1.+.37*bet)/(1.67+.28*bet)
        r2=1./(1.+.22*(ngj-8.)/ngj)
        r3=1./(1.+8.*alf/((6.28+alf)*ngj*ngj))
        z=z+(z-xgj(ngj-2))*r1*r2*r3
       else
        z=3.*xgj(i-1)-3.*xgj(i-2)+xgj(i-3)
       endif
       alfbet=alf+bet; its=0; z1=2*z
       do while(abs(z-z1).ge.tol) 
        its=1+its
        temp=2+alfbet; p1=(alf-bet+temp*z)/2; p2=1 
        do j=2,ngj
         p3=p2; p2=p1; temp=2.d0*j+alfbet
         a=2*j*(j+alfbet)*(temp-2); 
         b=(temp-1)*(alf*alf-bet*bet+temp*(temp-2)*z)
         c=2*(j-1+alf)*(j-1+bet)*temp
         p1=(b*p2-c*p3)/a
        enddo
        pp=(ngj*(alf-bet-temp*z)*p1+2*(ngj+alf)*(ngj+bet)*p2)/
     .                                                 (temp*(1-z**2))
        z1=z; z=z1-p1/pp
        if(its.gt.m)stop'Gauss Jacobi has failed'
       enddo
       xgj(i)=z
       wgj(i)=exp(gammln(alf+ngj)+gammln(bet+ngj)-gammln(ngj+1.d0)-
     .                  gammln(ngj+alfbet+1.))*temp*2.**alfbet/(pp*p2)
      enddo
      end subroutine gauss_jacobi
!---------------------------------------------------------------------!
!                                                                     !
!     Routine to calculate the log of the Gamma function - copied     !
!      verbatim from:                                                 !
!                                                                     !
!   - William H Press, Saul A Teukolsky, William T Vetterling, and    !
!     Brian P Flannery. Numerical recipes in Fortran 77: the art of   !
!     scientific computing, 1992.                                     !
!                                                                     !
!---------------------------------------------------------------------!
      function gammln(xx)
      implicit none
      integer::j
      double precision::gammln,xx,ser,stp,tmp,x,y
      double precision,dimension(6)::cof
      cof=(/76.18009172947146d0,-86.50532032941677d0,
     .      24.01409824083091d0,-1.231739572450155d0,
     .     .1208650973866179d-2,-.5395239384953d-5/)
      stp=2.5066282746310005d0; ser=1.000000000190015d0
      x=xx; y=x; tmp=x+5.5d0; tmp=(x+0.5d0)*log(tmp)-tmp
      do j=1,6; y=y+1.d0; ser=ser+cof(j)/y; enddo
      gammln=tmp+log(stp*ser/x)
      return
      end function gammln
!---------------------------------------------------------------------!
      end module functions
!---------------------------------------------------------------------!
