!---------------------------------------------------------------------!
      module basis
!---------------------------------------------------------------------!
      contains
!---------------------------------------------------------------------!
      subroutine initbase
      use rundata
      implicit none
      type(maxbase) b
      integer::i,j,k
c      double precision::con
c      con=(2/pi)**(3./4.)
      k=0
      do i=1,nion
        do j=1,nint(zion(i))
         if(j.eq.1.or.j.eq.2)then
          no=1
          b=basisfunction(ng,no,0)
          sto%a(k+1:k+ng)=b%a(1:ng); sto%d(k+1:k+ng)=b%d(1:ng)
          sto%i(k+1:k+ng)=i; sto%c(k+1:k+ng,:)=b%c(1:ng,:)
          k=k+ng
         elseif(j.eq.3.or.j.eq.4)then
          no=2
          b=basisfunction(ng,no,0)
          sto%a(k+1:k+ng)=b%a(1:ng); sto%d(k+1:k+ng)=b%d(1:ng)
          sto%i(k+1:k+ng)=i; sto%c(k+1:k+ng,:)=b%c(1:ng,:)
          k=k+ng
         elseif(j.eq.5.or.j.eq.6)then
          no=3
          b=basisfunction(ng,no,1)
          sto%a(k+1:k+ng)=b%a(1:ng); sto%d(k+1:k+ng)=b%d(1:ng)
          sto%i(k+1:k+ng)=i; sto%c(k+1:k+ng,:)=b%c(1:ng,:)
          k=k+ng
         elseif(j.eq.7.or.j.eq.8)then
          no=3
          b=basisfunction(ng,no,2)
          sto%a(k+1:k+ng)=b%a(1:ng); sto%d(k+1:k+ng)=b%d(1:ng)
          sto%i(k+1:k+ng)=i; sto%c(k+1:k+ng,:)=b%c(1:ng,:)
          k=k+ng
         elseif(j.eq.9.or.j.eq.10)then
          no=3
          b=basisfunction(ng,no,3)
          sto%a(k+1:k+ng)=b%a(1:ng); sto%d(k+1:k+ng)=b%d(1:ng)
          sto%i(k+1:k+ng)=i; sto%c(k+1:k+ng,:)=b%c(1:ng,:)
          k=k+ng
         endif
        enddo
      enddo
      if(allbase) then
       c(1:nb)=sto%d/(sqrt(2.)*ne)
       c(n1:n2)=-c(1:nb)
      else
       allocate(cfix(nb,n0)); cfix=0.d0
       do i=1,ne; j=(i-1)*ng; k=j+ng; j=j+1
        cfix(j:k,i)=sto%d(j:k)!*rrt2!/sqrt(dble(ne))
       enddo
       c=1.e-10; c(1)=1.; if(.not.hartfck)c(n1)=-1.
      endif
      end subroutine initbase
!---------------------------------------------------------------------!
      function basisfunction(ng,no,na) result(sto)
!---------------------------------------------------------------------!
!                                                                     !
!     Basis parameters taken from:                                    !
!                                                                     !
!   - J Hehre, et al., The Journal of Chemical Physics 51, 2657       !
!     (1969)                                                          !
!                                                                     !
!---------------------------------------------------------------------!
      use mytypes
      use rundata, only: pi
      implicit none
      type(maxbase) sto
      integer::no,ng,na
      if(no.eq.1) then
       if(ng.eq.2) then
        sto%a(1)=.151623;  sto%d(1)=.678914;
        sto%a(2)=.851819;  sto%d(2)=.430129;
       elseif(ng.eq.3) then
        sto%a(1)=.109818;  sto%d(1)=.444635;
        sto%a(2)=.405771;  sto%d(2)=.535328;
        sto%a(3)=2.22766;  sto%d(3)=.154329;
       elseif(ng.eq.4) then
        sto%a(1)=.0880187; sto%d(1)=.291626;
        sto%a(2)=.265204;  sto%d(2)=.532846;
        sto%a(3)=.954620;  sto%d(3)=.260141;
        sto%a(4)=5.21686;  sto%d(4)=.0567523;
       elseif(ng.eq.5) then
        sto%a(1)=.0744527; sto%d(1)=.193572;
        sto%a(2)=.197572;  sto%d(2)=.482570;
        sto%a(3)=.578648;  sto%d(3)=.331816;
        sto%a(4)=2.07173;  sto%d(4)=.113541;
        sto%a(5)=11.3056;  sto%d(5)=.0221406;
       elseif(ng.eq.6) then
        sto%a(1)=.0651095; sto%d(1)=.130334;
        sto%a(2)=.158088;  sto%d(2)=.416492;
        sto%a(3)=.407099;  sto%d(3)=.370563;
        sto%a(4)=1.18506;  sto%d(4)=.168538;
        sto%a(5)=4.23592;  sto%d(5)=.0493615;
        sto%a(6)=23.1030;  sto%d(6)=.00916360;
       endif
       sto%c=0
       sto%d=sto%d*(2*sto%a/pi)**(3./4.)
      elseif(no.eq.2) then
       if(ng.eq.2) then
        sto%a(1)=.0974545; sto%d(1)=.963782;
        sto%a(2)=.384244;  sto%d(2)=.0494718;
       elseif(ng.eq.3) then
        sto%a(1)=.0751386; sto%d(1)=.700115;
        sto%a(2)=.231031;  sto%d(2)=.399513;
        sto%a(3)=.994203;  sto%d(3)=-.0999672;
       elseif(ng.eq.4) then
        sto%a(1)=.0628104; sto%d(1)=.497767;
        sto%a(2)=.163541;  sto%d(2)=.558855;
        sto%a(3)=.502989;  sto%d(3)=.000029768;
        sto%a(4)=2.32350;  sto%d(4)=-.0622071;
       elseif(ng.eq.5) then
        sto%a(1)=.0544949; sto%d(1)=.346121;
        sto%a(2)=.127920;  sto%d(2)=.612290;
        sto%a(3)=.329060;  sto%d(3)=.128997;
        sto%a(4)=1.03250;  sto%d(4)=-.0653275;
        sto%a(5)=5.03629;  sto%d(5)=-.0294086;
       elseif(ng.eq.6) then
        sto%a(1)=.0485690; sto%d(1)=.240706;
        sto%a(2)=.105960;  sto%d(2)=.595117;
        sto%a(3)=.243977;  sto%d(3)=.250242;
        sto%a(4)=.634142;  sto%d(4)=-.0337854;
        sto%a(5)=2.04036;  sto%d(5)=-.0469917;
        sto%a(6)=.103087;  sto%d(6)=-.0132528;
       endif
       sto%c=0
       sto%d=sto%d*(2*sto%a/pi)**(3./4.)
      elseif(no.eq.3) then
       if(ng.eq.2) then
        sto%a(1)=.0974545; sto%d(1)=.612820;
        sto%a(2)=.384244;  sto%d(2)=.511541;
       elseif(ng.eq.3) then
        sto%a(1)=.0751386; sto%d(1)=.391957;
        sto%a(2)=.231031;  sto%d(2)=.607684;
        sto%a(3)=.994203;  sto%d(3)=.155916;
       elseif(ng.eq.4) then
        sto%a(1)=.0628104; sto%d(1)=.246313;
        sto%a(2)=.163541;  sto%d(2)=.583575;
        sto%a(3)=.502989;  sto%d(3)=.286379;
        sto%a(4)=2.32350;  sto%d(4)=.0436843;
       elseif(ng.eq.5) then
        sto%a(1)=.0544949; sto%d(1)=.156828;
        sto%a(2)=.127920;  sto%d(2)=.510240;
        sto%a(3)=.329060;  sto%d(3)=.373598;
        sto%a(4)=1.03250;  sto%d(4)=.107558;
        sto%a(5)=5.03629;  sto%d(5)=.0125561;
       elseif(ng.eq.6) then
        sto%a(1)=.0485690; sto%d(1)=.101708;
        sto%a(2)=.105960;  sto%d(2)=.425860;
        sto%a(3)=.243977;  sto%d(3)=.418036;
        sto%a(4)=.634142;  sto%d(4)=.173897;
        sto%a(5)=2.04036;  sto%d(5)=.0376794;
        sto%a(6)=.103087;  sto%d(6)=.0037597;
       endif
       sto%c=0; sto%c(:,na)=1
       sto%d=sto%d*(128*sto%a**5/pi**3)**(1./4.)
      else
       stop"Invalid value for n_orbital"
      endif
      end function basisfunction
!---------------------------------------------------------------------!
!                                                                     !
!     Routines to calculate the overlap matrix, S, and the            !
!      transformation matrix, X, which diagonalizes S. We use the     !
!      recurrance relations presented in:                             !
!                                                                     !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals I. Overlap Integrals. The          !
!     Mathematica Journal, 14, 2012.                                  !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine smatrix
      use mytypes
      use rundata
      use functions, only: transform
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,id,jd,i,j 
      double precision::sij
      call cpu_time(cputime(1))
      if(allbase) then
       do i=1,n0
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=i,n0
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g3%expo=g1%expo+g2%expo
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)/g3%expo
         g3%coef=exp(-g1%expo*g2%expo*sum((g1%cent-g2%cent)**2)/g3%expo)
         smat(i,j)=sqrt(pi/g3%expo)**3*g3%coef*
     .             olap(g1%cent(1),g1%cart(1),g1%expo,
     .                   g2%cent(1),g2%cart(1),g2%expo,g3%cent(1))*
     .               olap(g1%cent(2),g1%cart(2),g1%expo,
     .                     g2%cent(2),g2%cart(2),g2%expo,g3%cent(2))*
     .                 olap(g1%cent(3),g1%cart(3),g1%expo,
     .                       g2%cent(3),g2%cart(3),g2%expo,g3%cent(3))
         if(i.ne.j)smat(j,i)=smat(i,j)
        enddo
       enddo
      else
       smat=0.d0
       do i=1,nb
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=1,nb
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g3%expo=g1%expo+g2%expo
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)/g3%expo
         g3%coef=exp(-g1%expo*g2%expo*sum((g1%cent-g2%cent)**2)/g3%expo)
         sij=sqrt(pi/g3%expo)**3*g3%coef*
     .             olap(g1%cent(1),g1%cart(1),g1%expo,
     .                   g2%cent(1),g2%cart(1),g2%expo,g3%cent(1))*
     .               olap(g1%cent(2),g1%cart(2),g1%expo,
     .                     g2%cent(2),g2%cart(2),g2%expo,g3%cent(2))*
     .                 olap(g1%cent(3),g1%cart(3),g1%expo,
     .                       g2%cent(3),g2%cart(3),g2%expo,g3%cent(3))
         do ie=1,n0; do je=1,n0
          smat(ie,je)=smat(ie,je)+cfix(i,ie)*cfix(j,je)*sij
         enddo; enddo
        enddo
       enddo
      endif
      if(hartfck)then
       xmat=smat
      else
       xmat(1:n0,1:n0)=smat; xmat(n1:n2,n1:n2)=xmat(1:n0,1:n0)
       xmat(1:n0,n1:n2)=0.d0; xmat(n1:n2,1:n0)=0.d0
      endif
      call transform(nh,xmat,nx)
      end subroutine smatrix
!---------------------------------------------------------------------!
      recursive function olap(ba,la,ga,bb,lb,gb,bp) result(sx)
      implicit none
      integer::la,lb
      double precision::ba,ga,bb,gb,bp,sx
      intent(in) ba,la,ga,bb,lb,gb,bp
      if(lb.eq.0) then
       if(la.eq.0) then
        sx=1
       elseif(la.eq.1) then
        sx=bp-ba
       else
        sx=(bp-ba)*olap(ba,la-1,ga,bb,lb,gb,bp)+
     .           (la-1)*olap(ba,la-2,ga,bb,lb,gb,bp)/(2*(ga+gb))
       endif
      else
       sx=olap(ba,la+1,ga,bb,lb-1,gb,bp)+
     .     (ba-bb)*olap(ba,la,ga,bb,lb-1,gb,bp)
      endif
      end function olap
!---------------------------------------------------------------------!
!                                                                     !
!     Routines to calculate the kinetic matrix, K. We use the         !
!      recurrance relations presented in:                             !
!                                                                     !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals II. Kinetic Energy Integrals. The  !
!     Mathematica Journal, 15, 2013.                                  !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine kmatrix
      use mytypes
      use rundata
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,i,j,k 
      double precision::b,d,kij
      double precision,dimension(3)::s
      call cpu_time(cputime(2))
      if(allbase) then
       do i=1,nb
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=i,nb
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g3%expo=g1%expo+g2%expo
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)/g3%expo
         g3%coef=exp(-g1%expo*g2%expo*sum((g1%cent-g2%cent)**2)/g3%expo)
         do k=1,3
          s(k)=olap(g1%cent(k),g1%cart(k),g1%expo,
     .              g2%cent(k),g2%cart(k),g2%expo,g3%cent(k))
         enddo
         kmat(i,j)=2*g3%coef*sqrt(pi/g3%expo)**3*(
     .       krec(g1%cent(1),g1%cart(1),g1%expo,
     .             g2%cent(1),g2%cart(1),g2%expo,g3%cent(1))*s(2)*s(3)+
     .       krec(g1%cent(2),g1%cart(2),g1%expo,
     .             g2%cent(2),g2%cart(2),g2%expo,g3%cent(2))*s(1)*s(3)+
     .       krec(g1%cent(3),g1%cart(3),g1%expo,
     .             g2%cent(3),g2%cart(3),g2%expo,g3%cent(3))*s(1)*s(2))
         if(i.ne.j)kmat(j,i)=kmat(i,j)
        enddo
       enddo
      else
       kmat=0.d0
       do i=1,nb
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=1,nb
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g3%expo=g1%expo+g2%expo
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)/g3%expo
         g3%coef=exp(-g1%expo*g2%expo*sum((g1%cent-g2%cent)**2)/g3%expo)
         do k=1,3
          s(k)=olap(g1%cent(k),g1%cart(k),g1%expo,
     .              g2%cent(k),g2%cart(k),g2%expo,g3%cent(k))
         enddo
         kij=2*g3%coef*sqrt(pi/g3%expo)**3*(
     .       krec(g1%cent(1),g1%cart(1),g1%expo,
     .             g2%cent(1),g2%cart(1),g2%expo,g3%cent(1))*s(2)*s(3)+
     .       krec(g1%cent(2),g1%cart(2),g1%expo,
     .             g2%cent(2),g2%cart(2),g2%expo,g3%cent(2))*s(1)*s(3)+
     .       krec(g1%cent(3),g1%cart(3),g1%expo,
     .             g2%cent(3),g2%cart(3),g2%expo,g3%cent(3))*s(1)*s(2))
         do ie=1,n0; do je=1,n0
          kmat(ie,je)=kmat(ie,je)+cfix(i,ie)*cfix(j,je)*kij
         enddo; enddo
        enddo
       enddo
      endif
      end subroutine kmatrix
!---------------------------------------------------------------------!
      recursive function krec(ba,la,ga,bb,lb,gb,bp) result(kx)
      implicit none
      integer::la,lb
      double precision::ba,ga,bb,gb,bp,kx
      intent(in) ba,la,ga,bb,lb,gb,bp
      if(lb.eq.0) then
       if(la.eq.0) then
        kx=2*ga*gb*olap(ba,1,ga,bb,1,gb,bp)
       else
        kx=-la*gb*olap(ba,la-1,ga,bb,1,gb,bp)+
     .     2*ga*gb*olap(ba,la+1,ga,bb,1,gb,bp)
       endif
      elseif(la.eq.0) then
       kx=-ga*lb*olap(ba,1,ga,bb,lb-1,gb,bp)+
     .    2*ga*gb*olap(ba,1,ga,bb,lb+1,gb,bp)
      else
       kx=0.5*la*lb*olap(ba,la-1,ga,bb,lb-1,gb,bp)-
     .         ga*lb*olap(ba,la+1,ga,bb,lb-1,gb,bp)-
     .          la*gb*olap(ba,la-1,ga,bb,lb+1,gb,bp)+
     .         2*ga*gb*olap(ba,la+1,ga,bb,lb+1,gb,bp)
      endif
      end function krec
!---------------------------------------------------------------------!
!                                                                     !
!     Routines to calculate the nuclear matrix, N. We use the         !
!      recurrance relations presented in:                             !
!                                                                     !
!   - Minhhuy Ho and Julio Manuel Hernandex-Peresz. Evaluation of     !
!     Gaussian Molecular Integrals III. Nuclear-Electron Attraction   !
!     Integrals. The Mathematica Journal, 16, 2014.                   !
!                                                                     !
!     We also use Gauss-Jacobi integration as presented in:           !
!                                                                     !
!   - William H Press, Saul A Teukolsky, William T Vetterling, and    !
!     Brian P Flannery. Numerical recipes in Fortran 77: the art of   !
!     scientific computing, 1992.                                     !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine nmatrix
      use mytypes
      use rundata
      use functions, only: gauss_jacobi
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,i,j,k
      double precision::d,g12expo,nij
      double precision,dimension(ngj)::pl
      call cpu_time(cputime(3))
      call gauss_jacobi(-.5d0,-.5d0)
      nmat=0.d0
      if(allbase) then
       do i=1,nb
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=i,nb
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g12expo=1.d0/(g1%expo+g2%expo)
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)*g12expo
         g3%expo=g1%expo*g2%expo*g12expo
         g3%coef=-exp(-g3%expo*sum((g1%cent-g2%cent)**2))
         do k=1,ngj; pl(k)=npoly(g1,g2,xgj(k)); enddo
         nmat(i,j)=g3%coef*2*pi*sum(wgj*pl)*g12expo
         if(i.ne.j)nmat(j,i)=nmat(i,j)
        enddo
       enddo
      else
       do i=1,nb
        g1%expo=sto%a(i)
        g1%cent=ion(sto%i(i),:,iani)
        g1%cart=sto%c(i,:)
        do j=1,nb
         g2%expo=sto%a(j)
         g2%cent=ion(sto%i(j),:,iani)
         g2%cart=sto%c(j,:)
         g12expo=1.d0/(g1%expo+g2%expo)
         g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)*g12expo
         g3%expo=g1%expo*g2%expo*g12expo
         g3%coef=-exp(-g3%expo*sum((g1%cent-g2%cent)**2))
         do k=1,ngj; pl(k)=npoly(g1,g2,xgj(k)); enddo
         nij=g3%coef*2*pi*sum(wgj*pl)*g12expo
         do ie=1,n0; do je=1,n0
          nmat(ie,je)=nmat(ie,je)+cfix(i,ie)*cfix(j,je)*nij
         enddo; enddo
        enddo
       enddo
      endif
      end subroutine nmatrix
!---------------------------------------------------------------------!
      function npoly(g1,g2,t) result(np)
      use mytypes
      use rundata, only:ion,nion,zion,iani,pi
      implicit none
      type(gauss) g1,g2
      integer::i
      double precision::t,nx,ny,nz,np,p
      double precision,dimension(3)::bp,br
      p=g1%expo+g2%expo
      bp=(g1%expo*g1%cent+g2%expo*g2%cent)/p
      np=0
      do i=1,nion
       br=ion(i,:,iani)
       nx=nrec(t,g1%cent(1),g1%cart(1),g1%expo,
     .            g2%cent(1),g2%cart(1),g2%expo,bp(1),br(1))
       ny=nrec(t,g1%cent(2),g1%cart(2),g1%expo,
     .            g2%cent(2),g2%cart(2),g2%expo,bp(2),br(2))
       nz=nrec(t,g1%cent(3),g1%cart(3),g1%expo,
     .            g2%cent(3),g2%cart(3),g2%expo,bp(3),br(3))
       np=np+zion(i)*nx*ny*nz*exp(-p*t**2*sum((bp-br)**2))*sqrt(1-t**2)
      enddo
      end function npoly
!---------------------------------------------------------------------!
      recursive function nrec(t,ba,la,ga,bb,lb,gb,bp,br) result(nx)
      implicit none
      integer::la,lb
      double precision::t,ba,ga,bb,gb,bp,br,kx,nx
      if(lb.eq.0) then  
       if(la.eq.0) then
        nx=1
       elseif(la.eq.1) then
        nx=-(ba-bp+t**2*(bp-br))
       else
        nx=-(ba-bp+t**2*(bp-br))*nrec(t,ba,la-1,ga,bb,lb,gb,bp,br)+
     .      (la-1)*(1-t**2)*nrec(t,ba,la-2,ga,bb,lb,gb,bp,br)/(2*(ga+gb))
       endif
      else
       nx=nrec(t,ba,la+1,ga,bb,lb-1,gb,bp,br)+
     .                    (ba-bb)*nrec(t,ba,la,ga,bb,lb-1,gb,bp,br)
      endif
      end function nrec
!---------------------------------------------------------------------!
!                                                                     !
!     Routines to calculate the coulomb matrix, Q. We use the         !
!      recurance relations and auxilary integral moethod presented    !
!      in:                                                            !
!                                                                     !
!   - Shigeru Obara and A Saika. Efficient recursive computation of   !
!     molecular integrals over Cartesian Gaussian functions. The      !
!     Journal of chemical physics, 84(7):3963–3974, 1986.             !
!                                                                     !
!     We exploit the symmetry of the Q-matrix to avoid calculating    !
!      repeated parts. Namely,                                        !
!                                                                     !
!      (ab|dc)=(ab|cd) "qmat(ia,ib,id,ic)=qmat(ia,ib,ic,id)"          !
!      (ba|dc)=(ab|cd) "qmat(ib,ia,:,:)=qmat(ia,ib,:,:)"              !
!      (cd|ab)=(ab|cd) not implemented yet...                         !
!                                                                     !
!     These symmetry properties are given in:                         !
!                                                                     !
!     - David Yarkony. Modern electronic structure theory, volume 2.  !
!       World Scientific, 1995                                        !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine qmatrix
      use rundata, only: n0,nb,ng,ion,sto,qmat,allbase,cfix,iani,cputime
      use mytypes
      implicit none
      type(nrecval) pq
      integer::ae,be,ce,de,ia,ib,ic,id
      double precision::qabcd
      call cpu_time(cputime(4))
      if(allbase) then
       do ia=1,n0
        pq%ga=sto%a(ia); pq%RA=ion(sto%i(ia),:,iani); pq%a=sto%c(ia,:)
        do ib=ia,n0
         pq%gb=sto%a(ib); pq%RB=ion(sto%i(ib),:,iani); pq%b=sto%c(ib,:)
         pq%gp=pq%ga+pq%gb; pq%RP=(pq%ga*pq%RA+pq%gb*pq%RB)/pq%gp
         do ic=1,n0
          pq%gc=sto%a(ic); pq%RC=ion(sto%i(ic),:,iani); pq%c=sto%c(ic,:)
          do id=ic,n0
           pq%gd=sto%a(id); pq%RD=ion(sto%i(id),:,iani); pq%d=sto%c(id,:)
           pq%gq=pq%gc+pq%gd; pq%RQ=(pq%gc*pq%RC+pq%gd*pq%RD)/pq%gq
           pq%gw=pq%gp+pq%gq; pq%RW=(pq%gp*pq%RP+pq%gq*pq%RQ)/pq%gw
           qmat(ia,ib,ic,id)=qrec(pq,0)
           if(id.ne.ic) qmat(ia,ib,id,ic)=qmat(ia,ib,ic,id)
          enddo
         enddo
         if(ib.ne.ia) qmat(ib,ia,:,:)=qmat(ia,ib,:,:)
        enddo
       enddo
      else
       qmat=0.d0
       do ia=1,nb
        pq%ga=sto%a(ia); pq%RA=ion(sto%i(ia),:,iani); pq%a=sto%c(ia,:)
        do ib=ia,nb
         pq%gb=sto%a(ib); pq%RB=ion(sto%i(ib),:,iani); pq%b=sto%c(ib,:)
         pq%gp=pq%ga+pq%gb; pq%RP=(pq%ga*pq%RA+pq%gb*pq%RB)/pq%gp
         do ic=1,nb
          pq%gc=sto%a(ic); pq%RC=ion(sto%i(ic),:,iani); pq%c=sto%c(ic,:)
          do id=ic,nb
           pq%gd=sto%a(id); pq%RD=ion(sto%i(id),:,iani); pq%d=sto%c(id,:)
           pq%gq=pq%gc+pq%gd; pq%RQ=(pq%gc*pq%RC+pq%gd*pq%RD)/pq%gq
           pq%gw=pq%gp+pq%gq; pq%RW=(pq%gp*pq%RP+pq%gq*pq%RQ)/pq%gw
           qabcd=qrec(pq,0)
           if(id.ne.ic.and.ib.ne.ia) then
            do ae=1,n0; do be=1,n0; do ce=1,n0; do de=1,n0
             qmat(ae,be,ce,de)=qmat(ae,be,ce,de)+4*cfix(ia,ae)*
     .                        cfix(ib,be)*cfix(ic,ce)*cfix(id,de)*qabcd
            enddo; enddo; enddo; enddo 
           elseif(id.ne.ic.or.ib.ne.ia) then
            do ae=1,n0; do be=1,n0; do ce=1,n0; do de=1,n0
             qmat(ae,be,ce,de)=qmat(ae,be,ce,de)+2*cfix(ia,ae)*
     .                        cfix(ib,be)*cfix(ic,ce)*cfix(id,de)*qabcd
            enddo; enddo; enddo; enddo 
           else
            do ae=1,n0; do be=1,n0; do ce=1,n0; do de=1,n0
             qmat(ae,be,ce,de)=qmat(ae,be,ce,de)+cfix(ia,ae)*
     .                        cfix(ib,be)*cfix(ic,ce)*cfix(id,de)*qabcd
            enddo; enddo; enddo; enddo
           endif
          enddo
         enddo
        enddo
       enddo 
      endif
      end subroutine qmatrix
!---------------------------------------------------------------------!
      recursive function qrec(pq,m) result(Q)
      use rundata, only: pi52,ngj,xgj,wgj
      use mytypes
      implicit none
      type(nrecval) pq,m1,m2,p1,mm
      integer::i,j,m
      integer,dimension(3)::r1,r2,r3,r4
      double precision::T,Eab,Ecd,Del,Q
      double precision,dimension(ngj)::pl
      intent(in)pq,m
      if(pq%d(1).eq.0.and.pq%d(2).eq.0.and.pq%d(3).eq.0) then!    all-d
       if(pq%b(1).eq.0.and.pq%b(2).eq.0.and.pq%b(3).eq.0) then!   all-b
        if(pq%c(1).eq.0.and.pq%c(2).eq.0.and.pq%c(3).eq.0) then!  all-c
         if(pq%a(1).eq.0.and.pq%a(2).eq.0.and.pq%a(3).eq.0) then! all-a
          T=sum((pq%RP-pq%RQ)**2)*pq%gp*pq%gq/pq%gw
          Eab=exp(-pq%ga*pq%gb*sum((pq%RA-pq%RB)**2)/pq%gp)
          Ecd=exp(-pq%gc*pq%gd*sum((pq%RC-pq%RD)**2)/pq%gq)
          do i=1,ngj; pl(i)=qpoly(m,T,xgj(i)); enddo
          Q=2.*pi52*Eab*Ecd*sum(wgj*pl)/(pq%gp*pq%gq*sqrt(pq%gw))
         else!                                                    all-a
          do i=1,3; if(pq%a(i).gt.0) j=i; enddo
          m1=pq; m1%a(j)=m1%a(j)-1
          if(pq%a(j).gt.1) then!                                  big-a
           m2=pq; m2%a(j)=m2%a(j)-2
           Q=(pq%RP(j)-pq%RA(j))*qrec(m1,m)+
     .        (pq%RW(j)-pq%RP(j))*qrec(m1,m+1)+
     .         (pq%a(j)-1)*(qrec(m2,m)-
     .          .5*pq%gq*qrec(m2,m+1)/pq%gw)/(2*pq%gp)
          else!                                                   big-a
           Q=(pq%RP(j)-pq%RA(j))*qrec(m1,m)+
     .        (pq%RW(j)-pq%RP(j))*qrec(m1,m+1)
          endif
         endif!                                                   all-a
        else!                                                     all-c
         do i=1,3; if(pq%c(i).gt.0) j=i; enddo
         Del=pq%gb*(pq%RB(j)-pq%RA(j))+pq%gd*(pq%RD(j)-pq%RC(j))
         m1=pq; m1%c(j)=m1%c(j)-1
         p1=m1; p1%a(j)=p1%a(j)+1
         if(pq%c(j).gt.1) then!                                   big-c
          m2=pq; m2%c(j)=m2%c(j)-2
          if(pq%a(j).gt.0) then!                                c-and-a
           mm=m1; mm%a(j)=mm%a(j)-1
           Q=(Del*qrec(m1,m)-pq%gp*qrec(p1,m)+
     .          .5*(pq%a(j)*qrec(mm,m)+(pq%c(j)-1)*qrec(m2,m)))/pq%gq
          else!                                                 c-and-a
           Q=(Del*qrec(m1,m)-pq%gp*qrec(p1,m)+
     .          .5*(pq%c(j)-1)*qrec(m2,m))/pq%gq
          endif!                                                c-and-a
         else!                                                    big-c
          if(pq%a(j).gt.0) then!                                c-and-a
           mm=m1; mm%a(j)=mm%a(j)-1
           Q=(Del*qrec(m1,m)-pq%gp*qrec(p1,m)+
     .          .5*pq%a(j)*qrec(mm,m))/pq%gq
          else!                                                 c-and-a
           Q=(Del*qrec(m1,m)-pq%gp*qrec(p1,m))/pq%gq
          endif!                                                c-and-a
         endif!                                                   big-c
        endif!                                                    all-c
       else!                                                      all-b
        do i=1,3; if(pq%b(i).gt.0) j=i; enddo
         m1=pq; m1%b(j)=m1%b(j)-1
         p1=m1; p1%a(j)=p1%a(j)+1
         Q=qrec(p1,m)+(pq%RA(j)-pq%RB(j))*qrec(m1,m)
       endif!                                                     all-b
      else!                                                       all-d
       do i=1,3; if(pq%d(i).gt.0) j=i; enddo
        m1=pq; m1%d(j)=m1%d(j)-1
        p1=m1; p1%c(j)=p1%c(j)+1
        Q=qrec(p1,m)+(pq%RC(j)-pq%RD(j))*qrec(m1,m)
      endif!                                                      all-d
      end function qrec
!---------------------------------------------------------------------!
      function qpoly(m,T,x) result(q)
      implicit none
      integer::m
      double precision::q,x,T
      q=sqrt(1-x**2)*x**(2*m)*exp(-T*x**2)
      end function qpoly
!---------------------------------------------------------------------!
!                                                                     !
!     Routines to convert from the gaussian basis into a real-space   !
!      grid                                                           !
!                                                                     !
!---------------------------------------------------------------------!
      function plotbase(c) result(den)
      use rundata, only: xmin,dx,ion,sto,n,nh,n0,ng,allbase,iani
      implicit none
      integer::i,j,k,g,a,ag
      double precision::absr,ang(3),amom
      double precision,dimension(nh)::c
      double precision,dimension(n(1),n(2),n(3))::den
      double precision,dimension(3)::r,rH
      rH=0.
      if(allbase)then
       do i=1,n(1); r(1)=dx(1)*i+xmin
        do j=1,n(2); r(2)=dx(2)*j+xmin
         do k=1,n(3); r(3)=dx(3)*k+xmin
          den(i,j,k)=0.
          do g=1,n0
           absr=sum((r-ion(sto%i(g),:,iani))**2)
           ang=(r-ion(sto%i(g),:,iani))**(sto%c(g,:))
           amom=ang(1)*ang(2)*ang(3)
           den(i,j,k)=den(i,j,k)+(amom*c(g)*exp(-sto%a(g)*absr))**2
          enddo
         enddo
        enddo
       enddo
      else
       do i=1,n(1); r(1)=dx(1)*i+xmin
        do j=1,n(2); r(2)=dx(2)*j+xmin
         do k=1,n(3); r(3)=dx(3)*k+xmin
          den(i,j,k)=0.
          ag=0
          do a=1,n0
           do g=1,ng
            ag=ag+1
            absr=sum((r-ion(sto%i(ag),:,iani))**2)
            ang=(r-ion(sto%i(ag),:,iani))**(sto%c(ag,:))
            amom=ang(1)*ang(2)*ang(3)
            den(i,j,k)=den(i,j,k)+(amom*sto%d(ag)*exp(-sto%a(ag)*absr))**2
           enddo
           den(i,j,k)=den(i,j,k)*c(a)
          enddo
         enddo
        enddo
       enddo
      endif
      end function plotbase
!---------------------------------------------------------------------!
      end module basis
!---------------------------------------------------------------------!
