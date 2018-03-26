!---------------------------------------------------------------------!
      module renew_basis
!---------------------------------------------------------------------!
      contains
!---------------------------------------------------------------------!
!                                                                     !
!     Renew the overlap integrals - ignoring the on-site terms, which !
!      don't change                                                   !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine renew_smatrix
      use mytypes
      use rundata
      use functions, only: transform
      use basis
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,ig,jg,i,j 
      double precision::sij
      call cpu_time(cputime(1))
      do ie=1,ne; do ig=1,ng; i=(ie-1)*ng+ig 
       g1%expo=sto%a(i)
       g1%cent=ion(sto%i(i),:,iani)
       g1%cart=sto%c(i,:)
       do je=ie+1,ne; do jg=1,ng; j=(je-1)*ng+jg 
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
        smat(j,i)=smat(i,j)
       enddo; enddo
      enddo; enddo
      if(hartfck)then
       xmat=smat
      else
       xmat(1:n0,1:n0)=smat; xmat(n1:n2,n1:n2)=xmat(1:n0,1:n0)
       xmat(1:n0,n1:n2)=0.d0; xmat(n1:n2,1:n0)=0.d0
      endif
      call transform(nh,xmat,nx)
      end subroutine renew_smatrix
!---------------------------------------------------------------------!
!                                                                     !
!     Renew the kinetic integrals - ignoring the on-site terms, which !
!      don't change                                                   !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine renew_kmatrix
      use mytypes
      use rundata
      use basis
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,ig,jg,i,j,k 
      double precision::b,d,kij
      double precision,dimension(3)::s
      call cpu_time(cputime(2))
      do ie=1,ne; do ig=1,ng; i=(ie-1)*ng+ig
       g1%expo=sto%a(i)
       g1%cent=ion(sto%i(i),:,iani)
       g1%cart=sto%c(i,:)
       do je=ie+1,ne; do jg=1,ng; j=(je-1)*ng+jg
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
        kmat(j,i)=kmat(i,j)
       enddo; enddo
      enddo; enddo
      end subroutine renew_kmatrix
!---------------------------------------------------------------------!
!                                                                     !
!     Renew the nuclear integrals - ignoring the on-site terms, which !
!      don't change                                                   !
!                                                                     !
!---------------------------------------------------------------------!
      subroutine renew_nmatrix
      use mytypes
      use rundata
      use basis
      implicit none
      type(gauss) g1,g2,g3
      integer::ie,je,ig,jg,i,j,k
      double precision::d,g12expo,nij
      double precision,dimension(ngj)::pl
      call cpu_time(cputime(3))
      do ie=1,ne; do ig=1,ng; i=(ie-1)*ng+ig 
       g1%expo=sto%a(i)
       g1%cent=ion(sto%i(i),:,iani)
       g1%cart=sto%c(i,:)
       do je=1,ne; do jg=1,ng; j=(je-1)*ng+jg 
        g2%expo=sto%a(j)
        g2%cent=ion(sto%i(j),:,iani)
        g2%cart=sto%c(j,:)
        g12expo=1.d0/(g1%expo+g2%expo)
        g3%cent=(g1%expo*g1%cent+g2%expo*g2%cent)*g12expo
        g3%expo=g1%expo*g2%expo*g12expo
        g3%coef=-exp(-g3%expo*sum((g1%cent-g2%cent)**2))
        do k=1,ngj; pl(k)=npoly(g1,g2,xgj(k)); enddo
        nmat(i,j)=g3%coef*2*pi*sum(wgj*pl)*g12expo
        nmat(j,i)=nmat(i,j)
       enddo; enddo
      enddo; enddo
      end subroutine renew_nmatrix
!---------------------------------------------------------------------!
      end module renew_basis
!---------------------------------------------------------------------!
