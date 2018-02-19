!---------------------------------------------------------------------!
      module io
      contains
!---------------------------------------------------------------------!
      subroutine moses
      use rundata
      implicit none
      integer::i
      character(100),allocatable,dimension(:)::comlin
      ng=3; n=100; xmax=10.; xmin=-10.
      nion=-1; syslab=""; np=50; ngj=50;
      verbose=.false.
      allocate(comlin(command_argument_count()))
      do i=1,command_argument_count()
       call getarg(i,comlin(i))
      enddo
      do i=1,command_argument_count()
       if (index(comlin(i),'-ng').gt.0) then
        read(comlin(i+1),'(i10)') ng
       elseif(index(comlin(i),'-xmax').gt.0) then
        read(comlin(i+1),'(f10.5)') xmax
       elseif(index(comlin(i),'-xmin').gt.0) then
        read(comlin(i+1),'(f10.5)') xmin
       elseif(index(comlin(i),'-lab').gt.0) then
        read(comlin(i+1),'(a)') syslab
       elseif (index(comlin(i),'-v').gt.0) then
        verbose=.true.
       endif
      enddo
      if(trim(syslab).eq."")stop'you need provide a system label'
      if(ng.lt.2.or.ng.gt.6)stop'try again. unsupported basis'
      if(maxval(n).gt.1000)stop'stop being a cock! n is too big'
      dx=(xmax-xmin)/n
      write(*,100)ng,xmin,xmax,n
100   format("#",4x,"Using the STO-",i1,"G Basis Set.",/,
     .       "#",4x,"Cell Dimensions:",/,
     .       "#",5x,"xmax = ",f10.5,/,
     .       "#",5x,"xmin = ",f10.5,/,
     .       "#",5x,"n = ",3(i0,2x),/)
      end subroutine moses
!---------------------------------------------------------------------!
      subroutine input
      use rundata
      implicit none
      integer::i
      character(100)::fname
      character(3)::atom
      write(fname,'(a,".xyz")') trim(syslab)
      open(60,file=fname)
      read(60,*) nion
      read(60,*)
      allocate(ion(nion,3),zion(nion),slat(nion))
      do i=1,nion
       read(60,*)atom,ion(i,:)
       zion(i)=atomcharge(atom)       
      enddo
      close(60)
      ne=sum(zion); nb=ne*ng; n1=nb+1; n2=nb*2
      if(nb>91)stop'System too big for this Q matrix'
      end subroutine input
!---------------------------------------------------------------------!
      function atomcharge(atom) result(chrg)
      implicit none
      integer::i,chrg
      character(3)::atom
      character(3),dimension(118)::labs
      labs=(/ "H ", "He", "Li", "Be", "B ", "C ", "N ", "O ", "F ", 
     .  "Ne", "Na", "Mg", "Al", "Si", "P ", "S ", "Cl", "Ar", "K ", 
     .  "Ca", "Sc", "Ti", "V ", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", 
     .  "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y ", 
     .  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", 
     .  "Sn", "Sb", "Te", "I ", "Xe", "Cs", "Ba", "La", "Ce", "Pr", 
     .  "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", 
     .  "Yb", "Lu", "Hf", "Ta", "W ", "Re", "Os", "Ir", "Pt", "Au", 
     .  "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", 
     .  "Th", "Pa", "U ", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", 
     .  "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", 
     .  "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"/)
      do i=1,118
       if(trim(atom).eq.trim(labs(i))) chrg=i
      enddo
      end function atomcharge
!---------------------------------------------------------------------!
      subroutine energies
      use rundata, only: nb,smat,kmat,nmat,qmat,c,mu
      implicit none
      integer::i,j
      double precision::pop(3),etot
      pop(1:2)=pops(c)
      pop(3)=sum(pop(1:2))
      write(*,100)pop
      call finalenergy(kmat,nmat,qmat,c)
100   format("#",4x,"Mass Population  =",f10.5,/,
     .       "#",4x,"Spin Population  =",f10.5,/,
     .       "#",4x,"Total Population =",f10.5,/)
      end subroutine energies
!---------------------------------------------------------------------!
      function pops(c) result(pop)
      use rundata, only: nb,smat
      implicit none
      integer::i
      double precision,dimension(2)::pop
      double precision,dimension(nb*2)::c,b
      double precision,dimension(nb)::m,s
      double precision,dimension(nb,nb)::p
      intent(in)c
      m=c(1:nb); s=c(nb+1:nb*2)
      do i=1,nb; p(i,:)=m(i)*m; enddo; pop(1)=sum(p*smat)
      do i=1,nb; p(i,:)=s(i)*s; enddo; pop(2)=sum(p*smat)
      end function pops
!---------------------------------------------------------------------!      
      subroutine finalenergy(T,V,Q,c)
      use scf, only: coulomb,fullden,ionrep
      use rundata, only: etot,ry2ev,n1,n2,nb
      implicit none
      double precision::EK,EN,EH,Eion
      double precision,dimension(n2)::c
      double precision,dimension(nb,nb)::T,V,VH,RHO
      double precision,dimension(nb,nb,nb,nb)::Q
      double precision,dimension(n2,n2)::F,P
      P=fullden(c,nb); RHO=P(n1:n2,n1:n2)+P(1:nb,1:nb); F=0.d0
      VH=coulomb(RHO,Q); 
      F(1:nb,1:nb)=VH; F(n1:n2,n1:n2)=F(1:nb,1:nb); EH=sum(P*F)*ry2ev
      F(1:nb,1:nb)=2*T; F(n1:n2,n1:n2)=F(1:nb,1:nb); EK=sum(P*F)*ry2ev
      F(1:nb,1:nb)=2*V; F(n1:n2,n1:n2)=F(1:nb,1:nb); EN=sum(P*F)*ry2ev
      Etot=Etot*ry2ev; Eion=ionrep()*ry2ev
      write(*,100)EN,EK,EH,etot-(EN+EH+EK+Eion),Eion,etot
100   format("#",4x,"Nuclear Energy = ",f20.5,/,
     .       "#",4x,"Kinetic Energy = ",f20.5,/,
     .       "#",4x,"Coulomb Energy = ",f20.5,/,
     .       "#",4x,"Bivector Term  = ",f20.5,/,
     .       "#",4x,"Ion Repulsion  = ",f20.5,/,/,
     .       "#",4x,"Total Energy   = ",f20.5,/)
      end subroutine finalenergy
!---------------------------------------------------------------------!  
      subroutine output(fname,n1,n2,n3,dat,g)
      implicit none
      integer::n1,n2,n3,i,j,k
      double precision,intent(in),optional,dimension(n1)::g
      double precision,intent(in),dimension(n1,n2,n3)::dat
      character(*)::fname 
      intent(in)n1,n2,n3,fname
      open(60,file=fname)
      if (present(g)) then
       if(n3.gt.1) then
        do i=1,n1; do j=1,n2; do k=1,n3
         write(60,*) g(i),g(j),g(k),dat(i,j,k)
        enddo; write(60,*); enddo; write(60,*); enddo
       elseif(n2.gt.1) then
        do i=1,n1; do j=1,n2
         write(60,*) g(i),g(j),dat(i,j,1)
        enddo; write(60,*); enddo
       else
        do i=1,n1; write(60,*) g(i),dat(i,1,1); enddo
       endif
      else
       if(n3.gt.1) then
        do i=1,n1; do j=1,n2; do k=1,n3
         write(60,*) i,j,k,dat(i,j,k)
        enddo; write(60,*); enddo; write(60,*); enddo
       elseif(n2.gt.1) then
        do i=1,n1; do j=1,n2
         write(60,*) i,j,dat(i,j,1)
        enddo; write(60,*); enddo
       else
        do i=1,n1; write(60,*) i,dat(i,1,1); enddo
       endif
      endif
      close(60)
      end subroutine output
!---------------------------------------------------------------------!
      end module io
!---------------------------------------------------------------------!
