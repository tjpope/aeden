!---------------------------------------------------------------------!
      module io
      contains
!---------------------------------------------------------------------!
      subroutine moses
      use rundata
      implicit none
      integer::i
      character(100),allocatable,dimension(:)::comlin
      ng=6; n=100; xmax=10.; xmin=-10.
      nion=-1; syslab=""; np=50; ngj=50;
      verbose=.false.; allbase=.true.; hartfck=.false.
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
       elseif (index(comlin(i),'-cb').gt.0) then
        allbase=.false.
       elseif (index(comlin(i),'-hf').gt.0) then
        hartfck=.true.
       elseif (index(comlin(i),'-v').gt.0) then
        verbose=.true.
       elseif (index(comlin(i),'--help').gt.0) then
        write(*,100); write(*,200); stop
       endif
      enddo
      if(trim(syslab).eq."")stop'being a cock! you need provide a system label'
      if(ng.lt.2.or.ng.gt.6)stop'being a cock! unsupported basis'
      if(maxval(n).gt.1000)stop'being a cock! n is too big'
      dx=(xmax-xmin)/n
      write(*,100)
      write(*,300)ng,xmin,xmax,n
100   format("###############################################################",/,"#",t63,"#",/,
     .       "#       AA       EEEEEEEE   DDDDDDDD     EEEEEEEE  NN    NN   #",/,
     .       "#      AAAA      EE         DD      DD   EE        NNN   NN   #",/,
     .       "#     AA  AA     EE         DD      DD   EE        NNNN  NN   #",/,
     .       "#    AA    AA    EEEEEEEE   DD      DD   EEEEEEEE  NN NN NN   #",/,
     .       "#   AAAAAAAAAA   EE         DD      DD   EE        NN  NNNN   #",/,
     .       "#   AA      AA   EE         DD      DD   EE        NN   NNN   #",/,
     .       "#   AA      AA   EEEEEEEE   DDDDDDDD     EEEEEEEE  NN    NN   #",/,"#",t63,"#",/,
     .       "###############################################################")
200   format("#",t63,"#",/,"#",4x,"Usage:",t63,"#",/,"#",t63,"#",/,
     .       "#",4x,"./aeden.x -lab XXXX {options}",t63,"#",/,"#",t63,"#",/,
     .       "#",4x,"Essential Input",t63,"#",/,
     .       "#",4x,"-lab XXXX",t20,"Input label - XXXX.xyz must exist",t63,"#",/,
     .       "#",t63,"#",/,"#",4x,"Optional Input",t63,"#",/,
     .       "#",4x,"-hf",t20,"Run Hartree-Fock instead of Extended",
     .       t63,"#",/,"#",t20,"Electrons",t63,"#",/,
     .       "#",4x,"-cb",t20,"Contract Basis before SCF cycle",t63,"#",/,
     .       "#",4x,"-ng XX",t20,"Number of gaussian in STO-nG basis - XX",
     .       t63,"#",/,"#",t20,"must be an integer",t63,"#",/,"#",t63,"#",/,
     .       "###############################################################")
300   format("#",t63,"#",/,"#",4x,"Using the STO-",i1,"G Basis Set.",t63,"#",/,
     .       "#",4x,"Cell Dimensions:",t63,"#",/,
     .       "#",5x,"xmax = ",f10.5,t63,"#",/,
     .       "#",5x,"xmin = ",f10.5,t63,"#",/,
     .       "#",5x,"n = ",3(i0,2x),t63,"#",/,"#",t63,"#")
      end subroutine moses
!---------------------------------------------------------------------!
      subroutine input
      use rundata
      implicit none
      integer::i
      character(100)::fname
      character(3)::atom
      logical::ex
      write(fname,'(a,".xyz")') trim(syslab)
      inquire(file=fname,exist=ex)
      if(.not.ex)stop'being a cock! input file not there'
      open(60,file=fname)
      read(60,*) nion
      read(60,*)
      allocate(ion(nion,3),zion(nion),slat(nion))
      do i=1,nion
       read(60,*)atom,ion(i,:)
       zion(i)=atomcharge(atom)       
      enddo
      close(60)
      ne=sum(zion); nb=ne*ng;
      if(allbase)then; n0=nb; else; n0=ne; endif
      n1=n0+1; n2=n0*2
      if(hartfck)then; nh=n0; else; nh=n2; endif
      if(ne>91)stop'being a cock! System too big for this Q matrix'
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
      use rundata, only: smat,kmat,nmat,qmat,c,mu,hartfck
      implicit none
      double precision::pop(3),etot
      pop(1:2)=pops(c); pop(3)=sum(pop(1:2))
      if(hartfck)then
       write(*,100)pop(1)
      else
       write(*,200)pop
      endif
      call finalenergy(kmat,nmat,qmat,c)
100   format("#",4x,"Population       =",f10.5,t63,"#",/)
200   format("#",4x,"Mass Population  =",f10.5,t63,"#",/,
     .       "#",4x,"Spin Population  =",f10.5,t63,"#",/,
     .       "#",4x,"Total Population =",f10.5,t63,"#",/,"#",t63,"#")
      end subroutine energies
!---------------------------------------------------------------------!
      function pops(c) result(pop)
      use rundata, only: n0,n1,n2,nh,smat,hartfck
      implicit none
      integer::i
      double precision,dimension(2)::pop
      double precision,dimension(nh)::c
      double precision,dimension(n0)::m,s
      double precision,dimension(n0,n0)::p
      intent(in)c
      m=c(1:n0); do i=1,n0; p(i,:)=m(i)*m; enddo; pop(1)=sum(p*smat)
      if(.not.hartfck)then
       s=c(n1:n2); do i=1,n0; p(i,:)=s(i)*s; enddo; pop(2)=sum(p*smat)
      endif
      end function pops
!---------------------------------------------------------------------!      
      subroutine finalenergy(T,V,Q,c)
      use scf, only: coulomb,fullden,ionrep
      use rundata, only: etot,ry2ev,n0,n1,n2,nh,ne,hartfck
      implicit none
      double precision::EK,EN,EH,Eion
      double precision,dimension(nh)::c
      double precision,dimension(n0,n0)::T,V,VH,RHO
      double precision,dimension(n0,n0,n0,n0)::Q
      double precision,dimension(nh,nh)::F,P
      P=fullden(c); 
      if(hartfck)then
       RHO=P(1:n0,1:n0)
      else
       RHO=P(n1:n2,n1:n2)+P(1:n0,1:n0); 
      endif
      F=0.d0; VH=coulomb(RHO,Q);
      if(hartfck)then
       F(1:n0,1:n0)=VH; EH=sum(P*F)*ry2ev/ne
       F(1:n0,1:n0)=2*T; EK=sum(P*F)*ry2ev/ne
       F(1:n0,1:n0)=2*V; EN=sum(P*F)*ry2ev/ne
      else
       F(1:n0,1:n0)=VH; F(n1:n2,n1:n2)=F(1:n0,1:n0); EH=sum(P*F)*ry2ev/ne
       F(1:n0,1:n0)=2*T; F(n1:n2,n1:n2)=F(1:n0,1:n0); EK=sum(P*F)*ry2ev/ne
       F(1:n0,1:n0)=2*V; F(n1:n2,n1:n2)=F(1:n0,1:n0); EN=sum(P*F)*ry2ev/ne
      endif
      Etot=Etot*ry2ev/sqrt(dble(ne))**2; Eion=ionrep()*ry2ev
      if(hartfck)then
       write(*,100)EN,EK,EH,Eion,etot
      else
       write(*,200)EN,EK,EH,etot-(EN+EH+EK+Eion),Eion,etot
      endif
100   format("#",4x,"Nuclear Energy          = ",f20.5,t63,"#",/,
     .       "#",4x,"Kinetic Energy          = ",f20.5,t63,"#",/,
     .       "#",4x,"Coulomb+Exchange Energy = ",f20.5,t63,"#",/,
     .       "#",4x,"Ion Repulsion           = ",f20.5,t63,"#",/,"#",t63,"#",/,
     .       "#",4x,"Total Energy            = ",f20.5,t63,"#",/,"#",t63,"#",/,
     .       "###############################################################")
200   format("#",4x,"Nuclear Energy = ",f20.5,t63,"#",/,
     .       "#",4x,"Kinetic Energy = ",f20.5,t63,"#",/,
     .       "#",4x,"Coulomb Energy = ",f20.5,t63,"#",/,
     .       "#",4x,"Bivector Term  = ",f20.5,t63,"#",/,
     .       "#",4x,"Ion Repulsion  = ",f20.5,t63,"#",/,"#",t63,"#",/,
     .       "#",4x,"Total Energy   = ",f20.5,t63,"#",/,"#",t63,"#",/,
     .       "###############################################################")
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
      subroutine outmat(fname,n1,n2,dat)
      implicit none
      integer::n1,n2,i,j
      double precision,intent(in),dimension(n1,n2)::dat
      character(*)::fname 
      character(20)::frmt
      intent(in)n1,n2,fname
      write(frmt,'("(",i0,"f7.2)")') n2
      open(60,file=fname)
      do i=1,n1
       write(60,frmt) (dat(i,j),j=1,n2)
      enddo
      close(60)
      end subroutine outmat
!---------------------------------------------------------------------!
      end module io
!---------------------------------------------------------------------!
