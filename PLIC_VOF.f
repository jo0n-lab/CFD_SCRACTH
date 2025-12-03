
      
      
      
        module param
	        implicit none
	        integer NX,NY,NZ
	        real*8 PI, delta, dt
        endmodule


      
      
      !include 'param.h'
      use param
      
      character*30 tec, ch
      data tec /'out-t0000000000.vtk'/
      
      real*8,allocatable :: f(:,:,:)
      real*8,allocatable :: ux(:,:,:),uy(:,:,:),uz(:,:,:)
      real*8 t, dx,dy,dz,x,y,z,TT
      integer i,j,k,step,tswap

      real*8,allocatable :: fn(:,:,:)
      real*8 tmp,tmp1,init_f
      integer  nnx
      
      PI = 3.141592657589793d0
      nx = 100
      ny = 100
      nz = 100
      dt = 0.001
      delta = 1.0d0/NX
      
      allocate(ux(nx,ny,nz))
      allocate(uy(nx,ny,nz))
      allocate(uz(nx,ny,nz))
      allocate(f(nx,ny,nz))
      allocate(fn(nx,ny,nz))
      ux = 0.0
      uy = 0.0
      uz = 0.0
      f = 0.0
      fn = 0.0
      
      call init(f,ux,uy,uz,t,dx,dy,dz)
      step = 0 
      
      init_f = 0.0d0
      do k=1,NZ
      do j=1,NY
      do i=1,NX
         fn(i,j,k)=f(i,j,k)
         init_f = init_f+f(i,j,k)
      end do
      end do
      end do
      

      call PARAVIEW(tec,NX,NY,NZ,f) 

      
      open (11,file='mass.txt')
      
      TT = 3
      tswap = 0
      t = 0.0

      do step=1,500000 

c ***
c initializing of velocity field 
c ***

         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx-dx*0.5d0
           y = (j-1)*dy
           z = (k-1)*dz
        ux(i,j,k)=2*(sin(pi*x)**2)*sin(2*pi*y)*sin(2*pi*z)*cos(pi*t/TT) 
         enddo
         enddo
         enddo
         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx
           y = (j-1)*dy-dy*0.5d0
           z = (k-1)*dz 
        uy(i,j,k)= -sin(2*pi*x)*(sin(pi*y)**2)*sin(2*pi*z)*cos(pi*t/TT) 
         enddo
         enddo
         enddo
         do k=1,NZ
         do j=1,NY
         do i=1,NX
           x = (i-1)*dx
           y = (j-1)*dy
           z = (k-1)*dz-dz*0.5d0 
        uz(i,j,k)= -sin(2*pi*x)*sin(2*pi*y)*(sin(pi*z)**2)*cos(pi*t/TT)
         enddo
         enddo
         enddo
          
          tswap = tswap + 1
          if (tswap .gt. 3) tswap = 1
          
          if (MOD(tswap,3) .eq. 0) then
            call swp(uz,f,3)
            call swp(ux,f,1)
            call swp(uy,f,2)
          elseif (MOD(tswap,2) .eq. 0) then
            call swp(uy,f,2)
            call swp(uz,f,3)
            call swp(ux,f,1)
          else 
            call swp(ux,f,1)
            call swp(uy,f,2)
            call swp(uz,f,3)
          endif  
          
          call bc_c(f,nx,ny,nz)


         tmp = 0.0d0
         do k=1,NZ
         do j=1,NY
         do i=1,NX
            tmp = tmp +f(i,j,k)
         end do
         end do
         end do

         write(6,1001) step,t , 'Mass=',(tmp-init_f)*dx*dy*dz
 1001    format(i5,f,2x,a,d30.23)
         
         write(11,*)  step,t/TT, 'Mass=',(tmp-init_f)*dx*dy*dz
         

         if(mod(step,400*1).eq.0 .or. t/TT == 1.0 .or. t/TT == 2.0)then 
              call PARAVIEW(tec,NX,NY,NZ,f)
              !pause
         end if
         
         if(t/TT >= 1.0) then
         call PARAVIEW(tec,NX,NY,NZ,f)
         stop
         endif
         
         t = t+dt
         
      end do
      close(11)
      

      tmp=0.0d0
      tmp1=0.0d0
      do k=1,NZ
      do j=1,NY
      do i=1,NX
         tmp1 = tmp1 + f(i,j,k)
      end do
      end do
      end do
      tmp1 = (tmp1-init_f)/init_f
      write(6,*) '# Mass error=',tmp1
      
      deallocate(ux,uy,uz,f,fn) 

 1002 format(i5,2x,i5,2x,i5,2x,d30.23)

      stop
      end



c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c Split advection of the interface along the x (d=1), y (d=2) and z (d=3)
c direction
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      SUBROUTINE swp(us,c,d)
c***
      !include 'param.h'
      use param
      INTEGER i,j,k,invx,invy,invz,d
      DOUBLE PRECISION us(nx,ny,nz),c(nx,ny,nz),mx,my,mz,mm1,mm2
      DOUBLE PRECISION a1,a2,alpha,AL3D,FL3D
      DOUBLE PRECISION vof1(nx,ny,nz),vof2(nx,ny,nz),vof3(nx,ny,nz)
      INTRINSIC DMAX1,DMIN1
      EXTERNAL AL3D,FL3D
c***
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              a1 = us(i,j,k) *dt/delta
              if (d.eq.1) then
                a2 = us(i+1,j,k) *dt/delta
              elseif (d.eq.2) then
                a2 = us(i,j+1,k) *dt/delta
              elseif (d.eq.3) then
                a2 = us(i,j,k+1) *dt/delta
              endif
               
c***
c     3 cases: 1: DEFAULT (c=0. and fluxes=0.); 2: c=1.; 3:c>0.
c***
               vof1(i,j,k) = 0.0d0
               vof2(i,j,k) = 0.0d0
               vof3(i,j,k) = 0.0d0

               if (c(i,j,k) .EQ. 1.0d0) then
                  vof1(i,j,k) = DMAX1(-a1,0.d0)
                  vof2(i,j,k) = 1.d0 - DMAX1(a1,0.d0) + DMIN1(a2,0.d0)
                  vof3(i,j,k) = DMAX1(a2,0.d0)

               else if (c(i,j,k) .GT. 0.d0) then
c***
c     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
c     (3) get alpha;               (4) back to original plane;
c     (5) lagrangian advection;    (6) get fluxes
c*(1)*
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i-1,j+1,k-1)
     %                 +c(i-1,j+1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i-1,j+1,k)
     %                 +c(i-1,j,k-1)+c(i-1,j,k+1))+4.0d0*c(i-1,j,k)
                  mm2 = c(i+1,j-1,k-1)+c(i+1,j-1,k+1)+c(i+1,j+1,k-1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i+1,j-1,k)+c(i+1,j+1,k)
     %                 +c(i+1,j,k-1)+c(i+1,j,k+1))+4.0d0*c(i+1,j,k)
                  mx = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j-1,k+1)+c(i+1,j-1,k-1)
     %                 +c(i+1,j-1,k+1)+2.0d0*(c(i-1,j-1,k)+c(i+1,j-1,k)
     %                 +c(i,j-1,k-1)+c(i,j-1,k+1))+4.0d0*c(i,j-1,k)
                  mm2 = c(i-1,j+1,k-1)+c(i-1,j+1,k+1)+c(i+1,j+1,k-1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j+1,k)+c(i+1,j+1,k)
     %                 +c(i,j+1,k-1)+c(i,j+1,k+1))+4.0d0*c(i,j+1,k)
                  my = mm1 - mm2
                  
                  mm1 = c(i-1,j-1,k-1)+c(i-1,j+1,k-1)+c(i+1,j-1,k-1)
     %                 +c(i+1,j+1,k-1)+2.0d0*(c(i-1,j,k-1)+c(i+1,j,k-1)
     %                 +c(i,j-1,k-1)+c(i,j+1,k-1))+4.0d0*c(i,j,k-1)
                  mm2 = c(i-1,j-1,k+1)+c(i-1,j+1,k+1)+c(i+1,j-1,k+1)
     %                 +c(i+1,j+1,k+1)+2.0d0*(c(i-1,j,k+1)+c(i+1,j,k+1)
     %                 +c(i,j-1,k+1)+c(i,j+1,k+1))+4.0d0*c(i,j,k+1)
                  mz = mm1 - mm2
c*(2)*  
                  invx = 1
                  invy = 1
                  invz = 1
                  if (mx .LT. 0.0d0) then
                     mx = -mx
                     invx = -1
                  endif
                  if (my .LT. 0.0d0) then
                     my = -my
                     invy = -1
                  endif
                  if (mz .LT. 0.0d0) then
                     mz = -mz
                     invz = -1
                  endif
                  mm2 = mx+my+mz
                  mx = mx/mm2
                  my = my/mm2
                  mz = mz/mm2
c*(3)*  
                  alpha = AL3D(mx,my,mz,c(i,j,k))
c*(4)*  
                  mx = invx*mx
                  my = invy*my
                  mz = invz*mz
                  alpha = alpha + DMIN1(0.d0,mx) + DMIN1(0.d0,my) +
     %                 DMIN1(0.d0,mz)
c*(5)*  

                  mm1 = DMAX1(a1,0.0d0)
                  mm2 = 1.d0 - mm1 + DMIN1(0.d0,a2)
                  
                  if (d.eq.1) then
                    mx = mx/(1.0d0 - a1 + a2)
                    alpha = alpha + mx*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mx,my,mz,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mx,my,mz,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(mx,my,mz,alpha,mm1,mm2)
                  elseif (d.eq.2) then
                    my = my/(1.0d0 - a1 + a2)
                    alpha = alpha + my*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(my,mz,mx,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(my,mz,mx,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(my,mz,mx,alpha,mm1,mm2)
                  elseif (d.eq.3) then
                    mz = mz/(1.0d0 - a1 + a2)
                    alpha = alpha + mz*a1
                    if (a1 .LT. 0.d0) 
     %                 vof1(i,j,k) = FL3D(mz,mx,my,alpha,a1  ,-a1)
                    if (a2 .GT. 0.d0) 
     %                 vof3(i,j,k) = FL3D(mz,mx,my,alpha,1.d0,a2)
                       vof2(i,j,k) = FL3D(mz,mx,my,alpha,mm1,mm2)
                  endif
               endif
            enddo
         enddo
      enddo
c***  
c     (1) apply proper boundary conditions to fluxes
c     (2) new values of c and  clip it: 0. <= c <= 1.
c     (3) apply proper boundary conditions to c
c*(1)*  
      call bc_flux(vof1,vof3,us,nx,ny,nz,d)
c*(2)*  
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
              if (d.eq.1) then
                c(i,j,k) = vof3(i-1,j,k) + vof2(i,j,k) + vof1(i+1,j,k)
              elseif (d.eq.2) then
                c(i,j,k) = vof3(i,j-1,k) + vof2(i,j,k) + vof1(i,j+1,k)
              elseif (d.eq.3) then
                c(i,j,k) = vof3(i,j,k-1) + vof2(i,j,k) + vof1(i,j,k+1)
              endif
              c(i,j,k) = DMAX1(0.0d0,DMIN1(1.0d0,c(i,j,k)))
           enddo
         enddo
      enddo
c*(3)*
      call bc_c(c,nx,ny,nz)
c***
      return
      end





c***********************************************************************
      SUBROUTINE bc_c(scal,nx,ny,nz)
c***
 
      INTEGER nx,ny,nz, i, j, k
      DOUBLE PRECISION scal(nx,ny,nz),y0,z0,zz,yy,r,h
      DOUBLE PRECISION GP_DFETCH
      EXTERNAL GP_DFETCH
      h=1.d0/(nx-2)
c***  
c     periodicity in the x direction
c***  
!      IF(NEIGHB(1).LT.0) THEN
c      write(*,*)'BCC2 PPRANK',PPRANK,NEIGHB(1),NEIGHB(2)
      do k=2,nz-1
         do j=2,ny-1
            scal(1 ,j,k) = 0 !scal(2,j,k)
         enddo 
      enddo
!      ENDIF
c***  
!      IF(NEIGHB(2).LT.0) THEN
      do k=2,nz-1
         do j=2,ny-1
            scal(nx,j,k) = 0 !scal(nx-1,j,k)
         enddo 
      enddo
!      ENDIF
c***  
c     gradient equal to zero in the y and z directions
c***  
      do k=2,nz-1
         do i=1,nx
            scal(i,1,k ) = 0 !scal(i,2   ,k)
            scal(i,ny,k) = 0 !scal(i,ny-1,k)
         enddo
      enddo 
      do j=1,ny
         do i=1,nx
            scal(i,j,1 ) = 0 !scal(i,j,2)
            scal(i,j,nz) = 0 !scal(i,j,nz-1)
         enddo
      enddo 
c***
!      CALL BUPDAT3D(scal,NX,NY,NZ)
      return
      end




c***********************************************************************
      SUBROUTINE bc_flux(vof1,vof3,vel,nx,ny,nz,indx)
c***
 
      INTEGER nx, ny, nz, indx, i, j, k
      DOUBLE PRECISION vof1(nx,ny,nz), vof3(nx,ny,nz), vel(nx,ny,nz)
c***
c     along x direction
c***
      if (indx .eq. 1) then
!        IF(NEIGHB(1).LT.0) THEN
          do k=2,nz-1
            do j=2,ny-1
               vof3(1 ,j,k) = 0.d0
            enddo
          enddo
!        ENDIF
!        IF(NEIGHB(2).LT.0) THEN
          do k=2,nz-1
            do j=2,ny-1
               vof1(nx,j,k) = 0.d0
            enddo
          enddo
!        ENDIF
c***
c     along y direction
c***
      elseif (indx .eq. 2) then
         do k=2,nz-1
            do i=2,nx-1
               vof3(i,1 ,k) = 0 !DMAX1(0.0d0,vel(i,2,k))
               vof1(i,ny,k) = 0.d0
            enddo
         enddo
c***
c     along z direction
c***
      elseif (indx .eq. 3) then
         do j=2,ny-1
            do i=2,nx-1
               vof3(i,j,1) = 0 !DMAX1(0.0d0,vel(i,j,2))
               vof1(i,j,nz) = 0.0d0
            enddo
         enddo
c***  
c     wrong value
c***  
      else 
         stop 'wrong value for indx in one of the swaps'
      endif 
c***
!      CALL BUPDAT3D(vof1,nx,ny,nz)
!      CALL BUPDAT3D(vof3,nx,ny,nz)
      return
      end






c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND THE "CUT VOLUME" V0 GIVEN r0, dr0 AND
c m1 x1 + m2 x2 + m3 x3 = alpha
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION FL3D(m1,m2,m3,alpha,r0,dr0)
c***
      DOUBLE PRECISION m1,m2,m3,alpha,r0,dr0, vm1,vm2,vm3,vm12,a,v
      DOUBLE PRECISION al,al0,n1,n2,n3,b1,b2,b3,b12,bm,tmp,pr,CONST_TINY
      INTRINSIC DMAX1,DMIN1,DABS
      DOUBLE PRECISION, parameter :: ONE = 1.0d0, PB = 1.49d0, 
     &        PC2 = 0.239d0, PC1 = 0.132d0, 
     &        PC0 = (PB * (PB * PC2 + 4d0 * PC1 - 8d0) / 16d0),  
     &        PA = (PB * PB * (PB - 1d0))
      
      CONST_TINY = 1D-25
c***
c     (1) move origin to r0 along r ;  (2) reflect parallelepiped;
c     (3) limit alpha (0<= al0 <=0.5); (4) order coefficients: b1<b2<b3;
c     (5) calculate volume (NOTE: it is assumed:s0=t0=0; ds0=dt0=1.)
c*(1)*
      al = alpha - m1*r0
c*(2)*
      al = al + DMAX1(0.d0,-m1*dr0)+DMAX1(0.d0,-m2)+DMAX1(0.d0,-m3)
      tmp = DABS(m1)*dr0 + DABS(m2) + DABS(m3)
      n1 = DABS(m1)/tmp
      n2 = DABS(m2)/tmp
      n3 = DABS(m3)/tmp
      al = DMAX1(0.d0,DMIN1(1.d0,al/tmp))
c*(3)*
      al0 = DMIN1(al,1.d0-al)
c*(4)*
      b1 = DMIN1(n1*dr0,n2)
      b3 = DMAX1(n1*dr0,n2)
      b2 = n3
      if (b2 .LT. b1) then
         tmp = b1
         b1 = b2
         b2 = tmp
      else if (b2 .GT. b3) then
         tmp = b3
         b3 = b2
         b2 = tmp
      endif
      b12 = b1 + b2
      bm = DMIN1(b12,b3)
      pr = DMAX1(6.d0*b1*b2*b3,1.0d-50)
      
c*5*     

      
      
! Aoki PLIC method
        vm1 = b1 !min(vma, vmb, vmc) ! Sort vm1, vm2, vm3
        vm3 = b3 !max(vma, vmb, vmc) ! such that vm1 �� vm2 �� vm3 .
        vm2 = b2 !abs(1d0 - vm3 - vm1)
        vm12 = b12 !vm1 + vm2
        a = al0
        v = 0d0
        if (a > 0.0d0) then
            if (a < vm1) then
                v = a ** 3 / (6d0 * vm1 * vm2 * vm3)
            else if (a < vm2) then
                v = a * (a - vm1) / (2d0 * vm2 * vm3) +  
     &              vm1 ** 2 / (6d0 * vm2 * vm3 + CONST_TINY)
            else if (a < min(vm12, vm3)) then
                v = (a ** 2 * (3d0 * vm12 - a) + 
     &              vm1 ** 2 * (vm1 - 3d0 * a) + 
     &              vm2 ** 2 * (vm2 - 3d0 * a)) / 
     &              (6d0 * vm1 * vm2 * vm3)
            else if (vm3 < vm12) then
                v = (a ** 2 * (3d0 - 2d0 * a) + 
     &              vm1 ** 2 * (vm1 - 3d0 * a) +  
     &              vm2 ** 2 * (vm2 - 3d0 * a) + 
     &              vm3 ** 2 * (vm3 - 3d0 * a)) / 
     &              (6d0 * vm1 * vm2 * vm3)
            else
                v = (a - 0.5d0 * vm12) / vm3
            end if
        endif
        tmp = v
      
        
        
        
        
      
      if (al .LE. 0.5d0) then
         FL3D = tmp*dr0
      else
         FL3D = (1.d0-tmp)*dr0
      endif
c***  
      return
      end



c *** *** *** *** *** *** *** *** *** *** *** *** 
c    convert level set phi to VOF function 
c *** *** *** *** *** *** *** *** *** *** *** ***
      subroutine levelset2vof(nx,ny,nz,ls,cc)

      double precision cc(nx,ny,nz),ls(nx,ny,nz)
      double precision zero, one, normL1
      double precision mm1,mm2,mm3,mx,my,mz,alpha
      double precision fl3d
      integer i,j,k

      zero=0.d0
      one=1.d0
      do k=2,nz-1
         do j=2,ny-1
            do i=2,nx-1
c***
c     (1) normal vector: mx,my,mz; (2) mx,my,mz>0. and mx+my+mz = 1.;
c     (3) shift alpha to origin    (4) get volume from alpha;  
c***

               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j-1,k+1)+ls(i-1,j+1,k-1)
     %              +ls(i-1,j+1,k+1)+2.0d0*(ls(i-1,j-1,k)+ls(i-1,j+1,k)
     %              +ls(i-1,j,k-1)+ls(i-1,j,k+1))+4.0d0*ls(i-1,j,k)
               mm2 = ls(i+1,j-1,k-1)+ls(i+1,j-1,k+1)+ls(i+1,j+1,k-1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i+1,j-1,k)+ls(i+1,j+1,k)
     %              +ls(i+1,j,k-1)+ls(i+1,j,k+1))+4.0d0*ls(i+1,j,k)
               mx = (mm1 - mm2)/32.d0
               
               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j-1,k+1)+ls(i+1,j-1,k-1)
     %              +ls(i+1,j-1,k+1)+2.0d0*(ls(i-1,j-1,k)+ls(i+1,j-1,k)
     %              +ls(i,j-1,k-1)+ls(i,j-1,k+1))+4.0d0*ls(i,j-1,k)
               mm2 = ls(i-1,j+1,k-1)+ls(i-1,j+1,k+1)+ls(i+1,j+1,k-1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i-1,j+1,k)+ls(i+1,j+1,k)
     %              +ls(i,j+1,k-1)+ls(i,j+1,k+1))+4.0d0*ls(i,j+1,k)
               my = (mm1 - mm2)/32.d0
               
               mm1 = ls(i-1,j-1,k-1)+ls(i-1,j+1,k-1)+ls(i+1,j-1,k-1)
     %              +ls(i+1,j+1,k-1)+2.0d0*(ls(i-1,j,k-1)+ls(i+1,j,k-1)
     %              +ls(i,j-1,k-1)+ls(i,j+1,k-1))+4.0d0*ls(i,j,k-1)
               mm2 = ls(i-1,j-1,k+1)+ls(i-1,j+1,k+1)+ls(i+1,j-1,k+1)
     %              +ls(i+1,j+1,k+1)+2.0d0*(ls(i-1,j,k+1)+ls(i+1,j,k+1)
     %              +ls(i,j-1,k+1)+ls(i,j+1,k+1))+4.0d0*ls(i,j,k+1)
               mz = (mm1 - mm2)/32.D0
c     *(2)*  
               mx = DABS(MX)
               MY = DABS(MY)
               MZ = DABS(MZ)
               normL1 = mx+my+mz
               mx = mx/normL1
               my = my/normL1
               mz = mz/normL1

               alpha = ls(i,j,k)/normL1
c     *(3)*  
               alpha = alpha + 0.5d0
c     *(4)*  
               if(alpha.ge.1.d0) then 
                  cc(i,j,k) = 1.d0
               else if (alpha.le.0.d0) then
                  cc(i,j,k) = 0.d0
               else 
                  cc(i,j,k) = FL3D(mx,my,mz,alpha,zero,one)
c                  write(*,*) cc(i,j,k)
               end if
            enddo
         enddo
      enddo 
      return
      end












!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      SUBROUTINE PARAVIEW(tec,NX,NY,NZ,P)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!----------------------------------------------------------------------- 
!----------------------------------------------------------------------- 

      character*30 tec, ch
      real*8  ::  P(NX,NY,NZ)
      
      real*8,allocatable ::   tp(:)
      

      allocate(tp(nx*ny*nz))
      

      n = 0
      do k = 1,nz
      do j = 1,ny
      do i = 1,nx
        n = n + 1
        tp(n) = P(i,j,k)
      enddo
      enddo
      enddo
      
      open (unit=1, file=tec)
      write(1,'(a)') '# vtk DataFile Version 3.0'
      write(1,'(a)') 'Non-uniform Rectilinear - Rectilinear Grid' 
      write(1,'(a)') 'ASCII' 
      write(1,'(a)') 'DATASET RECTILINEAR_GRID' 
      write(1,'(a,3i)') 'DIMENSIONS',nx,ny,nz
      write(1,'(a,i,a)') 'X_COORDINATES',nx,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,nx)
      write(1,'(a,i,a)') 'Y_COORDINATES',ny,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,ny)
      write(1,'(a,i,a)') 'Z_COORDINATES',nz,'  float'
      write(1,'(f10.3)') (i*1.0,i=1,nz)

      write(1,'(a,i)') 'POINT_DATA',nx*ny*nz
      write(1,'(a)') 'SCALARS P float 1'
      write(1,'(a)') 'LOOKUP_TABLE default'
      write(1,'(F)')((tp(i)),i=1,n)
      close(1)
      

      deallocate (tp)
      ENDSUBROUTINE
!-----------------------------------------------------------------------



      subroutine init(f,ux,uy,uz,t,dx,dy,dz)
      !include 'param.h'
      use param
      real*8 f(NX,NY,NZ)
      real*8 ux(NX,NY,NZ),uy(NX,NY,NZ),uz(NX,NY,NZ)
      real*8 t,dx,dy,dz

      real*8  ls(NX,NY,NZ)
      integer i,j,k,is,ii
      real*8  x, y, z, x0, y0, z0, R, tmp 

      t  = 0.0d0
      dx = 1.0/NX !0.01d0
      dy = 1.0/NY !0.01d0
      dz = 1.0/NZ !0.01d0 
      
      R  = 0.15
      x0 = 0.35
      y0 = 0.35
      z0 = 0.35

      do k=1,NZ
      do j=1,NY
      do i=1,NX

         f(i,j,k)  = 0.0d0

         x = (i-1)*dx
         y = (j-1)*dy
         z = (k-1)*dz

c ***
c initializing of level set function 
c ***
         ls(i,j,k) = R - sqrt( (x-x0)**2 + (y-y0)**2 + (z-z0)**2 ) 
      end do
      end do
      end do

      call levelset2vof(nx,ny,nz,ls,f)

      
c ***

      return
      end


c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
c PROGRAM TO FIND alpha IN: m1 x1 + m2 x2 + m3 x3 = alpha,
c GIVEN m1+m2+m3=1 (all > 0) AND THE VOLUMETRIC FRACTION cc
c ****** 1 ******* 2 ******* 3 ******* 4 ******* 5 ******* 6 ******* 7 *
      DOUBLE PRECISION FUNCTION AL3D(b1,b2,b3,cc)
c***
C      INCLUDE   'param.h'
      use param
      DOUBLE PRECISION m1,m2,m3,cc,b1,b2,b3,tmp,pr,ch,mm,m12
      DOUBLE PRECISION p,p12,q,teta,cs
      DOUBLE PRECISION UNTIER,V1,V2,V3
      DOUBLE PRECISION alpha, w, vm1, vm2, vm3, vm12 
      DOUBLE PRECISION a0, a1, a2, q0, sp, th ,CONST_TINY,CONST_PI
      
      DOUBLE PRECISION xi, invp, vma, vmb, vmc
      DOUBLE PRECISION, parameter :: ONE = 1.0d0, PB = 1.49d0, 
     &   PC2 = 0.239d0, PC1 = 0.132d0, 
     &   PC0 = (PB * (PB * PC2 + 4d0 * PC1 - 8d0) / 16d0), 
     &   PA = (PB * PB * (PB - 1d0))
      
      PARAMETER (UNTIER=1.d0/3.d0)
      INTRINSIC DMAX1,DMIN1,DSQRT,DACOS,DCOS
      
      CONST_TINY = 1D-25
      CONST_PI = 3.14159265358979323846d0
c***  
c     (1) order coefficients: m1<m2<m3; (2) get ranges: V1<V2<v3;
c     (3) limit ch (0.d0 < ch < 0.5d0); (4) calculate alpha
c*(1)* 
      m1 = DMIN1(b1,b2)
      m3 = DMAX1(b1,b2)
      m2 = b3
      if (m2 .LT. m1) then
         tmp = m1
         m1 = m2
         m2 = tmp
      else if (m2 .GT. m3) then
         tmp = m3
         m3 = m2
         m2 = tmp
      endif
c*(2)*
      m12 = m1 + m2 
      pr  = DMAX1(6.d0*m1*m2*m3,1.d-50)
      V1  = m1*m1*m1/pr
      V2  = V1 + 0.5d0*(m2-m1)/m3
      if (m3 .LT. m12) then
         mm = m3
         V3 = (m3*m3*(3.d0*m12-m3) + m1*m1*(m1-3.d0*m3) +
     %        m2*m2*(m2-3.d0*m3))/pr
      else
         mm = m12
         V3 = 0.5d0*mm/m3
      endif
c*(3)*
      ch = DMIN1(cc,1.d0-cc)
c*(4)*      
      
!original method      
      if (ch .LT. V1) then
c***         AL3D = cbrt(pr*ch)
         AL3D = (pr*ch)**UNTIER
      else if (ch .LT. V2) then
         AL3D = 0.5d0*(m1 + DSQRT(m1*m1 + 8.d0*m2*m3*(ch-V1)))
      else if (ch .LT. V3) then
         p = 2.d0*m1*m2
         q = 1.5d0*m1*m2*(m12 - 2.d0*m3*ch)
         p12 = DSQRT(p)
         teta = DACOS(q/(p*p12))/3.d0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + m12
      else if (m12 .LT. m3) then
         AL3D = m3*ch + 0.5d0*mm
      else 
         p = m1*(m2+m3) + m2*m3 - 0.25d0
         q = 1.5d0*m1*m2*m3*(0.5d0-ch)
         p12 = DSQRT(p)
         teta = DACOS(q/(p*p12))/3.0
         cs = DCOS(teta)
         AL3D = p12*(DSQRT(3.d0*(1.d0-cs*cs)) - cs) + 0.5d0
      endif

      if (cc .GT. 0.5d0)  AL3D = 1.d0 - AL3D
      

      
c***
      return
      end