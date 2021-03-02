! #######################################################
! TO COMPUTE THE YANGLOM LAW AS IN HELLINGER PAPER
! #######################################################


subroutine yanglommix(ax,ay,az,bx,by,bz,jx,jy,jz,nx,ny,nz,lagx,lagy,lagz,nlagx,nlagy,nlagz,S,Su,Sb,Yxx,Yyy,Yzz,Hx,Hy,Hz)
   implicit none
   integer*8, intent(in) :: nlagx,nlagy,nlagz
   integer*8, intent(in) :: nx,ny,nz
   integer*8 :: xx,yy,zz,ix,iy,iz,x,y,z,kx,ky,kz
   double precision, intent(in), dimension(nx,ny,nz) :: ax,ay,az,bx,by,bz,jx,jy,jz
   double precision, intent(in), dimension(nlagx) :: lagx
   double precision, intent(in), dimension(nlagy) :: lagy
   double precision, intent(in), dimension(nlagz) :: lagz

   double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: S 
   double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Su 
   double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Sb
   double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Yxx,Yyy,Yzz,Hx,Hy,Hz
   !double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Yyy
   !double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Yzz
   !double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Hx
   !double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Hy
   !double precision, intent(out), dimension(nlagx,nlagy,nlagz) :: Hz
 
! LOCAL VARIABLES
   double precision :: dux, duy, duz, absu, dbx, dby, dbz, absb, djx, djy, djz, stu, stb, sstu,sstb,ss
   double precision :: dotub, dotjb
   double precision :: nn
   double precision :: St,Ytxx,Ytyy,Ytzz,Htxx,Htyy,Htzz


   nn=nx*ny*nz

   do kz=1,nlagz; do ky=1,nlagy; do kx=1,nlagx
      iz=lagz(kz)
      iy=lagy(ky)
      ix=lagx(kx)
          
      St=0.
      sstu=0.
      sstb=0.
      Ytxx=0.
      Ytyy=0.
      Ytzz=0.
      Htxx=0.
      Htyy=0.
      Htzz=0.
!      print *,kz,ky,kx
      do z=1,nz; do y=1,ny; do x=1,nx
         call roll_idx(x,nx,ix,xx)
         call roll_idx(y,ny,iy,yy)
         call roll_idx(z,nz,iz,zz)

         dux=ax(xx,yy,zz)-ax(x,y,z)
         duy=ay(xx,yy,zz)-ay(x,y,z)
         duz=az(xx,yy,zz)-az(x,y,z)
         absu=sqrt(dux**2+duy**2+duz**2)

         dbx=bx(xx,yy,zz)-bx(x,y,z)
         dby=by(xx,yy,zz)-by(x,y,z)
         dbz=bz(xx,yy,zz)-bz(x,y,z)
         absb=sqrt(dbx**2+dby**2+dbz**2)

         djx=jx(xx,yy,zz)-jx(x,y,z)
         djy=jy(xx,yy,zz)-jy(x,y,z)
         djz=jz(xx,yy,zz)-jz(x,y,z)
         
         stu=absu**2
            
         stb=absb**2
         
         sstu=sstu+stu
        
         sstb=sstb+stb    
         
         ss=stu+stb
     

         dotub= dux*dbx+duy*dby+duz*dbz      
         dotjb= dbx*djx+dby*djy+dbz*djz   
         
         St=St+ss

         Ytxx=Ytxx+(dux*ss-2*dbx*dotub)

         Ytyy=Ytyy+(duy*ss-2*dby*dotub)

         Ytzz=Ytzz+(duz*ss-2*dbz*dotub)
         
         Htxx=Htxx+(2*dbx*dotjb-djx*stb)

         Htyy=Htyy+(2*dby*dotjb-djy*stb)

         Htzz=Htzz+(2*dbz*dotjb-djz*stb)

      enddo; enddo; enddo
      S(kx,ky,kz)=St/nn

      Su(kx,ky,kz)=sstu/nn
      
      Sb(kx,ky,kz)=sstb/nn
     
      Yxx(kx,ky,kz)=Ytxx/nn

      Yyy(kx,ky,kz)=Ytyy/nn
     
      Yzz(kx,ky,kz)=Ytzz/nn

      Hx(kx,ky,kz)=Htxx/nn 

      Hy(kx,ky,kz)=Htyy/nn

      Hz(kx,ky,kz)=Htzz/nn


   enddo; enddo; enddo 
 
end subroutine yanglommix
