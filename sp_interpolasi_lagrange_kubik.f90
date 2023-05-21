program interpolasi_lagrange_kubik_prog
implicit none
integer, parameter :: dpr=kind(1.0D0)
integer :: i, j, k, ll, n, size, xmaks, xmin
real(dpr) :: b, m
real(dpr), allocatable :: D(:,:), X(:), Y(:), x0(:), p(:)

! membaca jumlah data
open(unit=1,file='input_matrix.txt',action='read')
read(1,*) size
close(1)

! jumlah titik data yang ingin dihasilkan
write(*,*) 'masukkan jumlah titik data yang ingin dihasilkan'
read(*,*) n

write(*,*) "Jumlah data = ", size

allocate(D(1:size,1:2),X(1:size),Y(1:size),x0(1:n),p(1:n))

! data untuk plot grafik
xmaks = 60
xmin = 0
m = n
b = (xmaks-xmin)/(m-1)
do i=1,n
   x0(i) = xmin + (i-1)*b 
end do

! membaca data    
open(unit=1,file='input_matrix.txt',action='read')
do i=1,size
    read(1,*,end=2) D(i,:)
    X(i)=D(i,1) !data X
    Y(i)=D(i,2) !data Y
2 end do

! menghitung pasangan titik dari x0 yaitu p

p = interpolasi_lagrange_kubik(X,Y,x0,n,size)

! menyimpan data ke file eksternal
write(*,*) 'Data tersimpan file eksternal'
open(unit=2,file='data_plot_x0.txt',status='replace',action='write')
do i=1,size
   write(2,*) X(i)
end do
close(2)

open(unit=3,file='data_plot_y0.txt',status='replace',action='write')
do i=1,size
   write(3,*) Y(i)
end do
close(3)

open(unit=4,file='data_plot_xlk.txt', status='replace',action='write')
do i=1,n
   write(4,*) x0(i)
end do
close(4)

open(unit=5,file='data_plot_ylk.txt', status='replace',action='write')
do i=1,n
   write(5,*) p(i)
end do
close(5)

deallocate(D,X,Y,x0,p)

stop

contains

function interpolasi_lagrange_kubik(X,Y,x0,n,size) result(ylk)
integer :: i, j, k, ll, in_min
integer, intent(in) :: n, size
real(dpr), intent(in) :: X(size), Y(size), x0(n)
real(dpr) :: b, m
real(dpr) :: u(n,size,size), L(n,size), ylk(n)

! melakukan interpolasi lagrange kubik
! inisiasi nilai

do k = 1,n 
   do i = 1,4
      L(k,i) = 1.0
      do j = 1,4
         u(k,i,j) = 0.0
      end do
   end do
end do

do k = 1,n
   ! mengecek posisi x0(k)
   do ll = 1, size-5
      if (x0(k) < X(3)) then
         in_min = 1
         exit
      else if (x0(k) > X(size-2)) then
         in_min = size-3
         exit
      else if (x0(k) >= X(ll+2) .and. x0(k) <= X(ll+3)) then
         in_min = ll+1
         exit
      end if
   end do
   ! didapat indeks minimum in_min

   ! menghitung L(k,i)
   do i = 1,4
      do j = 1,4
         if (i==j) then
            u(k,i,j) = 1
         else
            u(k,i,j) = (x0(k) - X(j+in_min-1))/(X(i+in_min-1)-X(j+in_min-1))
         end if
      end do
   end do

   do i = 1,4
      do j = 1,4
         L(k,i) = L(k,i) * u(k,i,j)
      end do
   end do

   ! menghitung pasangan titik dari x0 yaitu ylk
   ylk(k) = 0.0
   do i = in_min,in_min+3
      ylk(k) = ylk(k) + L(k,i+1-in_min) * Y(i)
   end do

end do

end function interpolasi_lagrange_kubik

end program interpolasi_lagrange_kubik_prog