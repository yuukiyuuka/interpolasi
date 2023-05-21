program interpolasi_hermite_kubik_prog
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

p = interpolasi_hermite_kubik(X,Y,x0,n,size)

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

open(unit=4,file='data_plot_xhk.txt', status='replace',action='write')
do i=1,n
   write(4,*) x0(i)
end do
close(4)

open(unit=5,file='data_plot_yhk.txt', status='replace',action='write')
do i=1,n
   write(5,*) p(i)
end do
close(5)

deallocate(D,X,Y,x0,p)

stop

contains

function interpolasi_hermite_kubik(X,Y,x0,n,size) result(yhk)
integer :: i, j, k, ll, in_min
integer, intent(in) :: n, size
real(dpr), intent(in) :: X(size), Y(size), x0(n)
real(dpr) :: b, m
real(dpr) :: h1(n,2), h2(n,2), h(n,4), yhk(n)

! melakukan interpolasi lagrange kubik
! inisiasi nilai

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
   
   ! menghitung h1(k,i) dan h2(k,i)
   h1(k,1) = (1-2*(x0(k)-X(in_min+1))/(X(in_min+1)-X(in_min+2)))*((x0(k)-X(in_min+2))/(X(in_min+1)-X(in_min+2)))**2 
   h1(k,2) = (1-2*(x0(k)-X(in_min+2))/(X(in_min+2)-X(in_min+1)))*((x0(k)-X(in_min+1))/(X(in_min+2)-X(in_min+1)))**2
   h2(k,1) = (x0(k)-X(in_min+1))*((x0(k)-X(in_min+2))/(X(in_min+1)-X(in_min+2)))**2
   h2(k,2) = (x0(k)-X(in_min+2))*((x0(k)-X(in_min+1))/(X(in_min+2)-X(in_min+1)))**2
   
   ! menghitung h(k,i)
   h(k,1) = h2(k,1)*((X(in_min+1)-X(in_min+2))/(X(in_min+0)-X(in_min+2)))/(X(in_min+0)-X(in_min+1))
   h(k,2) = h1(k,1) + h2(k,1)*(1/(X(in_min+1)-X(in_min+2)) + 1/(X(in_min+1)-X(in_min+0))) &
          + h2(k,2)*((X(in_min+2)-X(in_min+3))/(X(in_min+1)-X(in_min+3)))/(X(in_min+1)-X(in_min+2))
   h(k,3) = h1(k,2) + h2(k,2)*(1/(X(in_min+2)-X(in_min+1)) + 1/(X(in_min+2)-X(in_min+3))) &
          + h2(k,1)*((X(in_min+1)-X(in_min+0))/(X(in_min+2)-X(in_min+0)))/(X(in_min+2)-X(in_min+1))
   h(k,4) = h2(k,2)*((X(in_min+2)-X(in_min+1))/(X(in_min+3)-X(in_min+1)))/(X(in_min+3)-X(in_min+2))

   ! menghitung pasangan titik dari x0 yaitu yhk
   yhk(k) = 0.0
   do i = in_min,in_min+3
      yhk(k) = yhk(k) + h(k,i+1-in_min) * Y(i)
   end do

end do

end function interpolasi_hermite_kubik

end program interpolasi_hermite_kubik_prog