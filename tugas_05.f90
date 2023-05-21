program interpolasi_lagrange_prog
implicit none
integer, parameter :: dpr=kind(1.0D0)
integer :: i, j, k, n, size, xmaks, xmin, metode
real(dpr) :: b, m
real(dpr), allocatable :: D(:,:), X(:), Y(:), x0(:), p(:)
character(len=20), dimension(6) :: list_nama_file 
character(len=20) :: output_x, output_y

list_nama_file(1) = 'data_plot_xl.txt'
list_nama_file(2) = 'data_plot_yl.txt'
list_nama_file(3) = 'data_plot_xlk.txt'
list_nama_file(4) = 'data_plot_ylk.txt'
list_nama_file(5) = 'data_plot_xhk.txt'
list_nama_file(6) = 'data_plot_yhk.txt'

! membaca jumlah data
open(unit=1,file='input_matrix.txt',action='read')
read(1,*) size
close(1)

write(*,*) "Jumlah data input = ", size
write(*,*) ' '
! jumlah titik data yang ingin dihasilkan
write(*,*) 'masukkan jumlah titik data yang ingin dihasilkan'
write(*,*) 'jumlah titik data output ='
read(*,*) n

allocate(D(1:size,1:2),X(1:size),Y(1:size),x0(1:n),p(1:n))
write(*,*) ' '

! data untuk plot grafik
write(*,*) 'masukkan rentang nilai x yang diiginkan'
write(*,*) 'xmaks ='
read(*,*) xmaks
write(*,*) 'xmin ='
read(*,*) xmin
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

! memilih metode interpolasi
write(*,*) 'Pilih metode interpolasi yang ingin digunakan'
write(*,*) 'ketik angka numerik antara 1-3'
write(*,*) '1 = metode interpolasi lagrange'
write(*,*) '2 = metode interpolasi lagrange kubik'
write(*,*) '3 = metode interpolasi hermite kubik'
read(*,*) metode

! menghitung pasangan titik dari x0 yaitu p
if (metode==1) then 
   p = interpolasi_lagrange(X,Y,x0,n,size)
   output_x = list_nama_file(1)
   output_y = list_nama_file(2)
else if (metode==2) then
   p = interpolasi_lagrange_kubik(X,Y,x0,n,size)
   output_x = list_nama_file(3)
   output_y = list_nama_file(4)
else if (metode==3) then
   p = interpolasi_hermite_kubik(X,Y,x0,n,size)
   output_x = list_nama_file(5)
   output_y = list_nama_file(6)
else 
   write(*,*) 'Pilihan metode tidak valid, pilih metode antara 1-3'
end if

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

open(unit=4,file=output_x, status='replace',action='write')
do i=1,n
   write(4,*) x0(i)
end do
close(4)

open(unit=5,file=output_y, status='replace',action='write')
do i=1,n
   write(5,*) p(i)
end do
close(5)

write(*,*) ''
write(*,*) 'untuk melihat plot grafik, jalankan program matlab plot_tugas05.m di MatLab atau Octave'

deallocate(D,X,Y,x0,p)

stop

contains

function interpolasi_lagrange(X,Y,x0,n,size) result(yl)
integer :: i, j, k
integer, intent(in) :: n, size
real(dpr), intent(in) :: X(size), Y(size), x0(n)
real(dpr) :: b, m
real(dpr) :: u(n,size,size), L(n,size), yl(n)

! melakukan interpolasi lagrange
do k = 1,n 
   do i = 1,size
      do j = 1,size
         if (i==j) then
         u(k,i,j) = 1
         else
         u(k,i,j) = (x0(k) - X(j))/(X(i)-X(j))
         endif
      end do
   end do
end do

do k = 1,n 
   do i = 1,size 
   L(k,i) = 1.0
   end do
end do

do k = 1,n
   do i = 1,size
      do j = 1,size
         L(k,i) = L(k,i) * u(k,i,j)
      end do
   end do
end do

! menghitung pasangan titik dari x0 yaitu yl
do k = 1,n 
   yl(k) = 0.0
end do

do k = 1,n 
   do i = 1,size 
      yl(k) = yl(k) + L(k,i) * Y(i)
   end do
end do

end function interpolasi_lagrange

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

function interpolasi_hermite_kubik(XX,YY,x0,n,size) result(yhk)
integer :: i, j, k, ll, in_min
integer, intent(in) :: n, size
real(dpr), intent(in) :: XX(size), YY(size), x0(n)
real(dpr) :: b, m
real(dpr) :: h1(n,2), h2(n,2), h(n,4), yhk(n), X(size+2), Y(size+2)

! melakukan interpolasi lagrange kubik
! inisiasi nilai
do i = 1,size+2
   X(i) = 0.0
   Y(i) = 0.0
end do 

do i = 1,size
   X(i+1) = XX(i)
   Y(i+1) = YY(i)
end do

do k = 1,n
   ! mengecek posisi x0(k)
   do ll = 2, size+1
      if (x0(k) >= X(ll) .and. x0(k) <= X(ll+1)) then
         in_min = ll-1
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

end program interpolasi_lagrange_prog