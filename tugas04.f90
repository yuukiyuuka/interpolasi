program tugas04

implicit none
integer, parameter :: dpr=kind(1.0D0)
integer :: i,n,m,k,j,row
real(dpr) :: S
real(dpr), allocatable :: D(:,:),X(:),Y(:),Z(:),C(:,:),B(:),A(:)

!jumlah variabel yang tidak diketahui
write(*,*) 'masukkan ordo polinomial dan jumlah data yang ingin digunakan'
read(*,*) m,row

n = m+1

allocate(D(1:row,1:2),X(1:row),Y(1:row),Z(1:row),C(1:n,1:n),B(1:n),A(1:n))

!baca data matriks
open(unit=1,file='input_matrix.txt',action='read')
do i=1,row
    read(1,*,end=2) D(i,:)
    X(i)=D(i,1) !data X
    Y(i)=D(i,2) !data Y
end do

2 close(1)

!Inisialisasi matriks C dan B serta nilai S dengan nilai 0
S = 0.0
do i=1,row
   Z(i) = 0.0
end do

do i=1,n
   C(i,:) = 0.0
   B(i) = 0.0
end do

!menghitung matriks C
do k=1,n
   do j=1,n 
      do i=1,row
         C(k,j)=C(k,j)+X(i)**(k+j-2)
      end do
   end do
end do

!menghitung matriks B
do k=1,n
   do i=1,row
      B(k)=B(k)+Y(i)*X(i)**(k-1)
   end do
end do

!mencari matriks koefisien A
A=lu_solve(C,B,n)

!menghitung nilai S
do i=1,row
    do j=1,n 
       Z(i) = Z(i) + A(j)*X(i)**(j-1)
    end do 
end do

do i=1,row
   S = S + (Y(i) - Z(i))**2
end do

!save matriks C,B, dan A ke file eksternal
write(*,*) 'cek file eksternal'
open(unit=3,file='matriks_output.txt',status='replace',action='write')
write(3,*) 'Matriks koefisien A'
do i=1,n
   write(3,*) "a_",i-1,A(i)
end do
write(3,*) 'Nilai S'
write(3,*) S
close(3)

open(unit=4,file='koefisien_pol_orde_7.txt',status='replace',action='write')
do i=1,n 
   write(4,*) A(i)
end do
close(4)

deallocate(D,X,Y,Z,C,B,A)

stop

contains

function lu_solve(AA, BB, n) result(X)
  implicit none
  integer, intent(in) :: n
  real(dpr), intent(in) :: AA(n,n)
  real(dpr), intent(in) :: BB(n)
  real(dpr) :: X(n)
  integer :: i, j, k, maxind
  real(dpr) :: maxval, temp

  !copy AA and BB to temporary arrays
  real(dpr) :: A(n,n), B(n)
  A = AA
  B = BB

  ! Perform LU decomposition with partial pivoting
  do k = 1, n-1
    ! Find pivot element
    maxind = k
    maxval = abs(A(k,k))
    do i = k+1, n
      if (abs(A(i,k)) > maxval) then
        maxind = i
        maxval = abs(A(i,k))
      end if
    end do
    ! Swap rows if necessary
    if (maxind /= k) then
      do j = 1, n
        temp = A(k,j)
        A(k,j) = A(maxind,j)
        A(maxind,j) = temp
      end do
      temp = B(k)
      B(k) = B(maxind)
      B(maxind) = temp
    end if
    ! Perform elimination
    do i = k+1, n
      temp = A(i,k) / A(k,k)
      A(i,k) = temp
      do j = k+1, n
        A(i,j) = A(i,j) - temp * A(k,j)
      end do
      B(i) = B(i) - temp * B(k)
    end do
  end do

  ! Solve the system by back substitution
  X(n) = B(n) / A(n,n)
  do i = n-1, 1, -1
    X(i) = B(i)
    do j = i+1, n
      X(i) = X(i) - A(i,j) * X(j)
    end do
    X(i) = X(i) / A(i,i)
  end do

end function lu_solve

end program tugas04