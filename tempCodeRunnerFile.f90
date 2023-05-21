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