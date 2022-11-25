


set term wxt persist 


set view map
#set dgrid3d 60,60,2
#set pm3d at b  map
nt=1000


do for [k=3:nt+3:50]{
  set size square
  spl "dat5.out" i 0 u  1:2:(column(k)) w imag at b  t "CN"

  set title sprintf("%d",k-2)

  
  pause(0.001)
}

