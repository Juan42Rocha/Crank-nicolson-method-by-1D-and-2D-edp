


set term wxt persist size 999,333


set view map
#set dgrid3d 60,60,2
#set pm3d at b  map
nt=100


do for [k=3:nt+3]{
set multiplot layout 1,3
  set size square
  spl "dat2.out" i 0 u  1:2:(column(k)) w imag at b # w pm3d  t "E. explicito"

  spl "dat2.out" i 1 u  1:2:(column(k)) w imag at b#w pm3d  t "E. implicito"
  set title sprintf("%d",k-2)
  set key o t c

  spl "dat2.out" i 2 u  1:2:(column(k)) w imag at b#w pm3d  t "CM"
  
  pause(0.1)
}

