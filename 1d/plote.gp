#set term wxt persist

set term gif animate delay 5
set output "rig.gif"


u(x,t,a,L)=exp(-0.03*((x+a*t)-L/2)**2)


do for [k=1:401] {

pl "dat.out" i 0 u 0:(column(k)) w l lw 3 t "explicito",  "dat.out" i 1 u 0:(column(k)) w l lw 3 t  "implicito" ,  "dat.out" i 2 u 0:(column(k)) w l lw 3  t "CN" , [0:100] u(x,k,-0.25,100) t "Solucion analitica"  w p   pt 7 


set title sprintf("%i",k)
#pause(0.1)

}






