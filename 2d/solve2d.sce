// 
// Resolviendo caso 2d para la ecuacion de
// difusion con 3-4 metotodos diferentes
// el primero es euler explicito en el que 
// con diferencias finitas 
//
stacksize("max")


// Primero tratamos de visualizar las matrices 
// mapeando del "vector" de estados  nx x ny 
// al hipervector     nx*ny x 1
//

// Algunos parametos

D=1;

nx= 50;
ny= 60;

nt= 100;
dt= 0.01

lx=10;


// algunas variables

d=lx/(nx-1) 

al= D*dt/d;

// condicion inicilal

x=0:d:(nx-1)*d; 
y=0:d:(ny-1)*d;


[Y,X]=meshgrid(y,x);

u0=50*exp(-3*(X-0.5*d*nx).^2-3*(Y-0.3*d*ny).^2)+20;


uf0=u0*0;
uf0([1,nx],:)=1;
uf0(:,[1,ny])=1;




nxy= nx*ny;

// u(j,k,l) -> u(x,y,t)

// Pare euler explicito tenemos 


M1=zeros(nxy,nxy);

// crea matriz de unos uno arriba y uno abajo de la
// diagonal
//

id=eye(nxy,nxy);

on=ones(ny,ny);
t=triu(on,-1)-triu(on);

Iplyy=t+t';

on = ones(nx,nx);
t=triu(on,-1)-triu(on);

Iplxx=t+t';


Ixx=eye(nx,nx);
Iyy=eye(ny,ny);


// las matrices Izz es una matriz identidad de tamanio 
// nzxnz
//



//
// matriz para euler explicito
//

a=(1-2*al);
b=al/2;

u1=zeros(nx*ny,nt+1)
ufr=zeros(nx*ny,1)
xx=zeros(nx*ny,1)
yy=zeros(nx*ny,1)

//
// Calcula  vector u al super vector
// asi como las coordenadas
//

k=0;
for l=1:ny
  for j=1:nx
    k=k+1;
    u1(k,1)=u0(j,l);
    xx(k,1)=X(j,l);
    yy(k,1)=Y(j,l);
 
    ufr(k,1)=uf0(j,l);
  end
end

M1 = kron(Iyy,a*Ixx+b*Iplxx)+kron(Iplyy,b*Ixx);

M1(ufr==1,:)=id(ufr==1,:);
// resuelve 

for n=1:nt
  u1(:,n+1)=M1*u1(:,n)
 // u1(ufr==1,n+1)=u1(ufr==1,1);

end

// printea

for k=1:nx*ny
  printf("%f  %f  ",xx(k,1),yy(k,1)) 
   
  for n=1:nt+1
    printf("%f  ",u1(k,n))

  end
    printf("\n")

end


    printf("\n")
    printf("\n")
    printf("\n")


// inicia euler implicito 


u2=u1;

a=(1+2*al);
b=-al/2;

M2 = kron(Iyy,a*Ixx+b*Iplxx)+kron(Iplyy,b*Ixx);
M2(ufr==1,:)=id(ufr==1,:);

iM2 = inv(M2);

// resuelve

for n=1:nt
  u2(:,n+1)=iM2*u2(:,n)
 // u2(ufr==1,n+1)=u2(ufr==1,1);
  

end

//printeo
for k=1:nx*ny
  printf("%f  %f  ",xx(k,1),yy(k,1)) 
   
  for n=1:nt+1
    printf("%f  ",u2(k,n))

  end
    printf("\n")

end


    printf("\n")
    printf("\n")
    printf("\n")

// fin euler implicito


// inicia CM

u3=u1;


a=(1+al);
b=-al/4;


aa=(1-al);
bb=al/4;


M3r = kron(Iyy,a*Ixx+b*Iplxx)+kron(Iplyy,b*Ixx);
M3l = kron(Iyy,aa*Ixx+bb*Iplxx)+kron(Iplyy,bb*Ixx);

M3r(ufr==1,:)=id(ufr==1,:);
M3l(ufr==1,:)=id(ufr==1,:);

iM3r = inv(M3r);

M3 = iM3r*M3l;

// resuelve

for n=1:nt
  u3(:,n+1)=M3*u3(:,n)
 // u3(ufr==1,n+1)=u3(ufr==1,1);

end

//printeo

for k=1:nx*ny
  printf("%f  %f  ",xx(k,1),yy(k,1)) 
   
  for n=1:nt+1
    printf("%f  ",u3(k,n))

  end
    printf("\n")

end


// Fin de Crank Nicolson  en un paso 


// Inicia Crank Nicilson en dos pasos
// usando un paso intermeido de dt/2






//exit
