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

// nodos en la malla
nx= 40;
ny= 25;

// longuitud en x y tamanio de paso en x
//lx=20;
//d=lx/(nx-1);
d= 0.2; 

// Pasos en tiempo y tamanio de paso en tiempo
dt= 0.02;
nt= round(0.5/dt);

// Constante de difucion
D=1;


// algunas variables
nxy= (nx+1)*(ny+1);


al= D*dt/d^2;


// crea matriz de unos uno arriba y uno abajo de la
// diagonal
//

id=eye(nxy,nxy);

// matrices superiroes eh inferiores
on=ones(ny+1,ny+1);
t=triu(on,-1)-triu(on);
Iplyy=t+t';


on = ones(nx+1,nx+1);
t=triu(on,-1)-triu(on);
Iplxx=t+t';


Ixx=eye(nx+1,nx+1);
Iyy=eye(ny+1,ny+1);


// condicion inicilal

//x=0:d:(nx-1)*d; 
//y=0:d:(ny-1)*d;

 x= linspace(-0.5*d*nx,0.5*d*nx,nx+1);
 y= linspace(-0.5*d*ny,0.5*d*ny,ny+1);

// u(j,k,l) -> u(x,y,t)
[Y,X]=meshgrid(y,x);

u0=(100)*exp(-0.5.*(X+1).^2./0.3^2-0.5.*(Y-0.0).^2./0.6^2)+70;

// para obtener los indices de las fronteras
uf0=u0*0;
uf0([1,nx+1],:)=1;
uf0(:,[1,ny+1])=1;


k=0;
for l=1:ny+1
  for j=1:nx+1
    k=k+1;
    u3(k,1)=u0(j,l);
    xx(k,1)=X(j,l);
    yy(k,1)=Y(j,l);
 
    ufr(k,1)=uf0(j,l);
  end
end


// inicia CM


a=(1+al);
b=-al/4;


aa=(1-al);
bb=al/4;


M3r = kron(Iyy,a*Ixx+b*Iplxx)+kron(Iplyy,b*Ixx);
M3l = kron(Iyy,aa*Ixx+bb*Iplxx)+kron(Iplyy,bb*Ixx);


// aplica condicion de frontera
M3r(ufr==1,:)=id(ufr==1,:);
M3l(ufr==1,:)=id(ufr==1,:);

// calcula inverso
iM3r = inv(M3r);



M3 = iM3r*M3l;

// resuelve

for n=1:nt
  u3(:,n+1)=M3*u3(:,n)
 // u3(ufr==1,n+1)=u3(ufr==1,1);
 mprintf("%f",n/nt)
end

//printeo

for k=1:nxy
  printf("%f  %f  ",xx(k,1),yy(k,1)) 
   
  for n=1:nt+1
    printf("%f  ",u3(k,n))
   
  end
    printf("\n")

end


// Fin de Crank Nicolson  en un paso 


// Inicia Crank Nicilson en dos pasos
// usando un paso intermeido de dt/2






exit
