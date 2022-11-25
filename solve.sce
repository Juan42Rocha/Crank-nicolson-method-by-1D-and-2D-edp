// codigo para comparar diferentes esquemas en 1d
// como euler explicito, euler implicio y crank nicolson
//
//  tomando como edp 1d 
//    dp u        dp u
//   ------ =  a ----- 
//    dp t        dp x


// const
a = -0.25


// puntos de la malla espacial 
nx=100;

//paso en espacio
L=10;//%pi;
dx=L/(nx-1);

// paso de tiempo
dt=0.1;
// pasos en tiempo
nt=400


x= dx*[0:nx-1]';

//const del sist
al =  a*dt/dx/2.0
be =  al/2.0


// para el caso euler explicito

u1=zeros(nx,nt+1);

// condicion inical 
u1(:,1)=exp(-3.*(x-L/2.0).^2);//sin(dx*[0:nx-1]);

u1(1,:)=0;
u1(nx,:)=0;


for n=1:nt

   for j=2:nx-1
   	u1(j,n+1) = al*(u1(j+1,n)-u1(j-1,n))+u1(j,n);
   end
    
end

// printeo
for j=1:nx
   printf("%f  ",u1(j,:)')
   printf("\n")
end
   printf("\n")
   printf("\n")

// fin de euler explicito 



// Preparativos para el euler implicito
// se tiene un sitema de ecuaciones de la forma 
//         A2*u = b2
// cuya solucion es 
//  
//        u = inv(u)*b2



u2=zeros(nx,nt+1);

// condicion inical 
u2(:,1)=exp(-3.*(x-L/2).^2);//sin(dx*[0:nx-1]);

u2(1,:)=0;
u2(nx,:)=0;


A2=zeros(nx,nx);

b2=zeros(nx,1)

t2= [al,1,-al];


// crea matrices para el sistema de eq

A2(1,[1,2])= t2([2,3]);
A2(nx,[nx-1,nx])= [al,1];

for j=2:nx-1
  b2(j,1) = u2(j,1)+al*(u2(j+1,1)-u2(j-1,1));
  
  A2(j,[j-1,j,j+1])=t2;

end


// Calcula inveraso y aplica condicion de frontera
iA2=inv(A2)
iA2([1,nx],:)=0;


// resuelve 
for n=1:nt
  u2(:,n+1)=iA2*u2(:,n);

end


// printeo
for j=1:nx
   printf("%f  ",u2(j,:)')
   printf("\n")
end
   printf("\n")
   printf("\n")


// Fin euler implicito



// Preparativos para el CN
// se tiene un sitema de ecuaciones de la forma 
//         A3*u = b3
// cuya solucion es 
//  
//        u = inv(u)*b3


u3=zeros(nx,nt+1);

// condicion inical 
u3(:,1)=exp(-3.*(x-L/2).^2);//sin(dx*[0:nx-1]);

u3(1,:)=0;
u3(nx,:)=0;


A3=zeros(nx,nx);

b3=zeros(nx,1)

t3= [be,1,-be];


// crea matrices para el sistema de eq

A3(1,[2,3])= [1,be];
A3(nx,[nx-1,nx])= [be,1];

for j=2:nx-1
  b3(j,1) = be*(u1(j+1,1)-u1(j-1,1))+u1(j,1);
  
  A3(j,[j-1,j,j+1])=t3;

end


// Calcula inveraso y aplica condicion de frontera
iA3=inv(A3)
iA3([1,nx],:)=0;



// resuelve 
for n=1:nt
  u3(:,n+1)=iA3*b3;
  
  for j=2:nx-1
    b3(j,1) = be*(u3(j+1,n+1)-u3(j-1,n+1))+u3(j,n+1);
  end

end




// printeo
for j=1:nx
   printf("%f  ",u3(j,:)')
   printf("\n")
end
   printf("\n")
   printf("\n")












exit
