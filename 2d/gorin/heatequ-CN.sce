//
// We use the simplest discretization in space (5-point stencil 
// finite-difference formula) and time (Euler formula). Not here, instead we
// try the semi-implicite Crank-Nicholson scheme which involves super-matrix
// inversion.

// The grid in the coordinates (x,y). We require an equal grid-point spacing
// "h" in x- and y-direction. The grid area, however, may be rectangular. But
// in this first attempt, all boundary points must be grid points.

// Implementation of the Laplace operator with the 5-point stencil formula:
// [Laplace u]_{ij} 
//     = h^{-2} ( u_{i-1,j} + u_{i+1,j} + u_{i,j-1} + u_{i,j+1} - 4 u_{i,j}
//
 function [u2]= laplace5(u1)
   u2= zeros(u1);
   for i= 2:nx
     for j= 2:ny
       u2(i,j)= u1(i-1,j)+u1(i+1,j)+u1(i,j-1)+u1(i,j+1) - 4*u1(i,j);
     end
   end
 endfunction


// Convert internal part of the solution u into super-vector
//
 function [uv]= supervec(u,nx,ny)
   uv= zeros((nx-1)*(ny-1),1);
   for j=2:nx
     for l=2:ny
       a= (nx-1)*(l-2) + j-1;
       uv(a) = u(j,l);
     end
   end
 endfunction

 function [u]= isupervec(uv,nx,ny)   // do the inverse
   u= zeros(nx+1,ny+1);
   for j=2:nx
     for l=2:ny
       a= (nx-1)*(l-2) + j-1;
       u(j,l)= uv(a);
     end
   end
 endfunction


// Construct the matrix representation of X = del_x^2 + del_y^2 for the
// internal part of the solution:
//
 function [X]= laplace_cn(nx,ny)
   dim= (nx-1)*(ny-1);
   X= zeros(dim,dim);
   for j=2:nx
     for l=2:ny
       a= (nx-1)*(l-2) + j-1;
       X(a,a)= -4.0;
       if (j-1 > 1) then
         b= (nx-1)*(l-2) + j-2;
         X(a,b)= 1.0;
       end
       if (j+1 < nx+1) then
         b= (nx-1)*(l-2) + j;
         X(a,b)= 1.0;
       end
       if (l-1 > 1) then
         b= (nx-1)*(l-3) + j-1;
         X(a,b)= 1.0;
       end
       if (l+1 < ny+1) then
         b= (nx-1)*(l-1) + j-1;
         X(a,b)= 1.0;
       end
     end
   end
 endfunction


 h= 0.2   // spacing between grid points -- equal in x- and y-direction

 xcent= 0.0;       // This defines (nx+1) times (ny+1) grid points, and thereby
 nx= 40;           // nx times ny little squares.
 ycent= 0.0;       // Plotting a function u(x,y) with gnuplot and pm3d, i.e.
 ny= 25;           // splot 'test.ou' w pm3d color-codes the value of u inside
                   // each square -- interpolating the grid values.

 xa= linspace(-0.5*h*nx,0.5*h*nx,nx+1);
 ya= linspace(-0.5*h*ny,0.5*h*ny,ny+1);

// The initial state: we assume a Gaussian profile, which is separable in x-
// and y-coordinates: u0(x,y)= ux0(x)*uy0(y)

 alpha = 1.0;   // This is the diffusion constant in the heat equation

 T0 =  20.0;
 T1 = 120.0;
 x0 = -1.0;
 y0 =  0.0;
 sigx = 0.3;
 sigy = 0.6;

 deff('z= ux0(x)', 'z= (T1-T0)*exp( -0.5*(x-x0)^2/sigx^2)');
 deff('z= uy0(y)', 'z= exp( -0.5*(y-y0)^2/sigy^2)');

 u0= zeros(nx+1,ny+1);
 for i=1:nx+1
   for j=1:ny+1
     u0(i,j)= T0 + ux0(xa(i))*uy0(ya(j));
   end
 end

 ub0= zeros(nx+1,ny+1);
 ub0(1,:) = u0(1,:);
 ub0(nx+1,:) = u0(nx+1,:);
 ub0(:,1) = u0(:,1);
 ub0(:,ny+1) = u0(:,ny+1);

// The evolution in time. For that we need a grid in time.
// Well not really, we just need to know the time-increment for each Euler
// step, and the number of time-steps to be applied.

 tmax = 0.5;
 dt= 0.02;
 nt= round(tmax/dt);

 X= laplace_cn(nx,ny);
 L= eye(X) - alpha * dt/(2.0*h^2) * X;

 iL= inv(L);
 
// apply (del_x^2 + del_y^2) to the boundary of u(x,y)
// and convert the result to a supervector
//
 bv0= alpha * dt/h^2 * supervec(laplace5(ub0),nx,ny);

 u= u0;
 t= 0.0;

// pause

 uv= supervec(u,nx,ny);
 for it=1:nt
   t= t + dt;
   mfprintf(0," %g\n", t);

   bv= uv + alpha * dt/(2.0*h^2) * X*uv + bv0;
   uv= iL*bv;

//   u= u + alpha * dt/h^2 * laplace5(u);    Euler-scheme
 end

 u= ub0 + isupervec(uv,nx,ny); 

 for i=1:ny+1
   for j=1:nx+1
     printf(" %g   %g   %g\n", xa(j),ya(i),u(j,i));
   end
   printf("\n");
 end

 quit

