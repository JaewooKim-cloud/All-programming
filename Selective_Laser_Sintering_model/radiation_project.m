%% This is FDM calculation of Radiation transfer course project
Nx=60;
Ny=20;
Nz=20;
Nt=10;
T=zeros(Nx,Ny,Nz,Nt);
H=zeros(Nx,Ny,Nz,Nt);
fi=zeros(Nx,Ny,Nz,Nt);
 
X_minus=0.0003;
X_plus=0.0003;
Y=0.0002;
Z=0.0002;
x=zeros(Nx,1);
y=zeros(Ny,1);
z=zeros(Nz,1);
  
delx=(X_minus+X_plus)/Nx;
dely=Y/Ny;
delz=Z/Nz;
delt=10^(-6);
for i=0:Nx-1
   x(i+1)=-X_minus+delx/2+i*delx; 
end
for i=0:Ny-1
   y(i+1)=dely/2+i*dely; 
end
for i=0:Nz-1
   z(i+1)=delz/2+i*delz; 
end
 
Tm=1700;
Hm=2.18*10^9;
Cs=4.25*10^6;
Cl=5.95*10^6;
 
rho=0.7;
a=sqrt(1-rho);
D1=20*10^(-6);
L=50*10^(-6);
lamda=2;
beta=lamda/L;
v=20*10^(-2);
R=60*10^(-6); %%FWHM
W=30;
%%Initial conditions
T(:,:,:,1)=300;
H(:,:,:,1)=300*Cs;
for l=1:Nz
    if z(l)<L
            fi(:,:,l,1)=0;
    else
            fi(:,:,l,1)=1;
    end
end
%%define bdry conditions of fi
 
for l=1:Nz
    if z(l)<L
        fi(1,:,l,:)=0;
    else
        fi(1,:,l,:)=1;
    end
end
 
 
%%solver using minmod slope limiter
 
for t1=1:Nt
  %%from bdry conditions of fi, construct all the values of fi at time t1.
  
  for t2=1:1000
    for i=2:Nx-1 %%the reason of starting with 2 is that 1 stands for the bonundary.
      for j=2:Ny-1
        for l=2:Nz-1
            sigma=1/delx*minmod(fi(i+1,j,l,t1)-fi(i,j,l,t1),fi(i,j,l,t1)-fi(i-1,j,l,t1));
            fi(i,j,l,t1+1)=fi(i,j,l,t1)+v*delt/delx*sigma*(x(i+1)-x(i-1))/2;
            sigma1=1/delx*minmod(H(i+1,j,l,t1)-H(i,j,l,t1),H(i,j,l,t1)-H(i-1,j,l,t1));
            H(i,j,l,t1+1)=H(i,j,l,t1)+v*delt/delx*sigma1*(x(i+1)-x(i-1))/2+delt/delx*((k(fi(i+1,j,l,t1))+k(fi(i,j,l,t1)))/2*(T(i+1,j,l,t1)-T(i,j,l,t1))/delx-(k(fi(i,j,l,t1))+k(fi(i-1,j,l,t1)))/2*(T(i,j,l,t1)-T(i-1,j,l,t1))/delx)+delt/dely*((k(fi(i,j+1,l,t1))+k(fi(i,j,l,t1)))/2*(T(i,j+1,l,t1)-T(i,j,l,t1))/dely-(k(fi(i,j,l,t1))+k(fi(i,j-1,l,t1)))/2*(T(i,j,l,t1)-T(i,j-1,l,t1))/dely)+delt/delz*((k(fi(i,j,l+1,t1))+k(fi(i,j,l,t1)))/2*(T(i,j,l+1,t1)-T(i,j,l,t1))/delz-(k(fi(i,j,l,t1))+k(fi(i,j,l-1,t1)))/2*(T(i,j,l,t1)-T(i,j,l-1,t1))/delz)+delt*U(sqrt(x(i)^2+y(j)^2),beta*z(l),W,R,a,D(a,rho,lamda),rho,lamda,beta);
        
        end
      end
    end
    
      %%apply boundary condition(Neumann boundary condition of flux = 0)
    
  H(:,1,:,t1+1)=H(:,2,:,t1+1);
  H(:,Ny,:,t1+1)=H(:,Ny-1,:,t1+1);
  
  H(:,:,1,t1+1)=H(:,:,2,t1+1);
  H(:,:,Nz,t1+1)=H(:,:,Nz-1,t1+1);
  
  
  H(1,:,:,t1+1)=H(2,:,:,t1+1);
  H(Nx,:,:,t1+1)=H(Nx-1,:,:,t1+1);
 
  %%Define  temperature in each node
   for i=1:Nx
      for j=1:Ny
          for l=1:Nz
              if H(i,j,l,t1+1)<Cs*Tm
                  T(i,j,l,t1+1)=H(i,j,l,t1+1)/Cs;
              elseif Cs*Tm<H(i,j,l,t1+1) && H(i,j,l,t1+1)<Cs*Tm+Hm
                  T(i,j,l,t1+1)=Tm;
              else
                  T(i,j,l,t1+1)=Tm+(H(i,j,l,t1+1)-Cs*Tm-Hm)/Cl;
              end
              %%melting condition of fi
              if T(i,j,l,t1+1)>Tm
                  fi(i,j,l,t1+1)=1;
              end
          end
      end
   end
   H(:,:,:,t1)=H(:,:,:,t1+1);
   T(:,:,:,t1)=T(:,:,:,t1+1);
   fi(:,:,:,t1)=fi(:,:,:,t1+1);
  end
end
 
 
function answer=D(a,rho,lamda)
    answer=(1-a)*(1-a-rho*(1+a))*exp(-2*a*lamda)-(1+a)*(1+a-rho*(1-a))*exp(2*a*lamda);
end
 
function answer=q(zeta,a,D,rho,lamda)
   answer=rho*a/(D*(4*rho-3))*((1-rho^2)*exp(-lamda)*((1-a)*exp(-2*a*zeta)+(1+a)*exp(2*a*zeta))-(3+rho*exp(-2*lamda))*((1+a-rho*(1-a))*exp(2*a*(lamda-zeta))+(1-a-rho*(1+a))*exp(2*a*(zeta-lamda))))-3*(1-rho)*(exp(-zeta)-rho*exp(zeta-2*lamda))/(4*rho-3);
    
end
 
function answer=dqdz(zeta,a,D,rho,lamda)
    if zeta<=lamda
        answer=(a*rho*((2*a*exp(-2*a*(lamda - zeta))*(a + rho*(a + 1) - 1) + 2*a*exp(2*a*(lamda - zeta))*(a + rho*(a - 1) + 1))*(rho*exp(-2*lamda) + 3) - exp(-lamda)*(rho^2 - 1)*(2*a*exp(-2*a*zeta)*(a - 1) + 2*a*exp(2*a*zeta)*(a + 1))))/(D*(4*rho - 3)) - ((3*rho - 3)*(exp(-zeta) + rho*exp(zeta - 2*lamda)))/(4*rho - 3);
    else
        answer = 0;
    end
end
 
function answer=Q0(r,W,R)
   
    if r<R
        Qm=(W/(R^2))*(3/pi);
        answer=Qm*((1-r/R)^2)*(1+r/R)^2;
    else
        answer=0;
    end
end
 
function answer=k(x)
Kd=20;
Kp=0.3;
    answer=Kp+(Kd-Kp)*x;
end
function answer=U(r,zeta,W,R,a,D,rho,lamda,beta)
    answer=-beta*Q0(r,W,R)*dqdz(zeta,a,D,rho,lamda);
end
 
function answer=minmod(a,b)
    answer=0.5*(sign(a)+sign(b))*min(abs(a),abs(b));
end
