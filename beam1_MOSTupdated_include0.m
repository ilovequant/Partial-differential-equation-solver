%note: D = lamda^4 = (mb*omega^2)/(E*I) , so natural frequency omega =
%(E*I/mb)^(1/2)*(D^(1/4))^2
%parameters
%This program is to solve the eigenfrequencies and eigenvalues of the
%cantilever beam vibration
E = 2.8e9;
%E=2.8e9;
%E = 1;
%beam thickness
t = 5e-6;
%h1 = 1e-5;
%beam width
b = 0.5*(500e-6);
%upper cavity
h2 = 1*(500e-6);
%lower cavity
h1 =0.5*(500e-6);
%moment of inertia
%I = t*h1^3/12;
I=b*t^3/12;
%I = 1;
n = 500;
N = 500;
%L =1;
%beam length
%L = 2.5e-4;
L=0.5*(500e-6);
a=1*(500e-6);
%external force
G=1e-6;
%mesh interval
h = L/n;
ca = 1/h^4;
%mb = 1289*t*b;
mb=8e-6;
%liquid density
rho =1000;
%liquid viscosity
mu=(1e-3)*42.7;
%residual stress
%NS = 8e6*t*h1;
NS=8.2e-2;
cc = NS/(E*I*h^2);
%test
% L=0.5;
% a=1;
% E=1;
% b=0.5;
% I=1;
% h=L/n;
% ca=1/h^4;
% mb=1;
%create matrix A
A = zeros(N,N);
%A = zeros(N,N);
A(1,1) = (7*ca+2*cc);
A(N-1,N-2) = (-27*ca/7-cc);
A(N-1,N-1) = (33*ca/7+2*cc);
A(N-1,N) = (-13*ca/7-cc);
A(N,N-2) = (12*ca/7-cc/7);
A(N,N-1) = (-24*ca/7+2*cc/7);
A(N,N) = (12*ca/7-cc/7);

for i=1:N-2
    A(i,i+1) = -4*ca-cc;
end

for j=3:N-1
    A(j,j-2) = ca;
end

for j=1:N-2
    A(j,j+2) = ca;
end

for k=2:N-2
    A(k,k) = 6*ca+2*cc;
end

for m=2:N-2
    A(m,m-1) = -4*ca-cc;
end

[V D] = eig(A);
%initialization 
%vibration mode
k1= 12;
%fluid mode(fourier series mode)
k2= 12;
k3= 12;
s1 =zeros(1,k1);
s1 =zeros(1,k1);
omegan =zeros(1,k1);
f=zeros(1,k1);
c1=zeros(1,k1);
c2=zeros(1,k1);
c3=zeros(1,k1);
c4=zeros(1,k1);
phi = zeros(N+1,k1);
phi1 = zeros(N+1,k1);
phi2 = zeros(N+1,k1);
phi3 = zeros(N+1,k1);
phi4 = zeros(N+1,k1);
for mn =1:k1
s1(mn) = L*(NS/(2*E*I)+(NS^2/(2*E*I)^2 + D(N+1-mn,N+1-mn))^(1/2))^(1/2);
s2(mn) = L*((NS^2/(2*E*I)^2 + D(N+1-mn,N+1-mn))^(1/2)-NS/(2*E*I))^(1/2);
omegan(mn) = (E*I/mb)^(1/2)*D(N+1-mn,N+1-mn)^(1/2);
%omegatest=(1.8751)^2*sqrt(E*I/(mb*L^4));
%ftest=omegatest/(2*pi);
f(mn)=omegan(mn)/(2*pi);
%E = [0;0];
%H = C\E;
%coefficient for the mode shape of cantilever beam phi =
%c1*cosh(s1*x)+c2*sinh(s1*x)+c3*cos(s2*x)+c4*sin(s2*x)
% c1 = 1;
% c3= -c1;
% c4 = -c1*C(1,1)/C(1,2);
% c2 = -s2*c4/s1;
c1(mn) = -(s1(mn)^2*sinh(s1(mn))+s1(mn)*s2(mn)*sin(s2(mn)))/(s1(mn)^2*cosh(s1(mn))+s2(mn)^2*cos(s2(mn)));
c2(mn) = 1;
c3(mn) = (s1(mn)^2*sinh(s1(mn))+s1(mn)*s2(mn)*sin(s2(mn)))/(s1(mn)^2*cosh(s1(mn))+s2(mn)^2*cos(s2(mn)));
c4(mn)= -s1(mn)/s2(mn);
x= [0:h:L];
s1(mn) = s1(mn)/L;
s2(mn) = s2(mn)/L;
%C = [s1^2*cosh(s1*L)+s2^2*cos(s2*L) -s1*s2*sinh(s1*L)-s2^2*sin(s2*L);s1^3*sinh(s1*L)-s2^3*sin(s2*L)-NS*s1*sinh(s1*L)/(E*I)-NS*s2*sin(s2*L)/(E*I) -s1^2*s2*cosh(s1*L)-s2^3*cos(s2*L)+NS*s2*cosh(s1*L)/(E*I)-NS*s2*cos(s2*L)/(E*I)];
%disp(det(C));
phi(:,mn)= c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x);
phi1(:,mn) = c1(mn)*s1(mn)*sinh(s1(mn).*x)+c2(mn)*s1(mn)*cosh(s1(mn).*x)-c3(mn)*s2(mn)*sin(s2(mn).*x)+c4(mn)*s2(mn)*cos(s2(mn).*x);
phi2(:,mn) = c1(mn)*(s1(mn))^2*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^2*sinh(s1(mn).*x)-c3(mn)*(s2(mn))^2*cos(s2(mn).*x)-c4(mn)*(s2(mn))^2*sin(s2(mn).*x);
phi3(:,mn) = c1(mn)*(s1(mn))^3*sinh(s1(mn).*x)+c2(mn)*(s1(mn))^3*cosh(s1(mn).*x)+c3(mn)*(s2(mn))^3*sin(s2(mn).*x)-c4(mn)*(s2(mn))^3*cos(s2(mn).*x);
phi4(:,mn) = c1(mn)*(s1(mn))^4*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^4*sinh(s1(mn).*x)+c3(mn).*(s2(mn))^4*cos(s2(mn).*x)+c4(mn)*(s2(mn))^4*sin(s2(mn).*x);
%plot(x,phi) 
end
A= zeros(1,k2);
B= zeros(1,k3);
%mode of fluid flow
%p =zeros(3*N+1,N+1);
qa =zeros(1,k1);
qb =zeros(1,k1);
qc =zeros(1,k1);
qd =zeros(1,k1);
for mn=1:k1
%numerical integral
%funa = interal of phi*phi'''' from 0 to L
%funa1 = phi*phi''', funa2 = phi''*phi', funa = (phi'')^2
funa = @(x)(c1(mn)*(s1(mn))^2*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^2*sinh(s1(mn).*x)-c3(mn)*(s2(mn))^2*cos(s2(mn).*x)-c4(mn)*(s2(mn))^2*sin(s2(mn).*x))...
    .*(c1(mn)*(s1(mn))^2*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^2*sinh(s1(mn).*x)-c3(mn)*(s2(mn))^2*cos(s2(mn).*x)-c4(mn)*(s2(mn))^2*sin(s2(mn).*x));
funa1 = @(x)(c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x))...
    .*(c1(mn)*(s1(mn))^3*sinh(s1(mn).*x)+c2(mn)*(s1(mn))^3*cosh(s1(mn).*x)+c3(mn)*(s2(mn))^3*sin(s2(mn).*x)-c4(mn)*(s2(mn))^3*cos(s2(mn).*x));
funa2 = @(x)(c1(mn)*(s1(mn))^2*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^2*sinh(s1(mn).*x)-c3(mn)*(s2(mn))^2*cos(s2(mn).*x)-c4(mn)*(s2(mn))^2*sin(s2(mn).*x))...
    .*(c1(mn)*s1(mn)*sinh(s1(mn).*x)+c2(mn)*s1(mn)*cosh(s1(mn).*x)-c3(mn)*s2(mn)*sin(s2(mn).*x)+c4(mn)*s2(mn)*cos(s2(mn).*x));
qa(mn) = quad(funa,0,L,0.1)+funa1(L)-funa1(0)-(funa2(L)-funa2(0));
%funb = interal of phi*phi from 0 to L
funb = @(x)(c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x)).^2;
qb(mn) = quad(funb,0,L,0.1);
%func = interal of phi*phi'' from 0 to L
func = @(x)(c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x))...
             .*(c1(mn)*(s1(mn))^2*cosh(s1(mn).*x)+c2(mn)*(s1(mn))^2*sinh(s1(mn).*x)-c3(mn)*(s2(mn))^2*cos(s2(mn).*x)-c4(mn)*(s2(mn))^2*sin(s2(mn).*x));
qc(mn) = quad(func,0,L,0.1);
%fund = interal of phi from 0 to L
fund = @(x)(c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x));
qd(mn) = quad(fund,0,L,0.1);
end
%falpha =zeros(k2,k1);
%fsigma =zeros(1,k2);
qalpha =zeros(k2,k1);
qsigma =zeros(k2,k2);
qbeta =zeros(k2,k1);
% qc2 =zeros(k1,k1);
% %integral of phi(j)*phi''(i)
% for mn=1:k1
%     for mc=1:k1
%        func = @(x)(c1(mn)*cosh(s1(mn).*x)+c2(mn)*sinh(s1(mn).*x)+c3(mn)*cos(s2(mn).*x)+c4(mn)*sin(s2(mn).*x))...
%              .*(c1(mc)*(s1(mc))^2*cosh(s1(mc).*x)+c2(mc)*(s1(mc))^2*sinh(s1(mc).*x)-c3(mc)*(s2(mc))^2*cos(s2(mc).*x)-c4(mc)*(s2(mc))^2*sin(s2(mc).*x));
%        qc2(mn,mc) = quad(func,0,L,0.1);
%     end
% end
for nj= 1:k2
    for it=1:k1
      falpha = @(x)(c1(it)*cosh(s1(it).*x)+c2(it)*sinh(s1(it).*x)+c3(it)*cos(s2(it).*x)+c4(it)*sin(s2(it).*x)).*cos(nj*pi*x/a);
      qalpha(nj,it) = quad(falpha,0,L,0.1);
      fbeta =@(x)c1(it)*(s1(it))^2*cosh(s1(it).*x)+c2(it)*(s1(it))^2*sinh(s1(it).*x)-c3(it)*(s2(it))^2*cos(s2(it).*x)-c4(it)*(s2(it))^2*sin(s2(it).*x)...
          .*cos(nj*pi*x/a);
      qbeta(nj,it) = quad(fbeta,0,L,0.1);
    end
end
for nj=1:k2
    for nj1=1:k2
       fsigma = @(x)cos(nj*pi*x/a).*cos(nj1*pi*x/a);
       qsigma(nj,nj1) = quad(fsigma,L,a,0.1);
    end
end
%excitation frequency
  df=zeros(50,1);
for ome = 1:50
  omega=ome*0.5*2*pi;
  msab =zeros(k1+k2+k3,k1+k2+k3);
%input equation #1 into the matrix
  for ct1=1:k1
        msab(ct1,ct1)=E*I*qa(ct1)-mb*omega^2*qb(ct1)-NS*qc(ct1);
        for ct2=1:k2
            msab(ct1,ct2+k1)=-b*cosh(ct2*pi*h1/a)*qalpha(ct2,ct1);
        end
        for ct3=1:k3
            msab(ct1,k1+k2+ct3)=b*qalpha(ct3,ct1);
        end
        %msab(k1+ct1,ct1)=rho*omega^2*qb(ct1)+mu*(1i)*omega*qc(ct1);
        for ct4=1:k2
             msab(k1+ct4,ct1)=+rho*omega^2*qalpha(ct4,ct1)+mu*(1i)*omega*qbeta(ct4,ct1);
         end
        for ct5=1:k2
            msab(k1+ct5,k1+ct5)=-(a/2)*(ct5*pi/a)*sinh(ct5*pi*h1/a);
            for ct6=1:k3
              msab(k1+ct5,k1+k2+ct6)=-(ct6*pi/a)*tanh(ct6*pi*h2/a)*qsigma(ct5,ct6);
            end
        end
        %msab(k1+k2+ct1,ct1)=rho*omega^2*qb(ct1)+mu*(1i)*omega*qc(ct1);
        for ct7=1:k3
            msab(k1+k2+ct7,ct1)=+rho*omega^2*qalpha(ct7,ct1)+mu*(1i)*omega*qbeta(ct7,ct1);
        end
%             msab(k1+k2+ct1,ct7)=msab(k1+k2+ct1,ct7)+mu*(1i)*omega*qc(ct7);
%         end
        for ct8=1:k3
            for ct9=1:k2
              msab(k1+k2+ct8,k1+ct9)=(ct9*pi/a)*sinh(ct9*pi*h1/a)*qsigma(ct8,ct9);
            end
            msab(k1+k2+ct8,k1+k2+ct8)=(a/2)*(ct8*pi/a)*tanh(ct8*pi*h2/a);
        end
%         for ct9=1:k3
%             msab(k1+k2+ct1,k1+k2+ct9)=(ct9*pi/a)*tanh(ct9*pi*h2/a)*qalpha(ct9,ct1)-qsigma(ct9);
%         end 
  end
  bsab =zeros(k1+k2+k3,1);
  for it=1:k1
    bsab(it,1)=G*qd(it);
  end
  sab = zeros(k1+k2+k3,1);
  sab =msab\bsab;
%pressure field
  y1 =(0:h1/N:h1-h1/N)';
  x1 =(0:a/N:a)';
%pa = sab(2,1)*cos(pi*x1(:)./a)*cosh(pi*y1(:)'./a);
  y2 =(0:h2/(2*N):h2)';
  te=pi/omega;
%p(2*N+1:1:3*N+1,1:1:N+1) = sab(2,1)*cosh(pi*(h1-y1(:))./a)*cos(pi*x1(:)'./a);
   p=zeros(3*N+1,N+1);
  for njp=1:k3
    p(N+1:1:3*N+1,1:1:N+1)= p(N+1:1:3*N+1,1:1:N+1)+sab(k1+k2+njp,1)*(cosh(njp*pi*y2(:)./a)-tanh(njp*pi*h2./a)*sinh(njp*pi*y2(:)./a))*cos(njp*pi*x1(:)'./a);
  end
%plot(0:h2/(2*N):h2,p(N+1:1:3*N+1,1:1:N+1));
  for njp2=1:k2
   p(1:1:N,1:1:N+1)=p(1:1:N,1:1:N+1)+sab(k1+njp2,1)*cosh(njp2*pi*y1(:)./a)*cos(njp2*pi*x1(:)'./a);
  end
  W = zeros(1,N+1);
  for t=1:N/2+1
     diff(t)=p(N,N/2+t)-p(N+1,N/2+t);
  end
  %A0
  p(1:1:N,1:1:N+1)=p(1:1:N,1:1:N+1)-0.0016;
  for w1=1:k1
    W(1:1:N+1) = W(1:1:N+1)+sab(w1,1).*phi(:,w1)';
  end
  df(ome)=W(N+1);
%[X,Y]= meshgrid(0:a/N:a,0:h2/(2*N):h2);
%z=sab(3,1)*(cosh(pi*Y./a)-tanh(pi*h2)*sinh(pi*Y./a)).*cos(pi*X./a);
%contourf(X,Y,z);
%p(N+1:1:3*N+1,1:1:N+1) = abs(p(N+1:1:3*N+1,1:1:N+1)).*(-1);
%p(1:1:N+1,1:1:N+1) = abs(p(1:1:N+1,1:1:N+1));
end
 [C,z] = contourf(x1,0:h1/N:h2+h1,p,30);
 %clabel(C,z);
 colorbar;
 %hold on;
plot(1:1:50,abs(df));
