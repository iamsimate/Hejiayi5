clc; clear; format longE;

h  =                  1/128; % mesh size in [0,1]
xh =              (0:h:1)'; % nodes in [0,1]
N  =          length(xh)-1; % mesh number in [0,1]


Ahv  = gallery('tridiag', N-1,-1,2,-1);%v
Ahu  = gallery('tridiag', N ,-1,2,-1);%u

I1   = speye(N  ); % some identity matrix
I2  = speye(N-1 ); % some identity matrix
I3=sparse(N-1,N);
I4=sparse(N-1,N-1);

K=sparse(2*(N-1)*N+N^2,2*(N-1)*N+N^2);
%%�����Ҷ�����u1
ul=sparse(2*(N-1)*N+N^2,1);
for i =1:N-1
   ul(i*N)=2/(h^2); 
end
%%ѹ������
ul(2*(N-1)*N+1,1)=0.005;%Ϊp��һ����ֵ0.005

for i=1:N-1
   I3(i,i)=-1;
   I3(i,i+1)=1;
end

for i=1:N-2
I4(i,i)=-1;    
I4(i,i+1)=1;
end
I4(N-1,N-1)=-1;

K1 = kron(I2, Ahu) + kron(Ahv, I1); % difference matrix for u
K2 = kron(I1 ,Ahv) + kron(Ahu, I2); % difference matrix for v

for i =1:N-1
   K1(i*N,i*N)=5; 
   K1((i-1)*N+1,(i-1)*N+1)=5;
   K2(i*(N-1),i*(N-1))=5;
   K2((i-1)*(N-1)+1,(i-1)*(N-1)+1)=5;
end


B1=kron(I3,I1);
B2=kron(I1,I3);


P1=gallery('tridiag',N,-1,1,0);
A=kron(P1,I1);
A=A(:,1:(N-1)*N);


P2=gallery('tridiag',N,0,-1,1);
P2=P2(1:end-1,:);
P3=-P2';
B=kron(I1,P3);

K(1:(N-1)*N,1:(N-1)*N)=K1./(h^2);
K((N-1)*N+1:2*(N-1)*N,(N-1)*N+1:2*(N-1)*N)=K2./(h^2);
K(1:(N-1)*N,2*(N-1)*N+1:end)=B1./h;
K((N-1)*N+1:2*(N-1)*N,2*(N-1)*N+1:end)=B2./h;
K(2*(N-1)*N+1:end,1:(N-1)*N)=A./h;
K(2*(N-1)*N+1:end,(N-1)*N+1:2*(N-1)*N)=B./h;

K(2*(N-1)*N+1,:)=0;
K(2*(N-1)*N+1,2*(N-1)*N+1)=1;%��ѹ��P��һ����ֵ

%%��matlab�Դ�����ⷽ��
uh = K\ul; 

uhu=uh(1:(N-1)*N);%uֵ
uhv=uh(N*(N-1)+1:2*N*(N-1));%vֵ
uhp=uh(2*(N-1)*N+1:2*(N-1)*N+N^2);%pֵ

uu=zeros(N,N+1);%�洢u����ֵ,�������ұ߽�
uu(:,2:N)=reshape(uhu,N,N-1);
uv=zeros(N+1,N);%�洢v����ֵ���������±߽�
uv(2:N,:)=reshape(uhv,N-1,N);
up=reshape(uhp,N,N);%�洢p����ֵ

%����ƽ������u��v�ĵ��p�����
uu1=zeros(N);
uv1=zeros(N);
for i=1:N
    uu1(:,i)=(uu(:,i)+uu(:,i+1))./2;
    uv1(i,:)=(uv(i,:)+uv(i+1,:))./2;
end

%%��ͼ
x=[h/2:h:1-h/2];
y=[h/2:h:1-h/2];
[X,Y]=meshgrid(x,y);

figure(1)
pcolor(X,Y,up);
figure(2)
pcolor(X,Y,uu1);
figure(3)
pcolor(X,Y,uv1);
figure(4)
contourf(X,Y,uu1);
figure(5)
contourf(X,Y,uv1);
figure(6)
streamslice(uu1,uv1);