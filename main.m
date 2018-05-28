clear,clc;
%基站的真实坐标
% x0=[-2364337.4414;4870285.6211;-3360809.6724];
% x1=[-2364335.4220;4870281.4604;-3360816.7056];
x0 = [-2364337.4414;4870285.6211;-3360809.6724];
x1 = [-2364333.5346;4870287.3393;-3360809.5251];
wrong=0;
correct=0;
load('data16p5_27cut0cutb.mat')       %得到卫星的星历文件nav
%得到基站接收机星历文件base
%得到用户接收  机星历文件obs

a=6378137;
f = 1/298.257222101;
e=sqrt(f*(2-f));
lambda=atan2(x0(2),x0(1));
phi=0;
for i=1:4
    N=a/sqrt(1-e^2*(sin(phi))^2);
    p=sqrt(x0(1)^2+x0(2)^2);
    h=p/cos(phi)-N;
    phi=atan(x0(3)/(p*(1-(N/(N+h))*e^2)));
end
S=[-sin(lambda) cos(lambda) 0;...
    -sin(phi)*cos(lambda) -sin(phi)*sin(lambda) cos(phi);...
    cos(phi)*cos(lambda) cos(phi)*sin(lambda) sin(phi)];

h=waitbar(0,'请等待...');
group=2880;
for m=1:group
    [BASEprn,EWLb,basenum,xbs,ybs,zbs,rb,Base]=BASESateposAndC1c(navdata,basedata,x0,S,m);
    
    [OBSprn,EWLu,obsnum,Utheta,xus,yus,zus,ru,User]=BASESateposAndC1c(navdata,obsdata,x0,S,m);
     
%     [OBSprn,EWLu,obsnum,Utheta,xus,yus,zus,ru,User]=OBSSateposAndC1c(navdata,obsdata,x1,S,m);
    [Ix,Iy,Iz,prn,satnum,maxnum,fw,pw,EWL,FC,pc]=SD(OBSprn,BASEprn,obsnum,xus,yus,zus,Utheta,basenum,x0,EWLb,EWLu,Base,User);
    [N,d,Qxn,Qn,NEWL,NWL1,N1,N2,N3,DOP]=DD(satnum,Ix,Iy,Iz,maxnum,fw,pw,EWL,FC,pc);
    
    
    DOP1.HDOP(m) = sqrt(DOP(1,1)+DOP(2,2));
    DOP1.VDOP(m) = sqrt(DOP(3,3));
    DOP1.PDOP(m) = sqrt(DOP(1,1)+DOP(2,2)+DOP(3,3));
    %%       修正后的坐标值
    df=d-Qxn/Qn*(N-N1);
    posf.x(m) = df(1)+x0(1);
    posf.y(m) = df(2)+x0(2);
    posf.z(m) = df(3)+x0(3);
    %% 求取CEP
    %用户站的真实坐标
    x1=[-2364337.4414;4870285.6211;-3360809.6724];
    aberration.dx(m) = posf.x(m)-x1(1);
    aberration.dy(m) = posf.y(m)-x1(2);
    aberration.dz(m) = posf.z(m)-x1(3);
    env=S*[aberration.dx(m);aberration.dy(m);aberration.dz(m)];
    CEP.L(m) = sqrt(env(1)^2+env(2)^2) ;
    CEP.H(m) = env(3);
    
    st.num(m) = satnum;
    string = ['运行中  ',num2str(floor(m/group*100)),'%'];
    waitbar(m/group,h,string);
end
close(h);
cepl=sort(CEP.L,'ascend');
CEPL95 = cepl(round(2880*0.95));
fprintf('CEPL95为%.8f\n',CEPL95);
ceph=sort(CEP.H,'ascend');
CEPH95 = ceph(round(2880*0.95));
fprintf('CEPH95为%.8f\n',CEPH95);


figure(1)
plot(st.num,'.black')

figure(2)
plot(aberration.dx,'.blue')
hold on;
plot(aberration.dy,'.green')
hold on;
plot(aberration.dz,'.red')

figure(3)
plot(DOP1.HDOP,'blue')
hold on;
plot(DOP1.VDOP,'green')
hold on;
plot(DOP1.PDOP,'red')








