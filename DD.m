function[N,d,Qxn,Qn,NEWL,NWL1,N1,N2,N3,DOP]=DD(satnum,Ix,Iy,Iz,maxnum,fw,pw,EWL,FC,pc)
% ILS for AR and Solution for differenet combinations ambiguities

% *使用的Casading Rounding AR 方法中，由于建立的观测方程中默认使用的是-N这样的
%  与之前建立的观测方程中的不同，因此求得的模糊度均是和真实的模糊度相反*
% *式中的载波相位观测量是以米为单位的*
dnum=0;
c = 2.99792458e8;%光速
f1=1561.098e6;
f2=1207.14e6;
f3=1268.52e6;
lamda1 = c/f1;
lamda2 = c/f2;
lamda3 = c/f3;
% lamda = c/f1;        %波长λ
for id=1:(satnum-1)
    dnum=dnum+1;
    if(dnum==maxnum),dnum=dnum+1;end
    G(id,1)=-(Ix(dnum)-Ix(maxnum));
    G(id,2)=-(Iy(dnum)-Iy(maxnum));
    G(id,3)=-(Iz(dnum)-Iz(maxnum));
    %单频的伪距、载波双差计算量
    H1(id,1)=pc(dnum)-pc(maxnum);
    H2(id,1)=FC(dnum)-FC(maxnum);
    
    %% 双差后三频的伪距、载波观测量，1代表伪距，2代表载波
    H11(id,1)=EWL.pc1(dnum)-EWL.pc1(maxnum);
    H21(id,1)=EWL.FC1(dnum)-EWL.FC1(maxnum);
    
    H12(id,1)=EWL.pc2(dnum)-EWL.pc2(maxnum);
    H22(id,1)=EWL.FC2(dnum)-EWL.FC2(maxnum);
    
    H13(id,1)=EWL.pc3(dnum)-EWL.pc3(maxnum);
    H23(id,1)=EWL.FC3(dnum)-EWL.FC3(maxnum);
    
    
    %%  求解超宽巷组合的模糊度，为Cascading rounding方法的第一步
    % *式中HEWL为组合(0,-1,1)的载波相位观测量??(0,-1,1)*
    lamdaEWL = c/(-f2+f3);
    HEWL(id,1) = (-H22(id,1)*f2+H23(id,1)*f3)/(-f2+f3);
    NEWL(id,1) = round((H13(id,1)- HEWL(id,1))/lamdaEWL);
    
    %%  求解宽巷组合的模糊度，为Cascading rounding方法的第二步
    %  WL1、WL2分别为组合(1,-1,0)、(1,0,-1)
    %  式中HWL1、HWL2为组合(1,-1,0)、(1,0,-1)的载波相位观测量??(1,-1,0)、??(1,0，-1)
    lamdaWL1 = c/(f1-f2);
    HWL1(id,1) = (H21(id,1)*f1-H22(id,1)*f2)/(f1-f2);
    NWL1(id,1) = round((lamdaEWL*NEWL(id,1)+HEWL(id,1)- HWL1(id,1))/lamdaWL1);
    
    lamdaWL2 = c/(f1-f3);
    HWL2(id,1) = (H21(id,1)*f1-H23(id,1)*f3)/(f1-f3);
    NWL2(id,1) = round((lamdaEWL*NEWL(id,1)+HEWL(id,1)- HWL2(id,1))/lamdaWL2);
    
    %%  求解原始的模糊度，为Cascading rounding方法的第三步
    %  N1、N2、N3分别为三个频率上的最初始的模糊度
    N1(id,1) = round((lamdaWL1*NWL1(id,1)+HWL1(id,1)- H21(id,1))/lamda1);
    N2(id,1) = N1(id,1)-NWL1(id,1);
    N3(id,1) = N1(id,1)-NWL2(id,1);
    for j=1:satnum-1
        if(id==j)
            Q1(id,j) = 2*(pw(id)+pw(maxnum));%伪距加权
            Q2(id,j) = 2*(fw(id)+fw(maxnum));%载波相位
        else
            Q1(id,j) = 2*pw(maxnum);%伪距加权
            Q2(id,j) = 2*fw(maxnum);%载波相位
        end
    end
end

%% 将使用 Casading Rounding AR 方法得到的模糊度进行取反得到真实的模糊度
N1 = -N1;
N2 = -N2;
N3 = -N3;


B=eye(satnum-1);
b=B*lamda1;
G1=zeros(satnum-1);
Q=[Q1,G1;G1,Q2];
C=inv(Q);
A=[G,G1;G,b];
H=[H11;H21];
X=(A'*A)\A'*H;
d = X(1:3);
N = X(4:satnum-1+3);
Qx = inv(A'*C*A);
Qn = Qx(4:satnum-1+3,4:satnum-1+3);   %模糊度N的协方差矩阵
Qxn = Qx(1:3,4:satnum-1+3);           %基线向量与模糊度N之间的相关系数阵


DOP = (G'*G)';


end