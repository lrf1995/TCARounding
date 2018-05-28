function[BASEprn,EWLb,basenum,xbs,ybs,zbs,rb,Base]=BASESateposAndC1c(navdata,basedata,x0,S,m)
% 计算卫星的坐标及接收机伪距

a1=m;
basenum=0;
match = 0;%判断星历是否匹配

for a2=1:basedata.epoch(a1).gpsobs
    for a3=1:length(navdata.gps)
        if(basedata.epoch(a1).gps(a2).prn==navdata.gps(a3).prn)
            tk = (basedata.epoch(a1).gps(a2).bdst(1)-navdata.gps(a3).bdst(1))*604800 + ...
                basedata.epoch(a1).gps(a2).bdst(2)-navdata.gps(a3).toe - basedata.epoch(a1).gps(a2).C2I/299792458;
            if(tk > 302400)
                tk = tk - 604800;
            elseif (tk<-302400)
                tk = tk + 604800;
            end
            if abs(basedata.epoch(a1).gps(a2).gpst-navdata.gps(a3).gpst)<7200,break;end
        end
        
        %% 无匹配导航文件
        if(a3 == length(navdata.gps)),match = 1;end
    end
    if(match==1)
        match=0;
        continue;
    end
    
    num = navdata.gps(a3).prn;
    toe=navdata.gps(a3).toe;
    t=basedata.epoch(a1).gps(a2).gpst-basedata.epoch(a1).gps(a2).C2I/299792458;
    as=(navdata.gps(a3).sqrtas)^2;
    es=navdata.gps(a3).es;
    io=navdata.gps(a3).io;
    OMGAo=navdata.gps(a3).OMGAo;
    w=navdata.gps(a3).w;
    Mo=navdata.gps(a3).Mo;
    deltn=navdata.gps(a3).deltn;
    dti=navdata.gps(a3).dti;
    dtOMGA=navdata.gps(a3).dtOMGA;
    Cuc=navdata.gps(a3).Cuc;
    Cus=navdata.gps(a3).Cus;
    Crc=navdata.gps(a3).Crc;
    Crs=navdata.gps(a3).Crs;
    Cic=navdata.gps(a3). Cic;
    Cis=navdata.gps(a3).Cis;
    dtOMGAe = 7.2921150e-5;%地球自转角速度
    GM=3.986004418e+14;
    
    %  1.计算归化时间tk
    
    %  2.计算卫星的平均角速度
    no=sqrt(GM/(as^3));
    n=no+deltn;
    
    %  3.计算信号发射时的平近角Mk
    Mk=Mo+n*tk;
    while(Mk<0||Mk>2*pi)
        if(Mk<0)
            Mk=Mk+2*pi;
        else
            Mk=Mk-2*pi;
        end
    end
    
    %  4.计算信号发射时刻的偏近角E
    Ek=Mk;
    while (1)
        E0=Ek;
        Ek = Mk + es * sin(E0);
        if(abs(Ek-E0)<1e-12)
            break;
        end
    end
    
    %  5.计算信号发射时刻的真近点角vk
    cosvk=((cos(Ek)-es)/(1-es*cos(Ek)));
    sinvk=(sqrt(1-es^2))*sin(Ek)/(1-es*cos(Ek));
    vk=atan2(sinvk,cosvk);
    
    %  6.计算信号发射时刻的升交点角距Faik
    Faik=vk+w;
    
    %  7.计算信号发射时刻的摄动校正项Deltuk,Deltrk,Deltik
    Deltuk=Cus*sin(2*Faik)+Cuc*cos(2*Faik);
    Deltrk=Crs*sin(2*Faik)+Crc*cos(2*Faik);
    Deltik=Cis*sin(2*Faik)+Cic*cos(2*Faik);
    
    %  8.计算摄动校正后的升交点角距uk、卫星矢径长度rk、ik
    uk=Faik+Deltuk;
    rk=as*(1-es*cos(Ek))+Deltrk;
    ik=io+dti*tk+Deltik;
    
    %   9.计算信号发射时刻卫星在轨道平面的位置（xk1,yk1）
    xk1=rk*cos(uk);
    yk1=rk*sin(uk);
    
    % 10.计算信号发射时刻的升交点赤经OMGAk
    OMGAk=OMGAo+(dtOMGA-dtOMGAe)*tk-dtOMGAe*toe;
    
    % 11.计算卫星在WGS-84地心地固直角坐标系（Xt,Yt,Zt）中的坐标（xk,yk,zk）
    X=xk1*cos(OMGAk)-yk1*cos(ik)*sin(OMGAk);
    Y=xk1*sin(OMGAk)+yk1*cos(ik)*cos(OMGAk);
    Z=yk1*sin(ik);
    
    % 判断卫星是GEO卫星还是MEO/IGSO卫星
    if num<6
       n1 = 5/180*pi;  %geo旋转角度弧度
            pos = [cos(dtOMGAe*tk)  sin(dtOMGAe*tk)  0;....
                -sin(dtOMGAe*tk) cos(dtOMGAe*tk)  0;
                0   0  1 ]*[1 0 0;0 cos(-n1)  sin(-n1);0 -sin(-n1) cos(-n1)] ...
                *[cos(-dtOMGAe*tk)  sin(-dtOMGAe*tk)  0;....
                -sin(-dtOMGAe*tk) cos(-dtOMGAe*tk)  0;0   0  1 ]*[X;Y;Z];
            X=pos(1);
            Y=pos(2);
            Z=pos(3);
    end
    
    % 从发射时间转换到接收时间坐标系
    dw = dtOMGAe*(basedata.epoch(a1).gps(a2).C2I/299792458);%传播时间转过的角度
    cw = cos(dw);sw = sin(dw);
    anglepos=[cw sw 0;-sw cw 0;0 0 1]*[X;Y;Z];
    
    
    
    
    % 计算高度角 thet
    D = anglepos - x0;
    E = S * D;
    theta=asin(E(3)/sqrt(E(1)^2+E(2)^2+E(3)^2));
    
    if theta>(pi/18)
        basenum = basenum+1;
        
        % 计算FCub(ρu)即在概略点处据卫星的距离与基站与卫星距离差、Ru、卫星角度theta
        Rtheta(basenum) =theta;
        xbs(basenum)=X;        %根据用户接收机计算得到的卫星坐标(xus,yus,zus)
        ybs(basenum)=Y;
        zbs(basenum)=Z;
        BASEprn(basenum) = num;
        rb(basenum)=sqrt((xbs(basenum)-x0(1))^2+(ybs(basenum)-x0(2))^2+(zbs(basenum)-x0(3))^2);
        
        c = 2.99792458e8;%光速
        f1=1561.098e6;
        f2=1207.14e6;
        f3=1268.52e6;
        %% 不同频率下的伪距、载波观测量和进行RTK的计算量
        lamda = c/f1; %波长λ1
        % 伪距、载波观测量
        EWLb.FCb1(basenum)=basedata.epoch(a1).gps(a2).L2I*lamda;
        EWLb.pcb1(basenum) = basedata.epoch(a1).gps(a2).C2I;
        % RTK的计算量
        Base.FCb1(basenum)=basedata.epoch(a1).gps(a2).L2I*lamda-rb(basenum);
        Base.pcb1(basenum) = basedata.epoch(a1).gps(a2).C2I-rb(basenum);
        
        % 伪距、载波观测量
        lamda = c/f2; %波长λ2
        EWLb.FCb2(basenum)=basedata.epoch(a1).gps(a2).L7I*lamda;
        EWLb.pcb2(basenum) = basedata.epoch(a1).gps(a2).C7I;
        % RTK的计算量
        Base.FCb2(basenum)=basedata.epoch(a1).gps(a2).L7I*lamda-rb(basenum);
        Base.pcb2(basenum) = basedata.epoch(a1).gps(a2).C7I-rb(basenum);
        
        % 伪距、载波观测量
        lamda = c/f3; %波长λ3
        EWLb.FCb3(basenum)=basedata.epoch(a1).gps(a2).L6I*lamda;
        EWLb.pcb3(basenum) = basedata.epoch(a1).gps(a2).C6I;
        % RTK的计算量
        Base.FCb3(basenum)=basedata.epoch(a1).gps(a2).L6I*lamda-rb(basenum);
        Base.pcb3(basenum) = basedata.epoch(a1).gps(a2).C6I-rb(basenum);
        
        
    end
    
    
    
end
end




