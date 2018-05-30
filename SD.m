function[Ix,Iy,Iz,prn,satnum,maxsat,fw,pw,EWL,FC,pc]=SD(EWLb,basenum,Base,EWLu,obsnum,User,x0)
% 进行单差计算

satnum=0;
for  i=1:obsnum
    for  j=1:basenum
        if(User(i).prn == Base(j).prn)
            satnum=satnum+1;
            r=sqrt((User(i).xs-x0(1))^2+(User(i).ys-x0(2))^2+(User(i).zs-x0(3))^2);
            Ix(satnum) = (User(i).xs - x0(1))/r;
            Iy(satnum) = (User(i).ys - x0(2))/r;
            Iz(satnum) = (User(i).zs - x0(3))/r;
            theta1(satnum)= User(i).theta;
            fw(satnum) = 0.005^2*(1+(1/theta1(satnum)^2));
            pw(satnum) = 0.5^2*(1+(1/theta1(satnum)^2));
            prn(satnum) = User(i).prn;
            
            % 单差的伪距、载波观测量
            FC(satnum)=User(i).FC1-Base(j).FC1; 
            pc(satnum)=User(i).pc1-Base(j).pc1;
            
%% 单差后的伪距、载波观测量           
            EWL.FC1(satnum)=EWLu(i).FC1-EWLb(j).FC1;
            EWL.pc1(satnum)=EWLu(i).pc1-EWLb(j).pc1;
            
            EWL.FC2(satnum)=EWLu(i).FC2-EWLb(j).FC2;
            EWL.pc2(satnum)=EWLu(i).pc2-EWLb(j).pc2;
            
            EWL.FC3(satnum)=EWLu(i).FC3-EWLb(j).FC3;
            EWL.pc3(satnum)=EWLu(i).pc3-EWLb(j).pc3;
            

            
            
            break;
        end
    end
    
end
maxthet=0;
for i1=1:satnum
    if(theta1(i1)>maxthet)
        maxsat=i1;
        maxthet=theta1(i1);
    end
end
end