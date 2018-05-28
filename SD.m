function[Ix,Iy,Iz,prn,satnum,maxnum,fw,pw,EWL,FC,pc]=SD(EWLb,basenum,Base,EWLu,obsnum,User)
% 进行单差计算

satnum=0;
for  i=1:obsnum
    for  j=1:basenum
        if(OBSprn(i)==BASEprn(j))
            satnum=satnum+1;
            r=sqrt((xs(i)-x0(1))^2+(ys(i)-x0(2))^2+(zs(i)-x0(3))^2);
            Ix(satnum) = (xs(i) - x0(1))/r;
            Iy(satnum) = (ys(i) - x0(2))/r;
            Iz(satnum) = (zs(i) - x0(3))/r;
            theta1(satnum)= theta(i);
            fw(satnum) = 0.005^2*(1+(1/theta1(satnum)^2));
            pw(satnum) = 0.5^2*(1+(1/theta1(satnum)^2));
            prn(satnum)= OBSprn(i);
            
            % 单差的伪距、载波观测量
            FC(satnum)=User.FCu1(i)-Base.FCb1(j); 
            pc(satnum)=User.pcu1(i)-Base.pcb1(j);
            
%% 单差后的伪距、载波观测量           
            EWL.FC1(satnum)=EWLu.FCu1(i)-EWLb.FCb1(j);
            EWL.pc1(satnum)=EWLu.pcu1(i)-EWLb.pcb1(j);
            
            EWL.FC2(satnum)=EWLu.FCu2(i)-EWLb.FCb2(j);
            EWL.pc2(satnum)=EWLu.pcu2(i)-EWLb.pcb2(j);
            
            EWL.FC3(satnum)=EWLu.FCu3(i)-EWLb.FCb3(j);
            EWL.pc3(satnum)=EWLu.pcu3(i)-EWLb.pcb3(j);
            

            
            
            break;
        end
    end
    
end
thet=0;
for i1=1:satnum
    if(Utheta1(i1)>thet)
        maxnum=i1;
        thet=Utheta1(i1);
    end
end
end