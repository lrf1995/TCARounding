function[obsdata]=readobsfile
% 读取观测文件（用户接收机）
NAVfilepath = '.\cut01680.16o';
fid      = fopen(NAVfilepath);
while ~feof(fid)   %feof若未结束返回0值
    line = fgetl(fid);
    if strfind(line,'END OF HEADER'),break;end %% 跳过文件头
end


obsdata = NaN; 
epochnum = 0;
% 读取接收时间
while ~feof(fid)    
line = fgetl(fid);
epochnum = epochnum +1;
dataline = sscanf(line(2:end),'%f');
timeutc = dataline(1:6);
gpst = cal2gps(timeutc);
gpst(2)=gpst(2)-14; %北斗时
shu = dataline(8);  %该时刻卫星数目
gpsobs = 0;         %导航文件中每个历元中的卫星总个数
%% 根据卫星数目读取信息
for a=1:shu
% 补充空余数据
    line = fgetl(fid);
    linechang=length(line);
     if linechang<193
        line(linechang+1:193)=0;
     end
    
    c2 = str2double(line(5:17));
    l2 = str2double(line(20:33));
    c7 = str2double(line(53:65));
    l7 = str2double(line(68:81));
    c6 = str2double(line(101:113));
    l6 = str2double(line(117:129));

    cada = line(1);
    if (cada=='C')
        if(isnan(c2)||isnan(l2)||isnan(c7)||isnan(l7)||isnan(c6)||isnan(l6))
            continue;
        end
    end
    
   
%% 读取观测数据
    if line(1)=='C'
        gpsobs=gpsobs+1;
        obsdata.epoch(epochnum).gps(gpsobs).bdst =gpst;
        obsdata.epoch(epochnum).gps(gpsobs).gpst =(gpst(1)*604800+gpst(2));% 用户接收机接收机信号接收时间tu
        obsdata.epoch(epochnum).gps(gpsobs).prn = str2double(line(2:3));   % 用户接收的卫星prn号
        obsdata.epoch(epochnum).gps(gpsobs).C2I = str2double(line(5:17));  % 用户接收的卫星伪距量C1L
        obsdata.epoch(epochnum).gps(gpsobs).L2I = str2double(line(20:33)); % 用户接受的卫星载波相位L1C
               
        obsdata.epoch(epochnum).gps(gpsobs).C7I = str2double(line(53:65));
        obsdata.epoch(epochnum).gps(gpsobs).L7I = str2double(line(68:81));
        %obsdata.gps(obsnum).S2W = str2double(line(91:97));
        obsdata.epoch(epochnum).gps(gpsobs).C6I = str2double(line(101:113));
        obsdata.epoch(epochnum).gps(gpsobs).L6I = str2double(line(117:129));
        %obsdata.gps(obsnum).S2X = str2double(line(139:145));
        %obsdata.gps(obsnum).C5Q = str2double(line(149:161));
        %obsdata.gps(obsnum).L5Q = str2double(line(165:177));
        %obsdata.gps(obsnum).S5Q = str2double(line(187:193));
        
        
        %if (obsdata.epoch(epochnum).gps(gpsobs).prn == 13)
        %    gpsobs=gpsobs-1;
        %end
    elseif line(1)=='R'
        
    elseif line(1)=='C'
        
    elseif line(1)=='E'
        
    elseif line(1)=='L'
        
    end
obsdata.epoch(epochnum).gpsobs= gpsobs;
end
end
fclose(fid); 
              
end