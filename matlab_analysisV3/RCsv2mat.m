function [Re,dX,z_res,pixs] = RCsv2mat(R_file,frames)
%Get csv info
%R_file=strtrim(R_file); -Done in the main code-
R_file=[extractBefore(R_file,length(R_file)-3) '.csv'];
Re=readmatrix(R_file);
Re=Re(2:end,:);
x_res=Re(1,2);
y_res=Re(1,3);
Re=Re(:,4:5);

%Calculate values
dX=Re(:,2)/255;
dXmin=min(dX); %subtract the lowest value for background compensation.
dX=dX-dXmin;
%dX=abs(1000*dX); -For MuscleMotion-
dX(isnan(dX))=0;
z_res(:)=size(Re,1);
pixs=x_res*y_res;
Re=round(sum(Re(:,1)/100,'omitnan')*pixs);

%Add zeros to fill matrix
dXz=zeros(frames,1);
dXz(1:z_res)=dX;
dX=dXz;
