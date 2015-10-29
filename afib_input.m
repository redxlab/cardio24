function [ op ] = afib_input( input_file, sampFreq )
%AFIB_INPUT_MAY02 Summary of this function goes here
%   Detailed explanation goes here

op=sampFreq;
SampFreq=256;
file=input_file;
[pathstr,fname,ext] = fileparts(file);
%load(file)
Input_ECG=load(file);
[ x1, N, Morphs, Morphs_c,time_markRint, time_RR, P_abs, time_RR2 R_loc,R_value,S_loc,S_value,Q_loc,Q_value,T_loc,T_value,P_loc,P_value ] =afib_check_may02_fn(Input_ECG,sampFreq);
 
t = [0:N-1]/200;   % time index

OPdata=cell(length(R_loc)+1,5);
OPdata{1,1}='P_peak_loc';
OPdata{1,2}='Q_peak_loc';
OPdata{1,3}='R_peak_loc';
OPdata{1,4}='S_peak_loc';
OPdata{1,5}='T_peak_loc';
P_loc_n=[NaN P_loc];
T_loc_n=[T_loc NaN];
for i=2:length(P_loc)+1
    OPdata{i,1}=P_loc_n(i-1);
    OPdata{i,5}=T_loc_n(i-1);
end
for i=2:length(R_loc)+1
    OPdata{i,2}=Q_loc(i-1);
    OPdata{i,3}=R_loc(i-1);
    OPdata{i,4}=S_loc(i-1);
end

OPmData=zeros(length(R_loc)+1,5);
P_loc_n=[0 P_loc];
T_loc_n=[T_loc 0];
for i=2:length(P_loc)
    OPmData(i-1,1)=P_loc_n(i-1);
    OPmData(i-1,5)=T_loc_n(i-1);
end
for i=2:length(R_loc)+1
    OPmData(i-1,2)=Q_loc(i-1);
    OPmData(i-1,3)=R_loc(i-1);
    OPmData(i-1,4)=S_loc(i-1);
end

%EPISODES BEGIN & END POINTS
time_RR2_int=(abs(time_RR2)>0);
%plot(time_RR2_int)
left_epi_RR=find(diff([0 time_RR2_int])==1);
right_epi_RR=find(diff([time_RR2_int 0])==-1);
nEpi=length(left_epi_RR);
Epi= [left_epi_RR; right_epi_RR]';

%P wave Absences
P_abs_int=(abs(P_abs)>0);
P_abs_x=find(diff([0 P_abs_int])==1);
%plot(P_abs_int)

%ALL Parameters into time
% P_loc_n=P_loc_n/sampFreq;
% Q_loc=Q_loc/sampFreq;
% R_loc=R_loc/sampFreq;
% S_loc=S_loc/sampFreq;
% T_loc_n=T_loc_n/sampFreq;
% Epi=Epi/sampFreq;
% P_abs_x=P_abs_x/sampFreq;

%JSON Format saving
Output_JSON_txt_file=strcat(fname,'_op','.txt');
OP_S=struct('P_loc',P_loc_n,'Q_loc',Q_loc,'R_loc',R_loc,'S_loc',S_loc,'T_loc',T_loc_n,'Max_HR',Morphs(1,1),'Avg_HR',Morphs(1,2),'Min_HR',Morphs(1,3),'n_Beats',Morphs(1,4),'n_Irregular_Beats',Morphs(1,5),'Irregular_Beats_Percent',Morphs(1,6),'Avg_PR_Interval',Morphs(1,7),'Avg_QRS_Interval',Morphs(1,8),'Avg_QTc_Interval',Morphs(1,9),'n_Irregular_Beat_Episodes',Morphs(1,10),'n_P_Absences',Morphs(1,11),'n_P_Absence_Episodes',Morphs(1,12),'Irregular_RR_Episodes',Epi,'P_Absences',P_abs_x);
j=savejson('', OP_S,Output_JSON_txt_file);

%Output_file=strcat(fname,'_op','.xls')
%Output_csv_file=strcat(fname,'_op','.csv')
%xlswrite(Output_file,OPdata);
%csvwrite(Output_csv_file,OPmData)

end

