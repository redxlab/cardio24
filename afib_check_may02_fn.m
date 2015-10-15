function[x1, N, Morphs,Morphs_c, time_markRint, time_RR, P_abs, time_RR2, R_loc,R_value,S_loc,S_value,Q_loc,Q_value,T_loc,T_value,P_loc,P_value ] = afib_check_apr14_fn( Input_ECG, Fs )

%load('D:\matlab2011a1\matlab2011a\bin\Cardio24\MIT BIH Database Afib hour\04015m1h.mat')
%val=Input_ECG;
  val=Input_ECG';
    val=(val-1024)/200;
    flag=0;
    if(flag==0)
    x1=val(1,1:end);
 else
    x1=-val(1,1:end);
     end 

%  load('D:\matlab2011a1\matlab2011a\bin\Cardio24\r86.mat')
%   x1=r86;
%   x1=x1';

 %x1 = load('ECG3.dat');  % load the ECG signal from the file

fs = Fs;              % Sampling rate
N = length (x1);       % Signal length
t = [0:N-1]/fs;        % time index


% figure(1)
% subplot(2,1,1)
% plot(t,x1)
% xlabel('second');ylabel('Volts');title('Input ECG Signal')
% subplot(2,1,2)
% plot(t(200:1200),x1(200:1200))
% xlabel('second');ylabel('Volts');title('Input ECG Signal 1-3 second')
% xlim([1 6])


%xdum for P wave analysis
xdum=x1;
%Cancellation DC drift and normalization

x1 = x1 - mean (x1 );    % cancel DC conponents
mag_x1=max(abs(x1));
x1 = x1/ max( abs(x1 )); % normalize to one

% figure(2)
% subplot(2,1,1)
% plot(t,x1)
% xlabel('second');ylabel('Volts');title(' ECG Signal after cancellation DC drift and normalization')
% subplot(2,1,2)
% plot(t(200:1200),x1(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])

%Low Pass Filtering
% LPF (1-z^-6)^2/(1-z^-1)^2
b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
a=[1 -2 1];

h_LP=filter(b,a,[1 zeros(1,12)]); % transfer function of LPF
x2 = conv (x1 ,h_LP);
%x2 = x2 (6+[1: N]); %cancle delay
x2 = x2/ max( abs(x2 )); % normalize , for convenience .

% figure(3)
% subplot(2,1,1)
% plot([0:length(x2)-1]/fs,x2)
% xlabel('second');ylabel('Volts');title(' ECG Signal after Low Pass Filtering')
% xlim([0 max(t)])
% subplot(2,1,2)
% plot(t(200:1200),x2(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])

%High Pass Filtering
% HPF = Allpass-(Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
a = [1 -1];

h_HP=filter(b,a,[1 zeros(1,32)]); % impulse response iof HPF

x3 = conv (x2 ,h_HP);
%x3 = x3 (16+[1: N]); %cancle delay
x3 = x3/ max( abs(x3 ));

% figure(4)
% subplot(2,1,1)
% plot([0:length(x3)-1]/fs,x3)
% xlabel('second');ylabel('Volts');title(' ECG Signal after High Pass Filtering')
% xlim([0 max(t)])
% subplot(2,1,2)
% plot(t(200:1200),x3(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])


%Derivative Filter

% Make impulse response
h = [-1 -2 0 2 1]/8;
% Apply filter
x4 = conv (x3 ,h);
x4 = x4 (2+[1: N]);
x4 = x4/ max( abs(x4 ));

% figure(5)
% subplot(2,1,1)
% plot([0:length(x4)-1]/fs,x4)
% xlabel('second');ylabel('Volts');title(' ECG Signal after Derivative Filter')
% subplot(2,1,2)
% plot(t(200:1200),x4(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])


%Squaring
x5 = x4 .^2;
x5 = x5/ max( abs(x5 ));

% figure(6)
% subplot(2,1,1)
% plot([0:length(x5)-1]/fs,x5)
% xlabel('second');ylabel('Volts');title(' ECG Signal Squaring')
% subplot(2,1,2)
% plot(t(200:1200),x5(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])


%Moving Window Integration

% Make impulse response
h = ones (1 ,31)/31;
Delay = 15; % Delay in samples

% Apply filter
x6 = conv (x5 ,h);
x6 = x6 (15+[1: N]);
x6 = x6/ max( abs(x6 ));

% figure(7)
% subplot(2,1,1)
% plot([0:length(x6)-1]/fs,x6)
% xlabel('second');ylabel('Volts');title(' ECG Signal after Averaging (Moving Window Integration')
% subplot(2,1,2)
% plot(t(200:1200),x6(200:1200))
% xlabel('second');ylabel('Volts');title(' ECG Signal 1-3 second')
% xlim([1 6])

%Find QRS Points Which it is different than Pan-Tompkins algorithm

% figure(7)
% subplot(2,1,1)
max_h = max(x6);
thresh = mean (x6 );
poss_reg =(x6>thresh*max_h);

% figure (8)
% subplot(2,1,1)
% hold on
% plot (t(200:1200),x1(200:1200)/max(x1))
% box on
% xlabel('second');ylabel('Integrated')
% xlim([1 6])
% 
% subplot(2,1,2)
% plot (t(200:1200),x6(200:1200)/max(x6))
% xlabel('second');ylabel('Integrated')
% xlim([1 6])

left = find(diff([0 poss_reg])==1);
right = find(diff([poss_reg 0])==-1);

left=left-(6+16);  % cancle delay because of LP and HP
right=right-(6+16);% cancle delay because of LP and HP

R_value=zeros(1,length(left));
S_value=zeros(1,length(left));
Q_value=zeros(1,length(left));
for i=1:length(left)
    if(left(i)<1)
       left(i)=1;
    end
    if(right(i)<1)
       right(i)=1;
    end        
        [R_value(i) R_loc(i)] = max( x1(left(i):right(i)) );
        R_loc(i) = R_loc(i)-1+left(i); % add offset

        [Q_value(i) Q_loc(i)] = min( x1(left(i):R_loc(i)) );
        Q_loc(i) = Q_loc(i)-1+left(i); % add offset

        [S_value(i) S_loc(i)] = min( x1(R_loc(i):right(i)) );
        S_loc(i) = S_loc(i)-1+R_loc(i); % add offset

end

% there is no selective wave
% Q_loc=Q_loc(find(Q_loc~=0));
% R_loc=R_loc(find(R_loc~=0));
% S_loc=S_loc(find(S_loc~=0));


%Finding T wave %similar to previous algo 

% T_loc=ones(1,length(left)-1);
% for z=1:length(left)-1
%     len=floor(0.75*(left(z+1)-right(z)));
%     xt1=x1(right(z):right(z)+len);
%     %Derivative Filter
%     % Make impulse response
%     h = [-1 -2 0 2 1]/8;
%     xt2 = conv (xt1 ,h);
%     xt2 = xt2 (2+[1: length(xt1)]);
%     xt2 = xt2/ max( abs(xt2 ));
%     %squaring
%     xt3 = xt2 .^2;
%     xt3 = xt3/ max( abs(xt3 ));
%     %average filter
%     h = ones (1 ,31)/31;
%     Delay = 30; 
%     xt4 = conv (xt3 ,h);
%     xt4 = xt4 (15+[1: length(xt1)]);
%     xt4 = xt4/ max( abs(xt4 ));
%     %threshold
%     max_t = max(xt4);
%     thresh_t = mean (xt4 );
%     poss_reg_t =(xt4>thresh_t*max_t);
%     left_t = find(diff([0 poss_reg_t])==1);
%     right_t = find(diff([poss_reg_t 0])==-1);
%     xt5=zeros(1,length(xt4));
%     for k=1:length(left_t)    
%         xt5(left_t(k):right_t(k))=1;    
%     end
%      figure(10+z)
%        subplot(3,1,1)
%        plot(xt1)
%        subplot(3,1,2)
%        plot(xt4)
%        subplot(3,1,3)
%        plot(xt5)
%     if(length(left_t)==0)
%         disp('');
%     elseif(left_t(1)==1&&length(left_t)>1)
%          xt1_diff=diff(xt1(left_t(2):right_t(2)));
%          for q=left_t(2):right_t(2)-2
%              if (xt1_diff(q-left_t(2)+1)<0.001)
%                  T_loc(z)=q-1+left_t(2)+right(z);
%                  %T_loc(z)=q-1+right(z);
%                  break
%              end             
%          end         
% %          [T_value(z) T_loc(z)] = (max( xt1(left_t(2):right_t(2))));                
% %          T_loc(z)=T_loc(z)-1+left_t(2)+right(z);
%          T_value(z)=x1(T_loc(z));
%     else         
%          xt1_diff=diff(xt1(left_t(1):right_t(1)));
%          for q=left_t(1):right_t(1)-2
%              if (abs(xt1_diff(q-left_t(1)+1))<0.001)
%                  T_loc(z)=q-1+left_t(1)+right(z);
%                  %T_loc(z)=q-1+right(z);
%                  break
%              end             
%          end
% %         [T_value(z) T_loc(z)] = (max( xt1(left_t(1):right_t(1))));        
% %         T_loc(z)=T_loc(z)-1+left_t(1)+right(z);
%         T_value(z)=x1(T_loc(z));           
%     end    
%       
% end
%%
% %%
% %Finding T waves - algo similar to P waves
% T_loc=ones(1,length(left)-1);
% for z=2:length(left)-1
%     len=floor(0.6*(Q_loc(z+1)-right(z)));
%     xt1=x2(right(z):right(z)+len);
%     [T_value(z) T_loc(z)] = (max( xt1));
%     T_loc(z)=T_loc(z)-1+right(z);
%     T_value(z)=x1(T_loc(z));    
%     
%     figure(10+z)
%       subplot(3,1,1)
%         plot(x1(right(z):right(z)+len))
%         subplot(3,1,2)
%         plot(xt1)
%         subplot(3,1,3)
%         plot(xt1)
%end


%%
%Updated T wave algorithm
R_value(5);
S_value(5);
uplimit=(sum(R_value(1,:)))/length(R_value); % avg value of maxima
lowlimit=(sum(S_value(1,:)))/length(S_value); %avg value of minima
T_thresh=0.15*abs(uplimit-lowlimit); %threshold for detection of T

T_loc=ones(1,length(left)-1);
T_limit_left=zeros(1,length(left)-1);
T_limit_right=zeros(1,length(left)-1);
T_loc_left=ones(1,length(left)-1);
T_loc_right=ones(1,length(left)-1);
T_further_right=ones(1,length(left)-1);
T_further_right_value=ones(1,length(left)-1);
T_inv=zeros(1,length(left)-1);

for z=2:length(left)-1
    T_inv(z)=0; %T wave assumed to be not inverted
    len=floor(0.85*(Q_loc(z+1)-right(z)));
    xt1=x2(right(z):right(z)+len);
    [T_value(z) T_loc(z)] = (max( xt1));
    T_loc(z)=T_loc(z)-1+right(z);
    T_value(z)=x1(T_loc(z)); 
    
    %making further test to see if really T is upright
    
    right_max_T=T_loc(z)+floor(300/1000*fs); %taking a point right of maxima
    diff_T=0;
    if(right_max_T>length(x1)) %to avoid segmentation fault
        right_max_T=length(x1)-1;
    end
        diff_T=abs(x1(right_max_T)-x1(T_loc(z)));
    
    %if region to right of max is more or less flat then it is not a peak
    if(diff_T<0.03)
        T_inv(z)=1;
        %disp('cond1');
    end
    if(T_inv(z)==0)
        T_limit_left(z)=floor(T_loc(z)-(200/1000)*fs);  %locations of left of max
        T_limit_right(z)=floor(T_loc(z)+(200/1000)*fs);  %locations of right of max
        
        if(T_limit_left(z)>length(x1)) %to avoid segmentation fault
        T_limit_left(z)=length(x1)-1;
        end
        
        if(T_limit_right(z)>length(x1)) %to avoid segmentation fault
        T_limit_right(z)=length(x1)-1;
        end
        
        [T_value_left(z) T_loc_left(z)]=min(x1(T_limit_left(z):T_loc(z))); %finding left ending of T if exists
        [T_value_right(z) T_loc_right(z)]=min(x1(T_loc(z):T_limit_right(z)));%finding right ending of T if exists
        
        T_loc_left(z)=T_loc_left(z)+T_limit_left(z);  %adding offset
        T_loc_right(z)=T_loc_right(z)+T_loc(z); %adding offset
        
        T_value_left(z)=x1(T_loc_left(z));
        T_value_right(z)=x1(T_loc_right(z));
        
        T_further_right(z)=T_loc_right(z)+floor(100/1000*fs); %any point to right of right cross of T
       
        if(T_further_right(z)>length(x1)) %to avoid segmentation fault
        T_further_right(z)=length(x1)-1;
        end
        
        dif_T=abs(x1(T_further_right(z))-x1(T_loc_right(z))); %checking to see if it is peak or left part of inverted T wave
        T_further_right_value(z)=x1(T_further_right(z));
       
        if(dif_T>T_thresh) % means that the wave rises again indicating it was a part of inverted T
            T_inv(z)=1; %T is inverted
            %disp('cond2');
        end
        
%         %adjusting to avoid negative or zero subscript
%         if(T_value_left(z)<=0)
%             T_value_left(z)=1;
%         end
%         
%         if(T_loc_left(z)<=0)
%             T_loc_left(z)=1;
%         end
%         
%          if(T_value_right(z)<=0)
%             T_value_right(z)=1;
%         end
%         
%         if(T_loc_right(z)<=0)
%             T_loc_right(z)=1;
%         end
        
        
    end
    if (T_inv(z)==0)
      %disp('1');
    end
    if(T_inv(z)==1)
        [T_value(z) T_loc(z)] = (min( xt1));
        T_loc(z)=T_loc(z)-1+right(z);
        T_value(z)=x1(T_loc(z)); 
        %disp('-1');
    end
%     disp('T Value')
end

%%
%%
%%finding P Waves 
P_loc=ones(1,length(left)-1);
P_loc_left=ones(1,length(left)-1);
P_loc_right=ones(1,length(left)-1);
P_limit_left=zeros(1,length(left)-1);
P_limit_right=zeros(1,length(left)-1);
P_value_left=ones(1,length(left)-1);
P_value_right=ones(1,length(left)-1);
baseline_amp=zeros(1,length(left)-1);
amp_p=zeros(1,length(left)-1);
width_p=zeros(1,length(left)-1);
count_p=ones(1,length(left)-1);
T_inv_rt_val=zeros(1,length(left)-1);
T_inv_rt_loc=zeros(1,length(left)-1);

for z=1:length(left)-1
    len1=floor(0.3*(Q_loc(z+1)-right(z)));% 30 % of the region taken for P estimation
    len2=Q_loc(z+1)-T_loc(z); %taking length from T to next Q
    len=min(len1,len2);%taking whichever is minimum for finding P 
    xp1=x2(Q_loc(z+1)-len:Q_loc(z+1));
%     figure(20+z)
%     plot(xp1)
    if(numel(xp1)>5) 
        [P_value(z) P_loc(z)] = (max( xp1));
        P_loc(z)=P_loc(z)-1+Q_loc(z+1)-len-6;
        P_value(z)=x1(P_loc(z));
        
        P_limit_left(z)=floor(P_loc(z)-(55/1000)*fs);  %locations of left of P
        P_limit_right(z)=floor(P_loc(z)+(55/1000)*fs);  %locations of right of P
        
        
        [P_value_left(z) P_loc_left(z)]=min(x1(P_limit_left(z):P_loc(z)));
        [P_value_right(z) P_loc_right(z)]=min(x1(P_loc(z):P_limit_right(z)));
        
        P_loc_left(z)=P_loc_left(z)+P_limit_left(z);%adding offset
        P_loc_right(z)=P_loc_right(z)+P_loc(z);%adding offset
        
        [P_value(z) P_loc(z)] = max(x1(P_loc_left(z):P_loc_right(z)));
        P_loc(z)=P_loc(z)+P_loc_left(z)-1;
        P_value(z)=x1(P_loc(z));
        
        baseline_amp(z)=sum(x1(P_loc_left(z))+x1(P_loc_right(z)))/2; %calc. baseline of p
        amp_p(z)=abs(P_value(z)-baseline_amp(z)); %amplitude of p waves
        
       width_p(z)=abs(P_loc_right(z)-P_loc_left(z));%width of p waves
       
       if(width_p(z)<(60/1000*fs))  %||amp_p(z)<0.008||amp_p(z)>0.050
          % disp('BAD P');
           count_p(z)=0;
       else %disp('GOOD P');
       end
       
       if(T_inv(z)==1)
            xT_right=x1(T_loc(z):floor(T_loc(z)+(200/1000*fs))); %finding right limit of inv T
            [T_inv_rt_val(z) T_inv_rt_loc(z)]=max(xT_right);  %finding right edge of inv T wave
            T_inv_rt_loc(z)=T_inv_rt_loc(z)+T_loc(z)-1; %adding offset
            T_inv_rt_val(z)=x1(T_inv_rt_loc(z)); %updating value
            if(P_loc_left(z)-T_inv_rt_loc(z)<0.01*fs) % if left of P wave is very close to right of inv T wave
                count_p(z)=0;
            end
       end
       
        %||amp_p(z)<0.010||amp_p(z)>0.025
       
%        slope_thresh=1e-5;
%        if(T_inv(z)==1)
%             xT_right=x1(T_loc(z):floor(T_loc(z)+(200/1000*fs))); %finding right limit of inv T
%             [T_inv_rt_val(z) T_inv_rt_loc(z)]=max(xT_right);  %finding right edge of inv T wave
%             T_inv_rt_loc(z)=T_inv_rt_loc(z)+T_loc(z)-1; %adding offset
%             T_inv_rt_val(z)=x1(T_inv_rt_loc(z)); %updating value
%             slope_rt=abs(0.1*(T_inv_rt_val(z)-x1(T_loc(z)))/(T_inv_rt_loc(z)-T_loc(z))) %slope of right arm of inverted T
%             slope_T_p=abs(0.1*(x1(P_loc_left(z))-x1(T_loc(z)))/(P_loc_left(z)-T_loc(z))); %slope of left of P with T
%             delta_slope=abs(slope_rt-slope_T_p); %diff in slope
%             if(delta_slope<slope_thresh) %if p lies in the right arm of inv T then this will be very small
%                 count_p(z)=0;  %in that case reject P as proper
%             end
%        end       
      
        if (P_loc_left(z)<0)
            P_loc_left=1;
        end
        
%         [pks,locs] = findpeaks(xp1);
%         [sortpks, sortlocs]=sort(pks, 'descend');  
%         if(numel(sortlocs)~=0)
%             P_loc(z+1)=locs(sortlocs(1))-1+Q_loc(z+1)-len-6;
%             P_value(z+1)=x1(P_loc(z+1));
%         end
    end    
end

Avg_PR=0;
sum_PR=0;
for i=1:length(count_p)
if(count_p(i)==1)
    for j=1:length(R_loc)
        if(R_loc(j)>P_loc_left(i))
            break;
        end
    end  
      sum_PR=sum_PR+(R_loc(j)-P_loc(i));
      Avg_PR=(sum_PR/length(count_p)/fs)*1000;
end
end

P_pres=zeros(1,length(x1));
P_abs=zeros(1,length(x1));
P_pres_inv=0.5*ones(1,length(x1));
for i=1:length(count_p)
    if(count_p(i)==1)
        P_pres(P_loc(i))=count_p(i)*0.2;
    end
    if(count_p(i)==0)
        P_abs(P_loc(i))=0.3; 
    end
   P_pres_inv(i)=0;
end

num_p=sum(count_p);
plot_good_p=zeros(1,length(left)-1);
for i=1:length(left)-1
    if(count_p(i)==1)
        plot_good_p(i)=P_value(i);
       
    else
            plot_good_p(i)=NaN;
    end

end

% figure(9)
% subplot(2,1,1)
% title('ECG Signal with R points');
% plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o', t(T_loc),T_value, 'go',t(P_loc_ref), plot_good_p, 'k*',t(P_loc_left),P_value_left,'rX',t(P_loc_right),P_value_right,'rX');
% legend('ECG','R','S','Q');
%subplot(2,1,2)
%plot (t,x1/max(x1) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o', t(T_loc),T_value, 'go',t(P_loc), P_value, 'k*');
%xlim([2 8])

len_morph=length(x1);
QT=zeros(1,len_morph);
QRS=zeros(1,len_morph);
for h=1:length(R_loc)-1
   if(T_loc(h)~=1)
      QT(R_loc(h))=(T_loc(h)-Q_loc(h))/fs*1000;
   end
      QRS(R_loc(h))=(right(h)-Q_loc(h))/fs*1000;
end
Avg_QRS_Interval=sum(QRS)/length(R_loc);
Avg_QT_Interval=sum(QT)/length(R_loc);
%R loc stored in R_loc
%computing successive differences in the whole data
symp1=0; %total no of irr beats>5% of tot beat
symp2=zeros(1,length(R_loc)-20); %>3 consecutive beats outside interval mean(RR)+-10% --no of afibs in each window
symp3=zeros(1,length(R_loc)-20); %<95% of beats in the interval mean(RR)+-10%
Afib_dur=0;
diffdata=diff(R_loc(1:end));
time_RR=zeros(1,length(x1));
for i=1:length(diffdata)
    time_RR(R_loc(i))=(diffdata(i)/fs)*1000;    
end

avgR=((sum(diffdata)/length(R_loc))/fs)*1000;
markRint=zeros(1,length(R_loc));
time_markRint=zeros(1,length(x1));
Avg_HR = 60000/avgR;
Total_QRS=length(diffdata);
Avg_RR_Interval=sum(diffdata)/length(diffdata)/fs*1000;
Avg_QTc=Avg_QT_Interval/sqrt(Avg_RR_Interval/1000);
for i=2:length(R_loc)-1
    if((diffdata(i)-diffdata(i-1))/diffdata(i)>.10)
       markRint(i)=1; 
       time_markRint(R_loc(i))=0.9;
    end
end
num_irr=sum(markRint); %no of irr RR intervals
% figure(10)
% plot(t,time_markRint,'b')
% hold on
% plot(t,x1/max(x1),'g')
% hold on 
% plot(t,time_RR/max(time_RR),'r')

%condition 1 of Dr.MM's criteria for afib
%disp(' No. of isolated irregular beats: ')
%disp(num_irr)
Num_IRR=num_irr;
Perc_IRR=num_irr/length(diffdata)*100;
if(num_irr>=(0.05*length(R_loc)));
    symp1=1; %symp1=1 if condition 1 is satisfied
end
xf=zeros(1,length(x1));
time_RR2=zeros(1,length(x1));
Max_HR=0;
Min_HR=300;
Num_Consec_Beats=0;
for i=1:length(R_loc)-20
    window=(x1(R_loc(i):R_loc(i+20))); 
    time_win=[0:length(window)-1]/fs;
    beat_dur=((diff(R_loc(i:i+20))/fs)*1000); 
    avg_beat_dur=(sum(beat_dur)/20);  %this is in sec and fs=200 hopefully
    bpm(i)=60000/avg_beat_dur;
    if(bpm(i)>Max_HR && bpm(i)<=180)
        Max_HR=bpm(i);
    end
    if(bpm(i)<Min_HR && bpm(i)>=40)
        Min_HR=bpm(i);
    end
    %risk_bpm=0;
    start=0;
    count=0;
    thresh=avg_beat_dur;
    bool_window=zeros(1,20);    
    % checking for RR outside threshold region in the window
    for j=1:20
        if((beat_dur(j))<(thresh-.1*thresh)||(beat_dur(j))>(thresh+0.1*thresh))
            bool_window(j)=1;                    
        end          
    end 
    
%     figure(20+i)
%     subplot(3,1,1)
%     plot(time_win,window)
%     subplot(3,1,2)
%     plot(bool_window)
%     subplot(3,1,3)    
%     plot(beat_dur)    
    
    sums=0;%bool_window(1); %to store 1's for irregularities
    summax=sums;%max sum upto a certain iteration
    for j=1:19
        if(bool_window(j+1)>=bool_window(j))
            sums=sums+bool_window(j+1);
        end
        if(bool_window(j+1)<bool_window(j))
            if(sums>summax)
                summax=sums;
            end
            sums=0;
        end
    end
    len_afib=summax; %length of max sub-series of 1's (irr)
    symp2(i)=len_afib;    
    %testing if >3 consecutive irr occur in a window and thus 
    %marking the window 
    qre=0;
     if(len_afib>10)           
            if(xf(R_loc(i)+1)==1)
            %do nothing
            else
                xf(R_loc(i):R_loc(i+20))=1;
                time_RR2(R_loc(i):R_loc(i+20))=x1(R_loc(i));
                %P_pres(R_loc(i):R_loc(i+20))
                Afib_dur=Afib_dur+(R_loc(i+20)-R_loc(i))/fs;
                qre=qre+1;  
                Num_Consec_Beats=Num_Consec_Beats+1;
            end
     else            
            if(xf(R_loc(i)+1)==1)
            %do nothing
            else
                xf(R_loc(i):R_loc(i+20))=.5;
            end
     end       
     %condition 3 Dr.MM ECG diagnostic criterias
     tot_irr=sum(bool_window);
     if(tot_irr/20>.05)
         symp3(i)=1;  %increasing symp3 if condition 3 satisfied
     end     
end

%%
%Afib P and IRR
count_p_inv=ones(1,length(count_p));
for i=1:length(count_p)
    if(count_p(i)==0)
        count_p_inv(i)=1;
    elseif(count_p(i)==1)
        count_p_inv(i)=0;
    end
end
Num_P_abs=sum(count_p_inv);
Num_consec_P_abs=0;
for j=1:length(count_p)-4
   if(count_p_inv(j)==1&&count_p_inv(j+1)==1&&count_p_inv(j+2)==1&& count_p_inv(j+3)==1)
        Num_consec_P_abs=Num_consec_P_abs+1;
   end
end

% cpi_x_time_RR2=count_p_inv.*time_RR2;



%figure(300)
%plot (t,x1/max(x1), t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o', t(T_loc),T_value, 'go',t(P_loc), plot_good_p, 'k*',t(P_loc_left),P_value_left,'rX',t(P_loc_right),P_value_right,'rX');
% plot (t,x1/max(x1),'g',t(P_loc_ref), plot_good_p, 'k*');
% hold on 
% plot(t,time_RR/max(time_RR),'r')
% hold on
% plot(t,time_RR2,'b')

% figure(301)
% %plot (t,x1/max(x1), t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o', t(T_loc),T_value, 'go',t(P_loc), plot_good_p, 'k*',t(P_loc_left),P_value_left,'rX',t(P_loc_right),P_value_right,'rX');
% plot (t,x1/max(x1),'g',t(P_loc_ref), plot_good_p, 'k*');
% hold on 
% plot(t,time_RR/max(time_RR),'r')
% hold on
% plot(t,time_RR2,'b')

%P Presence
%figure(302)
%plot (t,x1/max(x1), t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o', t(T_loc),T_value, 'go',t(P_loc), plot_good_p, 'k*',t(P_loc_left),P_value_left,'rX',t(P_loc_right),P_value_right,'rX');
% plot (t,x1/max(x1),'g');
% hold on 
% plot(t,time_RR/max(time_RR),'r')
% hold on
% plot(t,P_abs,'k')
% hold on
% plot(t,time_RR2,'b')

Morphs=zeros(1,12);

for i=1:12
    Morphs_c(i,:)='  Abnormal';
end
%disp('Report')
%disp(' Max HR:')
%disp(Max_HR)
Morphs(1,1)=Max_HR;
if (Max_HR<100)
    Morphs_c(1,:)='  Normal  ';
end
%disp('Avg HR:')
%disp(Avg_HR);
Morphs(1,2)=Avg_HR;
if (60<Avg_HR<100)
    Morphs_c(2,:)='  Normal  ';
end
%disp('Min_HR:')
%disp(Min_HR);
Morphs(1,3)=Min_HR;
if (Min_HR>60)
    Morphs_c(1,:)='  Normal  ';
end
%disp('Total Number of QRS: ')
%disp(Total_QRS);
Morphs(1,4)=Total_QRS;
Morphs_c(4,:)='  ------  ';
%disp('Number of irregual beats:');
%disp(Num_IRR)
Morphs(1,5)=Num_IRR;
%disp('Percentage of irregular beats:')
%disp(Perc_IRR);
Morphs(1,6)=Perc_IRR;
if (Perc_IRR<10)
    Morphs_c(5,:)='  Normal  ';
    Morphs_c(6,:)='  Normal  ';
end
%disp('Number of episodes with consec beats:')
%disp(Num_Consec_Beats);
Morphs(1,7)=Num_Consec_Beats;
if (Num_Consec_Beats<25)
    Morphs_c(7,:)='  Normal  ';
end
%disp('Intervals')
%disp('PR:')
%disp(Avg_PR)
Morphs(1,8)=Avg_PR;
if (120<Avg_PR<200)
    Morphs_c(8,:)='  Normal  ';
end
%disp('QRS:')
%disp(Avg_QRS_Interval)
Morphs(1,9)=Avg_QRS_Interval;
if (60<Avg_PR<100)
    Morphs_c(9,:)='  Normal  ';
end
%disp('Avg_QTc:')
%disp(Avg_QTc);
Morphs(1,10)=Avg_QTc;
if (Avg_QTc<430)
    Morphs_c(10,:)='  Normal  ';
end
%disp('No of p absences');
%disp(Num_P_abs);
Morphs(1,11)=Num_P_abs;
if ((Num_P_abs*100)/Total_QRS<8)
    Morphs_c(11,:)='  Normal  ';
end
%disp('Num of consec > 4 P absences')
%disp(Num_consec_P_abs)
Morphs(1,12)=Num_consec_P_abs;
if ((Num_consec_P_abs*100)/Total_QRS<1)
    Morphs_c(12,:)='  Normal  ';
end

%%
% count_p_inverted=ones(1,length(count_p));
% loc=zeros(1,length(count_p));
% for i=1:length(count_p)
%    if(count_p(i)==1)
%        count_p_inverted(i)=0;
%    end
% end
% for j=1:length(count_p)
%    if(count_p_inverted(i)==1&&count_p_inverted(i+1)==1&&count_p_inverted(i+2)==1&& count_p_inverted(i+3)==1)
%         loc(j)=j;
%    end
% end
end

