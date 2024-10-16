%ExtractTDTDataAvgOnly_Penn

%This program extracts AEP, MUA, and CSD Data from the OpenEx data tank and
%saves the averaged data in individual files

%Note: some of the answers to initial prompt are or can be hardwired into the
%program. Check this before running...

%Note also that in this version of the program, scale factor is set to = 1.0,
%since there is no scaling of the data (data are stored as FP values).

%Important: Before running program, you must change or assign the stim
%codes and filter settings- see below:

clear all;
close all;
clc

Date = '20190416';
TankDir = 'E:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\CSD_Files\MrCassius\190416_CSD\MrCassius-190416';%hardwired answer
%TankName = strcat('Cassius-',Date(3:end));%hardwired answer
TankName = strcat('190416_CSD');
Block = 'CSD_24_24-190416-204908';

BLOCK_PATH = [TankDir,TankName,'\',Block];
BLOCK_PATH ='E:\03_Cohen_Lab\01_Top_Down_Coherence_Project\00_DATA\CSD_Files\MrCassius\190416_CSD\MrCassius-190416\CSD_24_24-190416-204908';

nChannel = 24; %3; % number of channels
EpochStart = -0.1;%Enter the start of epoch in seconds relative to stimulus onset
EpochEnd = 0.4;%Enter the end of epoch in seconds relative to stimulus onset
EpochDuration=EpochEnd-EpochStart;

LimitSwps=150;
ScaleFactor = 1;
ArtRej = 1000;

SampFreq = 12207.03125;
SampPeriod = 1/SampFreq;
time = EpochStart:SampPeriod:EpochEnd; 
time = time' * 1000; % convert time in ms

% ACCESSING TANK FILE
data = TDTbin2mat(BLOCK_PATH,'TYPE',{'epocs','streams'},'CHANNEL',0);
tempdata = TDTfilter(data, 'Tick', 'TIME', [EpochStart, EpochDuration]);
% GET AEP AND MUA OF EACH STIMULUS PRESENTATION    
ReadAEP = tempdata.streams.LFP1.filtered;
ReadMUA = tempdata.streams.MUA1.filtered;
t_AEP = linspace(EpochStart,EpochEnd,size(ReadAEP{1},2));
t_MUA = linspace(EpochStart,EpochEnd,size(ReadMUA{1},2));
t_AEP = t_AEP' * 1000; % in ms
t_MUA = t_MUA' * 1000; % in ms

disp(['Total Number of Sweeps: ',num2str(size(ReadAEP,2))]);
noisy_trial = [];
sweep=1;
while sweep<=size(ReadAEP,2)
%     AEP_trial = ReadAEP{sweep} * -1; % invert the sign of the signal
    AEP_trial = ReadAEP{sweep};
    AEP_trial = AEP_trial / ScaleFactor;
    AEP_trial = AEP_trial * 1000000;
    AEP(:,:,sweep) = AEP_trial;
    
    MUA_trial = abs(ReadMUA{sweep});
    MUA_trial = MUA_trial / ScaleFactor;
    MUA_trial = MUA_trial * 1000000;
    RectMUA(:,:,sweep) = MUA_trial;
    % artifact rejection
    if max(abs(AEP(1,:)))>ArtRej || max(abs(AEP(2,:)))>ArtRej || max(abs(AEP(3,:)))>ArtRej ||...
            max(abs(AEP(4,:)))>ArtRej || max(abs(AEP(5,:)))>ArtRej || max(abs(AEP(6,:)))>ArtRej ||...
            max(abs(AEP(7,:)))>ArtRej || max(abs(AEP(8,:)))>ArtRej || max(abs(AEP(9,:)))>ArtRej ||...
            max(abs(AEP(10,:)))>ArtRej || max(abs(AEP(11,:)))>ArtRej || max(abs(AEP(12,:)))>ArtRej ||...
            max(abs(AEP(13,:)))>ArtRej || max(abs(AEP(14,:)))>ArtRej || max(abs(AEP(15,:)))>ArtRej ||...
            max(abs(AEP(16,:)))>ArtRej || max(abs(AEP(17,:)))>ArtRej || max(abs(AEP(18,:)))>ArtRej ||...
            max(abs(AEP(19,:)))>ArtRej || max(abs(AEP(20,:)))>ArtRej || max(abs(AEP(21,:)))>ArtRej ||...
            max(abs(AEP(22,:)))>ArtRej || max(abs(AEP(23,:)))>ArtRej || max(abs(AEP(24,:)))>ArtRej || sweep>LimitSwps ;%if the absolute value of any channel exceeds ArtRej
        
        noisy_trial = [noisy_trial sweep];
        sweep=sweep+1;
    else 
        sweep=sweep+1;
    end
end
AEP(:,:,noisy_trial) = [];
RectMUA(:,:,noisy_trial) = [];

disp(['Number of Sweeps Accepted: ',num2str(size(AEP,3))]);
disp([' ']);

MeanAEP = transpose(mean(AEP,3));
MeanRectMUA = transpose(mean(RectMUA,3));

%derive CSD:
C1=(-1*MeanAEP(:,1))+(2*MeanAEP(:,2))+(-1*MeanAEP(:,3));
C2=(-1*MeanAEP(:,2))+(2*MeanAEP(:,3))+(-1*MeanAEP(:,4));
C3=(-1*MeanAEP(:,3))+(2*MeanAEP(:,4))+(-1*MeanAEP(:,5));
C4=(-1*MeanAEP(:,4))+(2*MeanAEP(:,5))+(-1*MeanAEP(:,6));
C5=(-1*MeanAEP(:,5))+(2*MeanAEP(:,6))+(-1*MeanAEP(:,7));
C6=(-1*MeanAEP(:,6))+(2*MeanAEP(:,7))+(-1*MeanAEP(:,8));
C7=(-1*MeanAEP(:,7))+(2*MeanAEP(:,8))+(-1*MeanAEP(:,9));
C8=(-1*MeanAEP(:,8))+(2*MeanAEP(:,9))+(-1*MeanAEP(:,10));
C9=(-1*MeanAEP(:,9))+(2*MeanAEP(:,10))+(-1*MeanAEP(:,11));
C10=(-1*MeanAEP(:,10))+(2*MeanAEP(:,11))+(-1*MeanAEP(:,12));
C11=(-1*MeanAEP(:,11))+(2*MeanAEP(:,12))+(-1*MeanAEP(:,13));
C12=(-1*MeanAEP(:,12))+(2*MeanAEP(:,13))+(-1*MeanAEP(:,14));
C13=(-1*MeanAEP(:,13))+(2*MeanAEP(:,14))+(-1*MeanAEP(:,15));
C14=(-1*MeanAEP(:,14))+(2*MeanAEP(:,15))+(-1*MeanAEP(:,16));
C15=(-1*MeanAEP(:,15))+(2*MeanAEP(:,16))+(-1*MeanAEP(:,17));
C16=(-1*MeanAEP(:,16))+(2*MeanAEP(:,17))+(-1*MeanAEP(:,18));
C17=(-1*MeanAEP(:,17))+(2*MeanAEP(:,18))+(-1*MeanAEP(:,19));
C18=(-1*MeanAEP(:,18))+(2*MeanAEP(:,19))+(-1*MeanAEP(:,20));
C19=(-1*MeanAEP(:,19))+(2*MeanAEP(:,20))+(-1*MeanAEP(:,21));
C20=(-1*MeanAEP(:,20))+(2*MeanAEP(:,21))+(-1*MeanAEP(:,22));
C21=(-1*MeanAEP(:,21))+(2*MeanAEP(:,22))+(-1*MeanAEP(:,23));
C22=(-1*MeanAEP(:,22))+(2*MeanAEP(:,23))+(-1*MeanAEP(:,24));

MeanCSD=[C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,C19,C20,C21,C22];

%Now we smooth the data using a 5-point average:

points=5; %specifies number of points to average for smoothing
%MeanAEP:

rows=size(MeanAEP,1);
columns=size(MeanAEP,2);
for c=1:columns
for r=(points+1):rows-(points+1) %start at row #points+1 and ends at row #end-(points+1) for n-point average smooth
MeanAEP(r,c)=mean(MeanAEP(r-points:r+points,c)); %n-point average smooth
end
end
clear c r;

%MeanRectMUA:
rows=size(MeanRectMUA,1);
columns=size(MeanRectMUA,2);
for c=1:columns
for r=(points+1):rows-(points+1) %start at row #points+1 and ends at row #end-(points+1) for n-point average smooth
MeanRectMUA(r,c)=mean(MeanRectMUA(r-points:r+points,c)); %n-point average smooth
end
end
clear c r;

%MeanCSD:
rows=size(MeanCSD,1);
columns=size(MeanCSD,2);
for c=1:columns
for r=(points+1):rows-(points+1) %start at row #points+1 and ends at row #end-(points+1) for n-point average smooth
MeanCSD(r,c)=mean(MeanCSD(r-points:r+points,c)); %n-point average smooth
end
end
clear c r;

%now baseline correct the data to the negative time segment:
p=1;
q=1;
while p==1
    if time(q,1)>0 %finds where the negative time segment terminates
        timezeroindex=q;
       p==0;
       break;
    end
    q=q+1;
end

for r=1:size(MeanAEP,2)
    MeanAEP(:,r)=MeanAEP(:,r)-mean(MeanAEP(1:timezeroindex,r),1);
end

for r=1:size(MeanRectMUA,2)
    MeanRectMUA(:,r)=MeanRectMUA(:,r)-mean(MeanRectMUA(1:timezeroindex,r),1);
end

for r=1:size(MeanCSD,2)
    MeanCSD(:,r)=MeanCSD(:,r)-mean(MeanCSD(1:timezeroindex,r),1);
end


% TimeMeanAEP=[time,MeanAEP];% append time axis
% TimeMeanRectMUA=[time,MeanRectMUA];% append time axis
% TimeMeanCSD=[time,MeanCSD];% append time axis
TimeMeanAEP=[t_AEP,MeanAEP];% append time axis
TimeMeanRectMUA=[t_MUA,MeanRectMUA];% append time axis
TimeMeanCSD=[t_AEP,MeanCSD];% append time axis


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Now We Plot the Data:

figure;
spacingCSD=75;%spacing in microvolts
for n=1:(nChannel-2) % data channels
   z=MeanCSD(:,n);
   points=size(MeanCSD,1);
   grid off
   axis on

%    plot(time,z-n*spacingCSD,'k','LineWidth',1);  %plots CSD
   plot(t_AEP,z-n*spacingCSD,'k','LineWidth',1);  %plots CSD
   hold on
   baseline=(zeros(1,points)-n*spacingCSD);
%    plot(time,baseline,'k','LineWidth',.5); %plots baseline
   plot(t_AEP,baseline,'k','LineWidth',.5); %plots baseline
end
   set(gcf,'Units','inches');%sets units of figure dimensions
   %D=[2 .75 5.75 8.3];%left,bottom,width,height
D=[9 .5 4.75 8.5];%left,bottom,width,height
   set(gcf,'Position',D);%sets figure position to values in D
   axis tight

xlabel('Time [ms]'), ylabel ('Amplitude [microvolts]')
% title(['CSD ',Stimulus,' Inter-Waveform Spacing= ',num2str(spacingCSD)])
title('CSD');
hold off;

% pause(3);
% close

figure;
% spacingRectMUA = 6;%spacing in microvolts
spacingRectMUA = 4;

for n=1:nChannel % data channels
   z=MeanRectMUA(:,n);
   points=size(MeanRectMUA,1);
   grid off
   axis on

%    plot(time,z-n*spacingRectMUA,'k','LineWidth',1);  %plots RectMUA
   plot(t_MUA,z-n*spacingRectMUA,'k','LineWidth',1);  %plots RectMUA
   hold on
   baseline=(zeros(1,points)-n*spacingRectMUA);
%    plot(time,baseline,'k','LineWidth',.5); %plots baseline
   plot(t_MUA,baseline,'k','LineWidth',.5); %plots baseline
end
   set(gcf,'Units','inches');%sets units of figure dimensions
   %D=[2 .75 5.75 8.3];%left,bottom,width,height
   D=[4.5 .5 4.75 8.5];%left,bottom,width,height
   set(gcf,'Position',D);%sets figure position to values in D
   axis tight

xlabel('Time [ms]'), ylabel ('Amplitude [microvolts]')
% title(['RectMUA ',Stimulus,' Inter-Waveform Spacing= ',num2str(spacingRectMUA)])
title('RectMUA');
hold off;

% pause(3);
% close

figure;
spacingAEP=100;%spacing in microvolts

for n=1:nChannel % data channels
   z=MeanAEP(:,n);
   points=size(MeanAEP,1);
   grid off
   axis on

%    plot(time,z-n*spacingAEP,'k','LineWidth',1);  %plots AEP
   plot(t_AEP,z-n*spacingAEP,'k','LineWidth',1);  %plots AEP
   hold on
   baseline=(zeros(1,points)-n*spacingAEP);
%    plot(time,baseline,'k','LineWidth',.5); %plots baseline
   plot(t_AEP,baseline,'k','LineWidth',.5); %plots baseline
end
   set(gcf,'Units','inches');%sets units of figure dimensions
   %D=[2 .75 5.75 8.3];%left,bottom,width,height
   D=[0 .5 4.75 8.5];%left,bottom,width,height
   set(gcf,'Position',D);%sets figure position to values in D
   axis tight

xlabel('Time [ms]'), ylabel ('Amplitude [microvolts]')
% title(['AEP ',Stimulus,' Inter-Waveform Spacing= ',num2str(spacingAEP)])
title('AEP');
hold off;

figure;
subplot(1,2,1);
% pcolor(time,(nChannel-2):-1:1,MeanCSD');
pcolor(t_AEP,(nChannel-2):-1:1,MeanCSD');
colormap jet;
set(gca,'YDir','normal','yTickLabel',(nChannel-2):-2:2);
shading interp;
xlabel('Time [ms]'); ylabel('electrode [ch]');
title('CSD');

subplot(1,3,3);
wn = [0 100]; % clipping window [ms]
% wn = [100 200]; % clipping window [ms] off response
% MeanCSD_clipped = MeanCSD(time>=wn(1)&time<=wn(2),:);
% pcolor(time(time>=wn(1)&time<=wn(2)),(nChannel-2):-1:1,MeanCSD_clipped');
MeanCSD_clipped = MeanCSD(t_AEP>=wn(1)&t_AEP<=wn(2),:);
pcolor(t_AEP(t_AEP>=wn(1)&t_AEP<=wn(2)),(nChannel-2):-1:1,MeanCSD_clipped');
colormap jet;
set(gca,'YDir','normal','yTickLabel',(nChannel-2):-2:2);
shading interp;
title('CSD clipped')
xlabel('Time [ms]'); ylabel('electrode [ch]');

% pause(3);
% close
pause(.1);

clear C1 C2 C3 C4 C5 C6 C7 C8 C9 C10 C11 C12 C13 C14 C15 C16 C17 C18 C19 C20 C21 C22;
clear Channel n MeanAEP MeanRectMUA MeanCSD ReadAEP ReadMUA AEP MUA RectMUA tranges

% end


