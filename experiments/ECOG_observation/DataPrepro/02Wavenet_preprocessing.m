%note that the motion should corresponds to the data set
%0628S2
%The monkeydata
fs=1000;
Fc=centfrq('morl');

%Fa=Fc*fs/a;

Fa=(2:11)*10;
SCALES=Fc*fs./Fa;

%S=1:100;

%COEFS = cwt(S,SCALES,'morl');

F = scal2frq(SCALES,'morl',1/fs);

%imagesc(S,F,abs(COEFS))

chall = load("chall_big_monkey1.csv")

%frequence filter.
sampleRate = 1000;
highEnd = 499; 
lowEnd = 0.3; 
filterOrder = 3; 
[b, a] = butter(filterOrder, [lowEnd highEnd]/(sampleRate/2)); 
for i=1:64
  i  
chall(i,:)= filtfilt(b, a, chall(i,:));
end


%reference
[~,n]=size(chall(1,:));

for j=1:n
    refer=-sum(chall(:,j))/64;
    for i=1:64
        chall(i,j)=chall(i,j)+refer;
    end
end

        
    




motion=load('Motion_1.mat');


%wavelet transformation
timeindex=15001:25000;

timeused=floor(1000*motion.MotionTime(timeindex));

%size(timeused)
chall_time_frenquen=zeros(64*10000,100);

tic;
for i=1:64
    for j=1:10000
    S=chall(i,floor((timeused(j)-1000+1):timeused(j)));
    COEFS = cwt(S,SCALES,'morl');
    COEFS_before=abs(COEFS(:,100:100:1000));
    COEFS_middle=zscore(COEFS_before');
    COEFS_middle=COEFS_middle';
    COEFS_middle_two=COEFS_middle;
    %COEFS_middle_two=COEFS_middle(:,100:100:1000);
    chall_time_frenquen((j-1)*64+i,:)=reshape(COEFS_middle_two,1,100);
 
end
end
toc;

    
csvwrite('chall_time_frenquen_1.csv',chall_time_frenquen)

%y=motion.MotionData{1}(timeindex,1)-motion.MotionData{4}(timeindex,1);

y=motion.MotionData{1}(timeindex,1);
csvwrite('y_1.csv',y)


