function [last,next,second,ratio,mean_r,std_r]=timeToLastNextPeak2(stressTime,cExperiment,lineageSource,posesToSort)

%input in the script:
%duration = length of experiment (for frequency calclulations)
%thresholds - for distinguishing between pre-START and post-START cells base on the corresponding 'last' 


%outputs:
%last - time from stress the last reporter peak
%next - duration of the first cell cycle during stress 
%second - duration of the second cell cycle during stress 
%ratio - frequency of divisions before stress / frequency of divisions after 
%mean_r - mean ratio 
%std_r - stamdard deviations of the ratios 
%r/mean/std_pre/post - last three measuremetns specifically for pre-START or post-START cells respectively 

duration=180
threshold1=5;
threshold2=15;
threshold3=20;
threshold4=50;
 
 
if nargin>3
if ischar(posesToSort)
    poses=find(~cellfun(@isempty,strfind(cExperiment.dirs,posesToSort)));
end
else
    poses=1:length(cExperiment.dirs);
end
 
 
birthTimesToSort=ismember(cExperiment.lineageInfo.motherInfo.(lineageSource).motherPosNum,poses);
 
allBirthTimes=cExperiment.lineageInfo.motherInfo.(lineageSource).birthTime(birthTimesToSort,:);

last=zeros(size(allBirthTimes,1),1);
next=zeros(size(allBirthTimes,1),1);
second=zeros(size(allBirthTimes,1),1);
freqA=zeros(size(allBirthTimes,1),1);
freqB=zeros(size(allBirthTimes,1),1);
ratio=zeros(size(allBirthTimes,1),1);
 
last_pre=zeros(size(allBirthTimes,1),1);
last_post=zeros(size(allBirthTimes,1),1);
next_pre=zeros(size(allBirthTimes,1),1);
next_post=zeros(size(allBirthTimes,1),1);
 

postST_A=zeros(size(allBirthTimes,1),1);
postST_B=zeros(size(allBirthTimes,1),1);
 
preST_A=zeros(size(allBirthTimes,1),1);
preST_B=zeros(size(allBirthTimes,1),1);
ratio_pre=zeros(size(allBirthTimes,1),1);
ratio_post=zeros(size(allBirthTimes,1),1);
 
 
for c=1:size(allBirthTimes,1)
 
    birthTimes=allBirthTimes(c,:);
    birthTimes(birthTimes==0)=[];
    birthsBefore=birthTimes(birthTimes<=stressTime);
    lastDivTime=max(birthsBefore);
    
    birthsAfter=birthTimes(birthTimes>stressTime);
    birthsAfterSorted=sort(birthsAfter); % to find 2nd next
    nextDivTime=min(birthsAfter);
    
    freqA(c)=length(birthsAfter)/((duration-stressTime)/12);
    freqB(c)=length(birthsBefore)/(stressTime/12);
    
    if isnan(freqA(c)) || isnan(freqB(c)) || ~any(freqA(c)) || ~any(freqB(c))
        ratio(c)=NaN;
    else 
        ratio(c)=freqB(c)./freqA(c);
        
    end    
    
    if length(birthsAfter)<2
        secondnext=[];
    else    
        secondnext=birthsAfterSorted(2); %2nd next
    end
    
    
    if isempty(lastDivTime)
        last(c)=nan;
    else
        last(c)=stressTime-lastDivTime;
    end
    
    if isempty(nextDivTime)
        next(c)=nan;
    else
        if isempty(lastDivTime)
            next(c)=nan;
        else    
        next(c)=nextDivTime-lastDivTime;
        
        end
    end
    
     if isempty(secondnext)
        second(c)=nan;
    else
        second(c)=secondnext-nextDivTime;
     end
     
end
 
last=last*5; %to conver from time points into minutes
next=next*5; %to conver from time points into minutes
second=second*5; %to conver from time points into minutes
 
for c=1:length(next)
    if next(c)>200 %thresholds for outliers
        next(c)=NaN;
        last(c)=NaN;
        second(c)=NaN;
        freqA(c)=NaN;
        freqB(c)=NaN;
        ratio(c)=NaN;
    end
    
    if last(c)>90 %thresholds for cells of unknown cell cycle stage
        next(c)=NaN;
        last(c)=NaN;
        second(c)=NaN;
        freqA(c)=NaN;
        freqB(c)=NaN;
        ratio(c)=NaN;
    end
end
 
for c=1:size(last)
        if last(c)>=threshold1 && last(c)<=threshold2 
            preST_A(c)=freqA(c);
            preST_B(c)=freqB(c);
            last_pre(c)=last(c);
            next_pre(c)=next(c);
            ratio_pre(c)=preST_B(c)./preST_A(c);
            
        elseif last(c)>=threshold3 && last(c)<=threshold4 
            postST_A(c)=freqA(c);
            postST_B(c)=freqB(c);
            last_post(c)=last(c);
            next_post(c)=next(c);
            ratio_post(c)=postST_B(c)./postST_A(c);
            
        else 
            continue   
        end
        b(c)=last(c)./next(c); 
end
 
ratio(ratio==0)=NaN;
ratio(ratio==Inf)=NaN;
mean_r=nanmean(ratio);
std_r=nanstd(ratio);
 
ratio_pre(ratio_pre==0)=NaN;
ratio_pre(ratio_pre==Inf)=NaN;
mean_rpre=nanmean(ratio_pre);
std_rpre=nanstd(ratio_pre);
 
ratio_post(ratio_post==0)=NaN;
ratio_post(ratio_post==Inf)=NaN;
mean_rpost=nanmean(ratio_post);
std_rpost=nanstd(ratio_post);
 
 
