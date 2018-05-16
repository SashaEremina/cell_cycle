function resultsStruct=measureFreq(cExperiment,poses,stressTime)

%measures the distance between the peaks and calculates it as a temporal
%fraction of the cell cycle taking the adjacent datapoints in the lineage2
%as the begining and the end of each cell cycle
 
%corresponds to one sliding window search method in Sasha's report
%(increments in bud lineage only)
 
%sources have to be defined in accordance with positions. if there is any
%discrepancy between the two - sources will be used as a reference

sources={'htb2Edited'};

if nargin>1
    if ischar(poses)
        poses=find(~cellfun(@isempty,strfind(cExperiment.dirs,poses)));        
    end
    cellsToPlot=find(ismember(cExperiment.cellInf(1).posNum,poses));
else
    cellsToPlot=1:length(cExperiment.cellInf(1).posNum);
end


source1=[cExperiment.lineageInfo.motherInfo.(sources{1}).motherPosNum;cExperiment.lineageInfo.motherInfo.(sources{1}).motherLabel;cExperiment.lineageInfo.motherInfo.(sources{1}).motherTrap];
cellIndices1=find(ismember(source1',source1','rows'));

durs=zeros(length(cellIndices1),8); %8 is arbitrary max number of cell cyles before stress
durs(durs==0)=NaN;

for c=1:length(cellIndices1)
   birthTimesL1=sort(cExperiment.lineageInfo.motherInfo.(sources{1}).birthTime(cellIndices1(c),:));
   cellInfoL1=source1(:,cellIndices1(c));
   birthTimesL1(birthTimesL1==0)=[];
   birthTimesL1(birthTimesL1>stressTime)=[];
   
   for cc=1:length(birthTimesL1)
   
       if isempty(birthTimesL1)    
           continue
       
       elseif cc>(length(birthTimesL1)-1)
           continue    
       
       else          
            st=birthTimesL1(cc);
            en=birthTimesL1(cc+1);
            cycleDur=en-st;
            
       end 
       
       if any(cycleDur)
           durs(c,cc)=cycleDur; 
       end    
   end
end  
resultsStruct.durs=durs;
