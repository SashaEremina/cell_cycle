function resultsStruct=measureMultipleLineage2(cExperiment,poses, stressTime)

%measures the distance between the peaks and calculates it as a temporal
%fraction of the cell cycle taking the adjacent datapoints in the lineage2
%as the begining and the end of each cell cycle

%corresponds to two sliding window search method in Sasha's report
%(increments in both lineages)

%sources have to be defined in accordance with positions. if there is any
%discrepancy between the two - sources will be used as a reference


sources={'whi5Edited', 'HMMEdited'};
gapThreshold=25;

if nargin>1
    if ischar(poses)
        poses=find(~cellfun(@isempty,strfind(cExperiment.dirs,poses)));        
    end
    cellsToPlot=find(ismember(cExperiment.cellInf(1).posNum,poses));
else
    cellsToPlot=1:length(cExperiment.cellInf(1).posNum);
end

source1=[cExperiment.lineageInfo.motherInfo.(sources{1}).motherPosNum;cExperiment.lineageInfo.motherInfo.(sources{1}).motherLabel;cExperiment.lineageInfo.motherInfo.(sources{1}).motherTrap];
source2=[cExperiment.lineageInfo.motherInfo.(sources{2}).motherPosNum;cExperiment.lineageInfo.motherInfo.(sources{2}).motherLabel;cExperiment.lineageInfo.motherInfo.(sources{2}).motherTrap];


cellIndices1=find(ismember(source1',source2','rows'));
cellIndices2=find(ismember(source2',source1','rows'));



eventGaps=zeros(length(cellIndices2), 18);%18 is arbitrary max number of cell cyles
for c=1:length(cellIndices1)
   birthTimesL1=sort(cExperiment.lineageInfo.motherInfo.(sources{1}).birthTime(c,:));
   cellInfoL1=source1(:,cellIndices1(c));
   indexL2=find(source2(1,:)==cellInfoL1(1) & source2(2,:)==cellInfoL1(2) & source2(3,:)==cellInfoL1(3));
   birthTimesL2=sort(cExperiment.lineageInfo.motherInfo.(sources{2}).birthTime(indexL2,:));
   %First fluorescence peak - lineage source 2
   birthTimesL1=birthTimesL1(birthTimesL1>=min(birthTimesL2(birthTimesL2>0)));
   birthTimesL1(birthTimesL1==0)=[];
   birthTimesL2(birthTimesL2==0)=[];
   
   %To only look at the subset of the experiment (e.g. to compare peak
   %distance before and after stress)- otherwise comment out 
   
   birthTimesL1(birthTimesL1>stressTime)=[];
   birthTimesL2(birthTimesL2>stressTime)=[];
   
   %Loop through the cell cycles
   ccL1=1;
   ccL2=2;
   previous=0;
   
   for cc=1:length(birthTimesL2)
       birthTimesL1
       birthTimesL2
       ccL1
       ccL2
      
   
       if isempty(birthTimesL1) || isempty(birthTimesL2)
           continue
           
       else    
           
           if birthTimesL2(ccL2-1)<=birthTimesL2(ccL1)<=birthTimesL2(ccL2)
               
               eventGap=birthTimesL2(ccL2)-birthTimesL1(ccL1);
               cycleDurs=birthTimesL2(ccL2)-birthTimesL2(ccL2-1);
               ratio=eventGap/cycleDurs;
      
           end
           
            if eventGap>cycleDurs && ccL1<=length(birthTimesL1-1)
                   ccL1=ccL1+1;
                   if ccL2<=length(birthTimesL2) && ccL1<=length(birthTimesL1)
                       eventGap=birthTimesL2(ccL2)-birthTimesL1(ccL1);
                   else
                       break;
                   end

            elseif eventGap<gapThreshold && ccL1<=length(birthTimesL1-1)
                   ccL2=ccL2+1;
                   if ccL2<=length(birthTimesL2) && ccL1<=length(birthTimesL1)
                   eventGap=birthTimesL2(ccL2)-birthTimesL1(ccL1);
                   else
                       break;
                   end
             end
         ccL1=ccL1+1;
         ccL2=ccL2+1;
           if ccL1>length(birthTimesL1) || ccL2>length(birthTimesL2-1)
               break
           end     
           eventGaps(c,cc)=eventGap;
       end
       end
end
resultsStruct.eventGaps=eventGaps;

average=zeros(length(resultsStruct.eventGaps),1);

for i=1:length(resultsStruct.eventGaps)
    average(i)=sum(abs(resultsStruct.eventGaps(i,:))) ./ length(nonzeros(resultsStruct.eventGaps(i,:)));    
end
 ans=nanmean(average);
 std=nanstd(average);


