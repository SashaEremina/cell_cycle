function plotMultipleLineage(cExperiment,cell_tracking_number,cell_position,trap_number,plot_field, sources)
% Inputs:
% cExperiment - object of class experimentTracking
% cell_tracking_number - vector of integers - tracking numbers of the cells
% to plot
% cell_position - vector of integers or cell array of strings. Position
% numbers or position names of the cells to plot
% trap_number - vector of integers - trap numbers of the cells to plot
% plot_field - string name of the field to plot
% for each cell (eg 'max5')
% sources - string or cell array of strings - sources of event data - ie
% fields of cExperiment.lineageInfo.motherInfo to be shown on the graph.
% Cell array should have one row per cell and one column per field to show
% eg sources={'htb2Edited', 'clbEdited' ; 'clbEdited' , 'HMMEdited'} means
% plot htb2Edited and clbEdited for cell 1 and clbEdited and HMMEdited for
% cell2

if ischar(cell_position)
    cell_position=find(strcmp(cExperiment.dirs,cell_position));
else
    if iscell(cell_position)
        cellpos=zeros(1,length(cell_position));
        for n=1:length(cell_position)
            cellpos(n)=find(strcmp(cExperiment.dirs,cell_position{n}));
        end
        cell_position=cellpos;
    end
end

smoothedOrNot=false;
plot_channel=3;%To plot GFP and mCherry change to [1 2]

if nargin<6
    sources={'htb2Edited'};
end

figure;

%% Get extracted data to plot
%Colourmap for labelling cell sources
cMap=[ 1 0 1 ; 0 0 1; 0 0 0 ; 1 0 1 ; 0 .5 0; 0 1 0 ;1 0 1 ;  ];
cNames={'magenta' 'blue' 'black' 'magenta' 'dark green'  'green' 'cyan' 'yellow' 'purple' };
legendText='';
positionNames={};
yMin=Inf;
yMax=0;
for n=1:length(cell_tracking_number)
    for ch=1:length(plot_channel)
        cell_data_index = (cExperiment.cellInf(1).posNum == cell_position(n)) &...
            (cExperiment.cellInf(1).trapNum == trap_number(n)) & ...
            (cExperiment.cellInf(1).cellNum == cell_tracking_number(n));
        cell_data(ch,:) = full(cExperiment.cellInf(plot_channel(ch)).(plot_field)(cell_data_index,:));
        xInc=cExperiment.cellInf(1).times(cell_data_index,:);
        xInc(cell_data==0)=[];
        cell_data(cell_data==0)=[];
        
        %% Plot the data
        if ~smoothedOrNot
            tData=cell_data(ch,:);
            tData(tData==0)=[];
            tDataMin=min(tData);
            tDataMax=max(tData);           
            tData=tData./tDataMin;
            positionName=cExperiment.dirs{cell_position(n)};
            positionName=positionName(1:end-4);
            positionNames{n}=positionName;
            plot(xInc,tData,'-','Color',cMap(n,:));
            
            yMin=min(yMin,min(tData));
            yMax=max(yMax,max(tData));
            
            
            ylim([yMin yMax]);
            %     xlim([min(xInc(:)) max(xInc(:))]);
        else
            [p,s,mu] = polyfit(1:numel(cell_data),cell_data,6);
            f_y = polyval(p,(1:numel(cell_data)),[],mu);
            dataToTest = cell_data - f_y;
            %Smooth the data
            
            
            dt=sgolayfilt(dataToTest,14,21);
            tData=dt;
            tData(tData==0)=[];
            yMin=min(tData);yMax=max(tData);
            tData=tData./yMax;
             positionName=cExperiment.dirs{cell_position(n)};
            positionName=positionName(1:end-4);
            plot(xInc,tData,'-','Color',cMap(n,:),'DisplayName',positionName);
            
            
            try
                %         ylim([yMin yMax]);
                %         xlim([min(xInc(:)) max(xInc(:))]);
            catch
            end
        end
        hold on;
        
    end
end
%legend(positionNames);

%% Add the peaks for two lineage sources


for n=1:length(cell_tracking_number)
    lineObject=zeros(1,length(sources));
    for s=1:length(sources)
        source=sources{s};
        
        if (strcmp('HMMEdited', source))
            legendText='Bud detected';
            sourceName='Bud detected';
            legendText=[legendText ': ' cNames{s} '   '];
            
        else
            legendText=[legendText source ': ' cNames{s} '   '];
            sourceName=source;
            
        end
        cell_mother_index = (cExperiment.lineageInfo.motherInfo.(source).motherPosNum == cell_position(n)) &...
            (cExperiment.lineageInfo.motherInfo.(source).motherTrap == trap_number(n)) & ...
            (cExperiment.lineageInfo.motherInfo.(source).motherLabel == cell_tracking_number(n));
        
        birth_times=cExperiment.lineageInfo.motherInfo.(source).birthTime(cell_mother_index,:);
        birth_times(birth_times==0) = [];
        ylim_plot = [yMin yMax+yMax*.01];
        for bi = 1:length(birth_times)            
            bt = birth_times(bi);
            if bi==1
                lineObject(s)=plot(xInc(bt)*[1 1],ylim_plot,'color',cMap(s+1,:),'DisplayName',source);
            end
            plot(xInc(bt)*[1 1],ylim_plot,'color',cMap(s+1,:),'DisplayName',source);
            set(gca,'Ylim',ylim_plot);
        end
        
                if isfield(cExperiment.lineageInfo.motherInfo.(source),'deathTime')
                    if length(cExperiment.lineageInfo.motherInfo.(source).deathTime)>=find(cell_mother_index)
                        deathTime=cExperiment.lineageInfo.motherInfo.(source).deathTime(cell_mother_index);
                        if deathTime>0
                            plot(xInc(round(deathTime))*[1 1],ylim_plot,'-b');
                        end
                    end
                end
        xInc=xInc+1;%to displace the lines slightly if they are on the same timepoint
    end
    legendPosX=get(gca,'Xlim');
    legendPosX=legendPosX(1)+0.05*legendPosX(2);
    legendPosY=get(gca,'Ylim');
    legendPosY=legendPosY(1)+0.95*(legendPosY(2)-legendPosY(1));
end
legend(positionNames);
    end
