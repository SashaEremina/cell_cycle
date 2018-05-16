function [sum]=nextavg(lasts, next)

%calculates average duration of 1st cell cycle during stress for all the
%cells with the same 'last' (same time from the stress to the previous
%reporter peak)

remove=isnan(lasts) | isnan(next);
lasts=last(~remove);
next=next(~remove);

sum=zeros(length(lasts),1);
sum(1)=next(1);
index=1;
[~,ind2,~]=unique(lasts);
refs=diff(ind2);
i=1;

for c=2:length(lasts)
 
    if isnan(lasts(c))
        break
    elseif lasts(c)==lasts(c-1)       
            if isnan(sum(index))
                sum(index)=0;
            end
            
            if ~isnan(next(c))     
                sum(index)=sum(index)+next(c);
            end 
            
            if c<length(lasts) && lasts(c)~=lasts(c+1) 
                sum(index)=sum(index)/refs(i);
                i=i+1;
            elseif c==length(lasts)
                sum(index)=sum(index)/refs(i);
                i=i+1;
            end    
    
    elseif lasts(c)~=lasts(c-1)
        index=index+1;
        sum(index)=next(c);
    end    
end