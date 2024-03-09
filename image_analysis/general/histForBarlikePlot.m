function [x,y] = histForBarlikePlot(bins,n)

    x = sort([bins bins(2:end-1)]);
    y = cat(1,n(1:end-1,:),n(1:end-1,:));
    
    y(1:2:end,:) = n(1:end-1,:);
    y(2:2:end,:) = n(1:end-1,:);
end