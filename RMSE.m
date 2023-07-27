function [rmse] = RMSE(y1,y2)

rmse = sqrt(nanmean((y1(:)-y2(:)).^2));

end

