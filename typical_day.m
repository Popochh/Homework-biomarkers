function [typical_day,typical_day_subj] = typical_day (activity_clean)

% The purpose of the function is to organize the activity in 24h, making a 
% typical day for each subject.
% 'typical_day_subj' contains in each column the 24H of a subject.
% 'typical day' contains only one column, derived by the mean of each 24H
% of the subjects.

N_subj = size(activity_clean,2);

daily_sample = 1440;

typical_day = zeros(24,1);
typical_day_subj = zeros(24, N_subj);
activity_clean_onesubj = mean(activity_clean,2,'omitnan');
temp_activity_clean_onesubj = reshape (activity_clean_onesubj,daily_sample,[]);

for i=1:24
     temp_activity_hour = temp_activity_clean_onesubj((i-1)*60+1: (i*60),:);
     typical_day(i) = mean(temp_activity_hour,'all','omitnan');

end


for i=1:N_subj
    temp_activity_clean = activity_clean(:,i);
    temp_activity_clean = reshape(temp_activity_clean,daily_sample,[]); % in each column we have the activity of i-th subject
    for j=1:24
        temp_activity_hour_subj = temp_activity_clean((j-1)*60+1: (j*60),:);
        typical_day_subj(j,i) = mean(temp_activity_hour_subj,'all','omitnan');
    end

end

% the algorithm does not suffers of Nan;


