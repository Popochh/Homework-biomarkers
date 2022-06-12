function [activity_clean,daily_mean, days] = mean_and_clean (activity, th)
% scrivi cosa fa la funtion ( media e toglie sottosoglis)
daily_sample = 24*60; % Number of elements in 24h
N_subj = size(activity,2); % Number of subject


%% Deleting the firt day, because of missing data in the first hours 

activity_new = activity((daily_sample+1):end,:);

% Temporary dimension of our dataset
days = size(activity_new,1)/(daily_sample);





%% Deleting the last day, because of missing data in the last hours

activity_clean = nan(size(activity_new));

for i=1:N_subj
    temp_activity = activity_new(:,i);
    pos = isnan(temp_activity);
    temp_activity(pos) = [];
    lunghezza = length(temp_activity);
    ultimo_giorno = floor(lunghezza/daily_sample);
    
    % Assign the last day as 'Nan'
    activity_new((ultimo_giorno*daily_sample+1):(ultimo_giorno+1)*daily_sample,i) = nan;

end

% Taking into account data without the last day
activity_clean = activity_new(1:(size(activity_new,1)-1440),:);

% Temporary dimension of our dataset without last day
days = size(activity_clean,1)/(daily_sample);

%% Evaluating daily mean of the activity

daily_mean = nan(days,N_subj);

for i = 1:N_subj

    temp_activity = activity_clean(:,i); % taking 1 subject into account
    
    
    for j=1:days

        temp_day = temp_activity((j-1)*daily_sample+1:(j*daily_sample));
        daily_mean(j,i) = mean(temp_day,'omitnan'); % pay attention to Nan

    end

end

clear temp_activity
clear temp_day;
clear i;
clear j;


%% Cleaning data under the Thresold, based on Daily activity
% the threshold is in input function

for i=1:N_subj
    temp_activity = activity_clean(:,i);
    temp_daily = daily_mean(:,i);
    pos = find(temp_daily<th); % finding days under the th
    
    if length(pos)>0
        
       % deleting data from the daily mean vector ( assign as Nan)
       daily_mean(pos,i) = nan;
    
        % deleting data, referring to the day under th, from the complete 
        % vector of activity ( assign as Nan)
    
        for j =1:length(pos)
    
            temp_delete = [((pos(j)-1)*daily_sample)+1:(pos(j)*daily_sample)]';
    
            
            temp_activity(temp_delete) = nan;
            activity_clean(:,i) = temp_activity;
    
        end
    end

end



clear i;
clear j;
clear temp_activity;
clear temp_daily;
clear temp_delete
clear pos;
