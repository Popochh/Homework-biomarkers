function [day,night] = day_night ( activity_clean,wake_up,sleeping_time,days)

daily_sample = 24*60; % numero elementi in 24h
N_subj = size(activity_clean,2); % salvo numero di soggetti

% pensa di mettere tutto in strutture

%% 22:00 fino 05:59
% sveglio : 06:00 fino 21.59

% ci sono dei cicli di controllo nel caso in cui si vada a letto dopo la
% mezzanotte; questo perchè sto considerando un sonno di 7h




    % cambiando il numero in wake up, cambio l'ora in cui mi sveglio
    wake_up = wake_up*60;
    
    % considero 16 ore di giornata ( 24-8 per dormire)
    % 24-17 = 7 ore di sonno


    go_bed = wake_up+((24-sleeping_time)*60); % in questo caso qui sono le 22:00, e dormirò fino le 05:59
    
    daily_vector = (1:daily_sample)';

if (wake_up/60) <9
    
    nighttime_vector = [(1:wake_up-1)'; (go_bed:daily_sample)']; % da 00:00 a 05:59 e poi da 22:00 fino 23.59
    daytime_vector = [wake_up:go_bed-1]';

else
    go_bed = wake_up+(16*60)-daily_sample; % in questo caso qui sono le 22:00, e dormirò fino le 05:59

    nighttime_vector = [(wake_up-8*60: wake_up-1)'];
    daytime_vector = [(1:go_bed-1)';(wake_up:daily_sample)'];


end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% divido in nighttime e daytime
% 
day.activity = nan(size(activity_clean));
night.activity = nan (size(activity_clean));

% activity_daytime = [];
% activity_nighttime = [];

for i = 1:N_subj

    temp_activity = activity_clean(:,i); % NB ricorda che ci sono i NAN
    
    % assegno la media giornaliera a un vettore
    for j=1:days

        temp_daytime = temp_activity(((j-1)*daily_sample)+daytime_vector);
        day.mean_subj(j,i) = mean(temp_daytime,'omitnan'); % OSS è la media di ogni giorno
        day.activity(((j-1)*daily_sample)+daytime_vector,i) = temp_daytime;
        
        temp_nighttime = temp_activity(((j-1)*daily_sample)+ nighttime_vector);
        night.mean_subject(j,i) = mean(temp_nighttime,'omitnan'); % OSS è la media di ogni giorno
        night.activity(((j-1)*daily_sample)+nighttime_vector,i) = temp_nighttime;


    end

end


night.mean_single= mean(night.mean_subject,1,'omitnan'); % media di attivita notturna di ogni soggetto
% nighttime_mean_single2= mean(activity_nighttime,1,'omitnan');


day.mean_single= mean(day.mean_subj,1,'omitnan'); % media di attivita giorno di ogni soggetto
%daytime_mean_single= mean(activity_daytime,1,'omitnan');

clear temp_activity;
clear temp_daytime;
clear temp_nighttime;
clear i;
clear j;