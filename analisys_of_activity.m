function [h_activity] = analisys_of_activity(dataset,n_subjects,name)
%% FUNCTION THAT HELPS IN THE DIVISION OF THE DIFFERENT DAYS 
% h_activity = is a struct in witch has been extracted for each subject
%              7 days of motor activity for 24 h each day

% dataset   = motor activity data that is being analyzed  
%             we will put the condition dataset or the control one
% n_subject = number of the patient or number of the control tha is being
%             analyzed by the for cycle
% name      = if the subject is a patient or a control

%% INPUT  take all the days


hours = 24;
week = 7;

% taking all the days and their specific code to perform the
% right assignment Monday, Tuesday ecc..

allDays = table2array(dataset(:,1));
allDays_unique = unique(yyyymmdd(allDays)); 
n_days = length(allDays_unique);

%% Eliminating who doesn't have 7 days in a row

if n_days <= 7
    disp(['L''analisi del soggetto: ' num2str(n_subjects) ' non si può fare'])
    disp('eliminalo !')
end

%% Transforming everything to serial codification 


% days in to code
serial_settimana = datenum((table2array(dataset(:,2))));
% eith unique I obtain the effective days of the analysis
serial_settimana_unique = unique(serial_settimana);

% i take monday = but is the start of the analysis not the real monday we
% will retrive it later
monday = serial_settimana_unique(1) + 1;
% Monday = end of the analysis
sunday = monday + 7;


% retriving the start and the end indexis
idx_mondays = find(serial_settimana == monday);
idx_monday = idx_mondays(1);
idx_sundays = find(serial_settimana == sunday);
idx_sunday = idx_sundays(1);

% trovo la durata dell'analisi che ora è in codice
% weekely_analisys = serial_settimana(idx_monday:idx_sunday);

%% AFTER THIS START AND END POINT, NOW I CHECK THE DAYS:

A = (allDays(idx_monday:1441:idx_sunday));
d = day(A,'name');

%% TAKING THE WEEK

sub = readtable([name '_' num2str(n_subjects) , '.csv'], 'Range','C2');

activity = table2array(sub);

weekely_analisys = activity(idx_monday:idx_sunday);

days = round(length(weekely_analisys)/7);


%% GOING DOWN TO THE DAYS

days_activity= zeros(days,week);
g = 1:days-1:length(weekely_analisys);

for i = 1:week  
    days_activity(:,i) = weekely_analisys(g(i):g(i+1));
end

x = categorical(d);


%% GOING DOWN TO THE HOURS

hours = 24;
week = 7;

% 60 idx of the length of the hours
idx_hours = round((length(weekely_analisys)/week)/hours); %60

os = 1:idx_hours:round(length(weekely_analisys)/7)-1;

weekely_hours = 1:idx_hours:length(weekely_analisys);

for i = 0:week-1  
    for j = 1:hours
        activity_time(j,i+1) = mean(weekely_analisys(os(j)+1440*i:(os(j)+60)+1440*i));      
    end
end


%% SAVE

h_activity = cell2table(num2cell(activity_time),'VariableNames',d);



end