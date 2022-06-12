%% BIOMARKERS, PRECISION MEDICINE AND DRUG DEVELOPMENT, HOMEWORK

% CARRANZA MELLANA M.
% BAGANTE A.
% ALUSHANI T.
% LADDOMADA P.


clear all
close all
clc

%%
N_cond = 23;
N_contr = 32;

% number of elements in 24 hours (one sample of activity each minute)
daily_sample = 24*60;
%% Add path to load data
base_path=pwd; 
addpath(genpath(fullfile(base_path, 'DepressionDataset')))
%% Evaluating Age and Sex of the two population

scores = readtable('scores.csv');
patients_sex.condition_gender = table2array( scores(1:23,'gender'));
patients_age.condition_age = table2array (scores(1:23,'age'));

patients_sex.control_gender =table2array( scores(24:end,'gender'));
patients_age.control_age = table2array (scores(24:end,'age'));

% Mean and Std of condition patients' age
for i=1: length(patients_age.condition_age)
    temp= patients_age.condition_age(i);
    aux = cell2mat(temp);
    age_aux(i) = mean([str2double(aux(:,1:2)),str2double(aux(:,4:5))]);

    clear temp;
    clear aux;

end

patients_age.average_age_condition = mean(age_aux);
patients_age.sd_age_condition = std(age_aux);
clear age_aux;

% Mean and Std of control patients' age
for i=1: length(patients_age.control_age)
    temp= patients_age.control_age(i);
    aux = cell2mat(temp);
    age_aux(i) = mean([str2double(aux(:,1:2)),str2double(aux(:,4:5))]);

    clear temp;
    clear aux;

end

patients_age.average_age_control = mean(age_aux);
patients_age.sd_age_control = std(age_aux);
clear age_aux;

% Patients's Sex (1= female, 2 = male)

patients_sex.male_condition = length(find(patients_sex.condition_gender==2));
patients_sex.female_condition = length(find(patients_sex.condition_gender==1));


patients_sex.male_control = length(find(patients_sex.control_gender==2));
patients_sex.female_control = length(find(patients_sex.control_gender==1));




%% VISUAL INSPECTION
% we just played with the data to see if we can just visualize some
% differences
% We decided to take a week a for each subject and compare the activity
% between the different days
% The depressed patient week vs The control patient week
% We took the raw data eliminating the first day from the analysis and then
% retrive the different days of the week in which they have been analysed

addpath(genpath('DepressionDataset'))

% using the function analisys_of_activity we extract the activity of the
% patients in the two datasets, 
% the activity extracted 24 h for 7 days, the days are not in order yet

% CONTROL
for c1 = 1:32
    control = readtable(['control_' num2str(c1) '.csv']);
    h_activity_control(c1).patient = analisys_of_activity(control,c1,'control');
end

% CONDITION
for c2 = 1:23
    condition = readtable(['condition_' num2str(c2) '.csv']);
    h_activity_condition(c2).patient = analisys_of_activity(condition,c2,'condition');
end


clear c1
clear c2
clear condition
clear control
%%
% here we create a new struct in order to divide the week from the weekend
% and retrive the days of the week in the right order

% CONTROL
for i1 = 1 : 32
    weekend_control(i1).patient= [h_activity_control(i1).patient(:,'Friday'),h_activity_control(i1).patient(:,'Saturday'),h_activity_control(i1).patient(:,'Sunday')];
    week_control(i1).patient = [h_activity_control(i1).patient(:,'Monday'),h_activity_control(i1).patient(:,'Tuesday'),h_activity_control(i1).patient(:,'Wednesday'),h_activity_control(i1).patient(:,'Thursday')];
end

% CONDITION
for i2 = 1 : 23
    weekend_condition(i2).patient= [h_activity_condition(i2).patient(:,'Friday'),h_activity_condition(i2).patient(:,'Saturday'),h_activity_condition(i2).patient(:,'Sunday')];
    week_condition(i2).patient = [h_activity_condition(i2).patient(:,'Monday'),h_activity_condition(i2).patient(:,'Tuesday'),h_activity_condition(i2).patient(:,'Wednesday'),h_activity_condition(i2).patient(:,'Thursday')];
end

clear i1
clear i2
clear condition
clear control
clear h_activity_condition
clear h_activity_control


%% Here we extrapolate the activity count in order to visualize it 
% as difference of activity between depressed and patients days

% CONTROL
h_sum_control = zeros(24,7);
for p = 1:32
    for h = 1:24
        % week
        for w1 = 1:4
            h_sum_control(h,w1) = h_sum_control(h,w1) + table2array(week_control(p).patient(h,w1));
        end
        % weekend
        for w2 = 1:3
            h_sum_control(h,w1+w2) = h_sum_control(h,w1+w2) + table2array(weekend_control(p).patient(h,w2));
        end
    end
end
h_avg_control = h_sum_control/p;
clear h_sum_control

% CONDITION
h_sum_condition = zeros(24,7);
for p = 1:23
    for h = 1:24
        % week
        for w1 = 1:4
            h_sum_condition(h,w1) = h_sum_condition(h,w1) + table2array(week_condition(p).patient(h,w1));
        end
        % weekend
        for w2 = 1:3
            h_sum_condition(h,w1+w2) = h_sum_condition(h,w1+w2) + table2array(weekend_condition(p).patient(h,w2));
        end
    end
end
h_avg_condition = h_sum_condition/p;

clear h_sum_condition
clear weekend_condition
clear weekend_control
clear week_condition
clear week_control
clear p
clear h
clear w1
clear w2

%% Weekely activity
% Creation of 2 graphs
% Differences between the depressed patients and the healthy subject days 
% Differences between the depressed patients and the healthy subject week 

week_labels = {'Monday','Tuesday','Wednesday','Thoursday','Friday','Saturday','Sunday'};

Subject = [ones(size(h_avg_control,1),1); zeros(size(h_avg_condition,1),1)];
Subject = categorical(Subject, [0 1],{'Condition','Control'});

h_avg = [h_avg_control; h_avg_condition];

t = tiledlayout(2,4);

for w=1:7
    ax = ['ax' num2str(w)];
    ax = nexttile;
    boxchart(ax, h_avg(:,w), 'GroupByColor',Subject)
    legend
    title(week_labels(:,w))
end

clear w
clear h_avg
clear week_labels

ylabel(t,'Activity Count')
figure,
tiledlayout(1,2);
nexttile
imagesc(h_avg_control)
% colorbar
title('average weak of an healthy control')
ylabel('Activity Count')

nexttile
imagesc(h_avg_condition)
% colorbar[0 600]
title('average weak of a depressed patient')
cb = colorbar;
cb.Layout.Tile = 'east';
colormap jet


clear h_avg_control
clear h_avg_condition
clear ax
clear cb
clear t
clear subject

%% END OF THE VISUALIZATION AND BEGGING OF THE DATA ANALYSIS

%% Creating vector to assign the Activity of each Condition subject

length_max_cond =0;


for i=1:N_cond

    % Evaluating max length
    name = ("condition_" + num2str(i)) + '.csv';
    data = table2array(readtable(name,'Range','C2'));
    time = (readtable(name,'Range','A2'));
    time = table2array(time(:,1));

    first_minute_cond(i) = hour(time(1))*60 + minute(time(1));
    last_minute_cond(i) = hour(time(end))*60 + minute(time(end));

    aux_begin = nan(first_minute_cond(i)-1,1);
    aux_close = nan(24*60-last_minute_cond(i),1);

    temp_activity = [aux_begin;data;aux_close];
    if length(temp_activity)> length_max_cond
        length_max_cond = length(temp_activity);

    end

end

% Creating the vector to contain the Activity of the patients
activity_condition = nan(length_max_cond,N_cond);

clear time;
clear temp_activity;
clear aux_close;
clear aux_begin;
clear data;
clear i;
clear name;


% Assign Activity Condition


for i=1:N_cond

    name = ("condition_" + num2str(i)) + '.csv';
    data = table2array(readtable(name,'Range','C2'));
    time = (readtable(name,'Range','A2'));
    time = table2array(time(:,1));

    aux_begin = nan(first_minute_cond(i)-1,1);
    temp_activity = [aux_begin;data];
    aux_close = nan(length_max_cond- length(temp_activity),1);

    temp_activity = [aux_begin;data;aux_close];

    activity_condition(:,i) = temp_activity;


end

clear time;
clear temp_activity;
clear aux_close;
clear aux_begin;
clear data;
clear i;
clear name;



%% Creating vector to assign the Activity of each Control subject
length_max_contr = 0;


for i=1:N_contr

    % Evaluating max length
    name = ("control_" + num2str(i)) + '.csv';
    data = table2array(readtable(name,'Range','C2'));
    time = (readtable(name,'Range','A2'));
    time = table2array(time(:,1));

    first_minute_contr(i) = hour(time(1))*60 + minute(time(1));
    last_minute_contr(i) = hour(time(end))*60 + minute(time(end));

    aux_begin = nan(first_minute_contr(i)-1,1);
    aux_close = nan(24*60-last_minute_contr(i),1);

    temp_activity = [aux_begin;data;aux_close];
    if length(temp_activity)> length_max_contr
        length_max_contr = length(temp_activity);

    end

end

% Creating the vector to contain the Activity of the patients 
activity_control = nan(length_max_contr,N_contr);

clear time;
clear temp_activity;
clear aux_close;
clear aux_begin;
clear data;
clear i;
clear name;


% Assign Activity Control


for i=1:N_contr

    name = ("control_" + num2str(i)) + '.csv';
    data = table2array(readtable(name,'Range','C2'));
    time = (readtable(name,'Range','A2'));
    time = table2array(time(:,1));

    aux_begin = nan(first_minute_contr(i)-1,1);
    temp_activity = [aux_begin;data];
    aux_close = nan(length_max_contr- length(temp_activity),1);

    temp_activity = [aux_begin;data;aux_close];

    activity_control(:,i) = temp_activity;


end

clear time;
clear temp_activity;
clear aux_close;
clear aux_begin;
clear data;
clear i;
clear name;


%% Cleaning unuseful variables

clear first_minute_contr;
clear first_minute_cond;
clear last_minute_contr;
clear last_minute_cond;
clear length_max_contr;
clear length_max_cond;



%% Cleaning the data - Creating Daily mean activity
% A daily mean under the treshold will be assigned as Nan, using 
% 'mean_and_clean' function (go inside the function for more info)

th = 40; % thershold

% Condition
[activity_clean_cond,daily_mean_cond, days_cond] = mean_and_clean (activity_condition, th);
    
% Control
[activity_clean_contr,daily_mean_contr, days_contr] = mean_and_clean (activity_control, th);





%% Evaluating a typical day of Control and Condition patients

[typical_day_cond,typical_day_cond_subj] = typical_day (activity_clean_cond);
[typical_day_contr,typical_day_contr_subj] = typical_day (activity_clean_contr);



% Plot, comparing Condition vs Control subject
figure; plot(typical_day_cond,'--r');
hold on
plot (typical_day_contr,'k');
grid on;
xlabel ('Time [hours]');
ylabel ('Activity Count');
legend(['Condition'],['Control']);
title (' Typical day Condition vs Control')



%% Dividing the activity in Night-time and Day-time

% Based on literature, considering Night time from 23.00 to 06.00

wake_up = 6; % time to wake up
sleeping_time = 7; % how many hours of sleep?

[day_cond,night_cond] = day_night ( activity_clean_cond,wake_up,sleeping_time,days_cond);

[day_contr, night_contr] = day_night ( activity_clean_contr,wake_up,sleeping_time,days_contr);



%% All Day
% verifing Normality Distribution
% Comparison test between 2 distibution


daily_mean_contr_reshape = reshape (daily_mean_contr,[],1);
pos = isnan(daily_mean_contr_reshape);
daily_mean_contr_reshape(pos) = [];

[h,p,kstat,critval] = lillietest(daily_mean_contr_reshape) % h=0


%
daily_mean_cond_reshape = reshape (daily_mean_cond,[],1);
pos = isnan(daily_mean_cond_reshape);
daily_mean_cond_reshape(pos) = [];

[h,p,kstat,critval] = lillietest(daily_mean_cond_reshape) % h=1


% Mann-Whithney Test
[p,h,stats] = ranksum(daily_mean_contr_reshape,daily_mean_cond_reshape)% h=1, p = 3.8735e-32 <0.001

mean(daily_mean_contr_reshape)
std(daily_mean_contr_reshape)

mean(daily_mean_cond_reshape)
std(daily_mean_cond_reshape)


% Plot the results ( All day)
daily_mean_visual = [daily_mean_cond_reshape; nan(length(daily_mean_contr_reshape)-length(daily_mean_cond_reshape),1)];
vect_daily_mean = [daily_mean_visual daily_mean_contr_reshape];
str_d = {'p-value <0.001 *'};
x = [ones(size(daily_mean_contr_reshape,1),1) 2*ones(size(daily_mean_contr_reshape,1),1)];
xx = string(x);
xx(xx=="1") = "Condition";
xx(xx=="2") = "Control";
xx = categorical(xx);

figure
swarmchart(xx,vect_daily_mean); 
hold on
plot(nanmean(vect_daily_mean,1),'-ok')
hold off
legend(["Condition","Control","Mean"])


text(0.5,700,str_d)
ylim([-10 750])
ylabel('Activity Count')
title('Daily activity ')

clear daily_mean_visual
clear vect_daily_mean
clear str_d
clear x
clear xx



%% Daytime
%
day_contr.mean_reshape = reshape (day_contr.mean_subj,[],1);
pos = isnan(day_contr.mean_reshape);
day_contr.mean_reshape(pos) = [];
[h,p,kstat,critval] = lillietest(day_contr.mean_reshape) % h=0


%
day_cond.mean_reshape = reshape (day_cond.mean_subj,[],1);
pos = isnan(day_cond.mean_reshape);
day_cond.mean_reshape(pos) = [];

[h,p,kstat,critval] = lillietest(day_cond.mean_reshape) % h=1 not Gaussian


% Mann-Whitney test
[p,h,stats] = ranksum(day_cond.mean_reshape,day_contr.mean_reshape)% h=1, p = 1.7415e-30

mean(day_contr.mean_reshape)
std(day_contr.mean_reshape)

mean(day_cond.mean_reshape)
std(day_cond.mean_reshape)


%% Nighttime
% the vector contains the mean values of the night activity for each day

night_contr.mean_reshape = reshape (night_contr.mean_subject,[],1);
pos = isnan(night_contr.mean_reshape);
night_contr.mean_reshape(pos) = [];
[h,p,kstat,critval] = lillietest(night_contr.mean_reshape) % h=1


%
night_cond.mean_reshape = reshape (night_cond.mean_subject,[],1);
pos = isnan(night_cond.mean_reshape);
night_cond.mean_reshape(pos) = [];
[h,p,kstat,critval] = lillietest(night_cond.mean_reshape) % h=1

%
[p,h,stats] = ranksum(night_cond.mean_reshape,night_contr.mean_reshape) % h=1, p = 1.6117e-08

mean(night_contr.mean_reshape)
std(night_contr.mean_reshape)

mean(night_cond.mean_reshape)
std(night_cond.mean_reshape)
% Subplot between Daytime and Nighttime distribution
figure;

% Day-time

subplot 121

daytime_mean_visual = [day_cond.mean_reshape; nan(length(day_contr.mean_reshape)-length(day_cond.mean_reshape),1)];
vect_daily_mean = [daytime_mean_visual day_contr.mean_reshape];
str_all = {'p-value < 0.001*'}; 
x = [ones(size(day_contr.mean_reshape,1),1) 2*ones(size(day_contr.mean_reshape,1),1)];
xx = string(x);
xx(xx=="1") = "Condition";
xx(xx=="2") = "Control";
xx = categorical(xx);


swarmchart(xx,vect_daily_mean)
hold on
plot(nanmean(vect_daily_mean,1),'-ok')
hold off
legend(["Condition","Control","Mean"])
text(0.5,700,str_all)
% text(0.1,750,{'A'})
ylim([-10 750])
% xlim([0 3])
ylabel('Activity')
title('Daytime activity ')

clear daytime_mean_visual
clear vect_daily_mean
clear str_all
clear x
clear xx



% Night-time

nigth_mean_visual = [night_cond.mean_reshape; nan(length(night_contr.mean_reshape)-length(night_cond.mean_reshape),1)];
vect_nigth_mean = [nigth_mean_visual night_contr.mean_reshape];
str_n = {'p-value < 0.001 *'};
x = [ones(size(vect_nigth_mean,1),1) 2*ones(size(vect_nigth_mean,1),1)];
xx = string(x);
xx(xx=="1") = "Condition";
xx(xx=="2") = "Control";
xx = categorical(xx);


subplot 122

swarmchart(xx,vect_nigth_mean) 
hold on
plot(nanmean(vect_nigth_mean,1),'-ok')
hold off
legend(["Condition","Control",'Mean'])
text(0.5,700,str_n)
ylim([-10 750])
% xlim([0 3])
ylabel('Activity Count')
title('Nigthtime activity ')


clear nigth_mean_visual
clear vect_nigth_mean
clear str_n
clear x



%% Evaluating the differences of activity between day-time and night-time
difference_cond = day_cond.mean_single- night_cond.mean_single;
difference_contr = day_contr.mean_single- night_contr.mean_single;

mean (difference_cond,'omitnan')
std (difference_cond,'omitnan')

mean (difference_contr,'omitnan')
std (difference_contr,'omitnan')


% Evaluating Normal distribution
[h,p,kstat,critval] = lillietest(difference_contr) % h=0 %p_value = 0.3687 Not statistically significant
[h,p,kstat,critval] = lillietest(difference_cond) % h=0 %p_value = 0.5 

% Mann-Whitney Test
[p,h,stats] = ranksum(difference_contr,difference_cond) % h=1 p_value = 9.0414e-04


%
vect_diff = [[difference_cond, nan(1,(N_contr-N_cond))]; difference_contr]';
str_diff = {'p-value = 2.4969e-04 * '};
x = [ones(size(N_contr,2),1) 2*ones(size(N_contr,2),1)];
xx = string(x);
xx(xx=="1") = "Condition";
xx(xx=="2") = "Control";
xx = categorical(xx);


aux = mean(vect_diff,'omitnan');

vect_diff_cond = aux(:,1);
vect_diff_contr = aux(:,2);

vect_daytime_cond = mean(day_cond.mean_single,'omitnan')
vect_daytime_contr = mean(day_contr.mean_single,'omitnan')

vect_nighttime_cond = mean(night_cond.mean_single,'omitnan')
vect_nighttime_contr = mean(night_contr.mean_single,'omitnan')

vect_24h_cond = mean(daily_mean_cond_reshape,'omitnan')
vect_24h_contr = mean(daily_mean_contr_reshape,'omitnan');


% Barplot
figure;

str_diff = {'p-value < 0.001*'};

xx2 = ([xx(1), xx(1), xx(1), xx(1); ...
    xx(2), xx(2), xx(2), xx(2)])
bar(xx2,[vect_24h_cond, vect_daytime_cond, vect_nighttime_cond, vect_diff_cond ; ...
    vect_24h_contr, vect_daytime_contr, vect_nighttime_contr, vect_diff_contr])

ylabel('Activity Count')
text(0.25,290,str_diff)

title('Comparison between Condition and Control ')
legend(["Daily Activity","Day-time","Night-time","Difference"],'Location','northwest')
 
%%
clear daytime_mean_visual
clear vect_diff
clear str_diff
clear aux
clear xx
clear xx2
clear x




%%
%% PART 1 - GROUP DIFFERENCES
%% are the motor activity or any derived scores associated to depression score severity?

%% <> <> <> <> <> <> <> <> <> <> CORRELATIONS <> <> <> <> <> <> <> <> <> <> 



%% 
% madrs1 -> severity depression score assessed before the analysis 
% madrs2 -> severity depression score assessed after the analysis 
% for the correlation analisys it was used this scores separatly and then their
% difference and than their mean

% MILD DEPRESSION:    7-19
% DEPRESSION:        20-34
% SEVERE DEPRESSION: 35-more

daily_mean = mean(daily_mean_cond,1,'omitnan'); % we will use only the condition vector

madrs1 = str2num(cell2mat(table2array( scores(1:23,'madrs1') )));
madrs2 = str2num(cell2mat(table2array( scores(1:23,'madrs2') )));

dMADRS = madrs2 - madrs1;

%% TAKING ALL THE COVARIATES

marriage = str2num(cell2mat(table2array( scores(1:23,'marriage') )));
work = str2num(cell2mat(table2array( scores(1:23,'work') )));
gender = table2array( scores(1:23,'gender') );
edu = [0 0 0 1 1 0 1 1 0 0 0 0 1 0 1 1 0 1 2 0 0 0 2]'; % created manually: a 0 was insert for the 22th record
afftype = str2num(cell2mat(table2array( scores(1:23,'afftype') )));
% melanch = table2array( scores(1:23,'melanch') )); this covariate was not considered in the analysis
inpatient = str2num(cell2mat(table2array( scores(1:23,'inpatient') )));
age =   [4  5  6  2  7  4  1  2  6  6  6  5  3  9  8  6  7  5  7  3  4  10 3]'; % created manually



%% fit with linear regression
% To asses manually the weights of the covariates, we manually removed the
% less significant iteratively to assess the values
% We stopped when the results were worsening and then we decided the best
% model watching the result and also using also BIC and AIC

avgMADRS = round((madrs1+madrs2)/2);

tbl = table(daily_mean',work,marriage,gender,afftype,inpatient,edu,age,madrs1,madrs2,dMADRS,avgMADRS, 'VariableNames',{'Dailymean','Work','Marriage','Gender','Afftype','Inpatient','Edu','Age','Madrs1','Madrs2','dMADRS','avgMADRS'});


% MADRS1
lm1 = fitlm(tbl,'Madrs1~Dailymean+Work+Marriage+Gender+Afftype+Inpatient+Edu+Age'); % R2 = 0.562 || p-value = 0.0892
lm2 = fitlm(tbl,'Madrs1~Dailymean+Marriage+Inpatient+Edu+Age');                     % R2 = 0.542 || p-value = 0.0137
lm2.ModelCriterion                                                                  %% --> AIC = 130.4 | BIC = 137.2
lm3 = fitlm(tbl,'Madrs1~Dailymean+Marriage+Inpatient+Edu');                         % R2 = 0.534 || p-value = 0.00598  %% Best
lm3.ModelCriterion                                                                  %% --> AIC = 128.8 | BIC = 134.5   %% confirmed
lm4 = fitlm(tbl,'Madrs1~Marriage+Inpatient+Edu');                                   % R2 = 0.475 || p-value = 0.00579
lm4.ModelCriterion                                                                  %% --> AIC = 129.6 | BIC = 138.1





% MADRS2 
% We did start the analysis with also Madrs1 but then removed because of
% its obvious dependency on Madrs2
lm1 = fitlm(tbl,'Madrs2~Dailymean+Work+Marriage+Gender+Afftype+Inpatient+Edu+Age+Madrs1'); % R2 = 0.595
lm1 = fitlm(tbl,'Madrs2~Dailymean+Work+Marriage+Gender+Afftype+Inpatient+Edu+Age');        % R2 = 0.595 
lm2 = fitlm(tbl,'Madrs2~Dailymean+Work+Gender+Inpatient+Edu+Age');                         % R2 = 0.338 
lm3 = fitlm(tbl,'Madrs2~Dailymean+Inpatient+Edu+Age');                                     % R2 = 0.296
lm4 = fitlm(tbl,'Madrs2~Inpatient');                                                       % R2 = 0.168

% dMADR
lm1 = fitlm(tbl,'dMADRS~Dailymean+Work+Marriage+Gender+Afftype+Inpatient+Edu+Age'); % R2 = 0.381 || p-value = 0.431
lm2 = fitlm(tbl,'dMADRS~Work+Marriage+Gender+Afftype');                             % R2 = 0.322 || p-value = 0.118


% average MADRS
lm1 = fitlm(tbl,'avgMADRS~Dailymean+Work+Marriage+Gender+Afftype+Inpatient+Edu+Age'); % R2 = 0.504 || p-value = 0.166
lm2 = fitlm(tbl,'avgMADRS~Dailymean+Marriage+Inpatient+Edu+Age');                     % R2 = 0.499 || p-value = 0.0264    
lm2.ModelCriterion                                                                    %% --> AIC = 128.0 | BIC = 134.8    
lm3 = fitlm(tbl,'avgMADRS~Dailymean+Marriage+Inpatient+Edu');                         % R2 = 0.468 || p-value = 0.0178    
lm3.ModelCriterion                                                                    %% --> AIC = 127.4 | BIC = 138.1    


% --> The best model
lm3 = fitlm(tbl,'Madrs1~Dailymean+Marriage+Inpatient+Edu');                         % R2 = 0.536 || p-value = 0.00576  
figure, plot(lm3)

%% QUANTITAVELY MEASURE THE CORRELATION VALUES BETWEEN THE COVARIATES AND THE MADRS SCORE
% WE RUN IT IN A VISUAL WAY TO ASSES DIRECTLY IF THERE WERE SIGNIFICANT
% RESULTS TO REPORT
% It was done for some variables a previous check to understand if we had
% to use a Pearson's or a Spearman's Correlation



%% - MARRIAGE 

[rho pval_dMADRS] = corr(dMADRS, marriage,'type','Spearman')
R2_dMADRS = rho.^2
[rho pval_madrs1] = corr(madrs1, marriage,'type','Spearman')
R2_MADRS1 = rho.^2
[rho pval_madrs2] = corr(madrs2, marriage,'type','Spearman')
R2_MADRS2 = rho.^2
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, marriage,'type','Spearman')
R2_meanMADRS = rho.^2


if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end
%% GENDER

[h,p] = lillietest(gender)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, gender,'type','Spearman')
R2_dMADRS = rho.^2 
[rho pval_madrs1] = corr(madrs1, gender,'type','Spearman')
R2_MADRS1 = rho.^2
[rho pval_madrs2] = corr(madrs2, gender,'type','Spearman')
R2_MADRS2 = rho.^2
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, gender,'type','Spearman')
R2_meanMADRS = rho.^2

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end
%% DAILY ACTIVITY

[h,p] = lillietest(daily_mean)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, daily_mean','type','Spearman')
R2_dMADRS = rho.^2; 
[rho pval_madrs1] = corr(madrs1, daily_mean','type','Spearman')
R2_MADRS1 = rho.^2;
[rho pval_madrs2] = corr(madrs2, daily_mean','type','Spearman')
R2_MADRS2 = rho.^2;
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, daily_mean','type','Spearman')
R2_meanMADRS = rho.^2;

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end
%% WORK

[h,p] = lillietest(work)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, work,'type','Spearman')
R2_dMADRS = rho.^2; 
[rho pval_madrs1] = corr(madrs1, work,'type','Spearman')
R2_MADRS1 = rho.^2;
[rho pval_madrs2] = corr(madrs2, work,'type','Spearman')
R2_MADRS2 = rho.^2;
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, work,'type','Spearman')
R2_meanMADRS = rho.^2;

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end

%% AFFTYPE

[h,p] = lillietest(afftype)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, afftype,'type','Spearman')
R2_dMADRS = rho.^2; 
[rho pval_madrs1] = corr(madrs1, afftype,'type','Spearman')
R2_MADRS1 = rho.^2;
[rho pval_madrs2] = corr(madrs2, afftype,'type','Spearman')
R2_MADRS2 = rho.^2;
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, afftype,'type','Spearman')
R2_meanMADRS = rho.^2;

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end

%% INPATIENT significant value with madrs1 and meanMADRS 

[h,p] = lillietest(inpatient)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, inpatient,'type','Spearman')
R2_dMADRS = rho.^2 
[rho pval_madrs1] = corr(madrs1, inpatient,'type','Spearman')
R2_MADRS1 = rho.^2 % R2 = 0.3504                      
[rho pval_madrs2] = corr(madrs2, inpatient,'type','Spearman')
R2_MADRS2 = rho.^2
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, inpatient,'type','Spearman')
R2_meanMADRS = rho.^2 % R2 = 0.2297

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end
%% EDU 

[h,p] = lillietest(edu)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, edu,'type','Spearman')
R2_dMADRS = rho.^2; 
[rho pval_madrs1] = corr(madrs1, edu,'type','Spearman')
R2_MADRS1 = rho.^2; 
[rho pval_madrs2] = corr(madrs2, edu,'type','Spearman')
R2_MADRS2 = rho.^2;
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, edu,'type','Spearman')
R2_meanMADRS = rho.^2; 

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end
%% AGE

[h,p] = lillietest(age)

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

[rho pval_dMADRS] = corr(dMADRS, age,'type','Spearman')
R2_dMADRS = rho.^2; 
[rho pval_madrs1] = corr(madrs1, age,'type','Spearman')
R2_MADRS1 = rho.^2;
[rho pval_madrs2] = corr(madrs2, age,'type','Spearman')
R2_MADRS2 = rho.^2;
[rho pval_meanmadrs] = corr((madrs1+madrs2)/2, age,'type','Spearman')
R2_meanMADRS = rho.^2;

if pval_madrs1 <= 0.05 
    disp('significativo con madrs1')
elseif pval_madrs2 <= 0.05 
    disp('significativo con madrs2')
elseif pval_dMADRS <= 0.05 
    disp('significativo con madrs2')
elseif pval_meanmadrs <= 0.05 
    disp('significativo con madrs2')
end


tbl1 = table(inpatient,madrs1,'VariableNames',{'Inpatient','Madrs1'});
figure;
subplot(121)
[R pValue] =  corrplot(tbl1,'type','Spearman') 
title('Correlation matrix between Inpatient and Madrs1')

tbl2 = table(inpatient,avgMADRS,'VariableNames',{'Inpatient','avgMADRS'});
subplot(122)
[R pValue] =  corrplot(tbl2,'type','Spearman') 
title('Correlation matrix between Inpatient and avgMadrs')

%% <<<<<<<<<<<<<<<<<<<<<<<<<<<    PREDICTION   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% AS VARIABLES THAT COULD HELP US IN THE PREDICTION WE CHOSE:
% THE MEANS OF THE DAILY ACTIVITY, THE STANDARD DEVIATION AND THE COUNT OF
% THE ZEROS
% WE USED THE MATLAB TOOLBOX TO TAKE THE BEST PREDICTION MODEL

%% 
perc_cond = zeros(N_cond,1);


for i = 1:N_cond
    tmp_activity = activity_clean_cond(:,i);
    Nan_pos = isnan(tmp_activity);
    tmp_activity(Nan_pos) = [];
    zero_n = length(find(tmp_activity == 0));
    perc_cond(i) = zero_n/length(tmp_activity);
end


perc_contr = zeros(N_contr,1);


for j = 1:N_contr
    tmp_activity = activity_clean_contr(:,j);
    Nan_pos = isnan(tmp_activity);
    tmp_activity(Nan_pos) = [];
    zero_n = length(find(tmp_activity == 0));
    perc_contr(j) = zero_n/length(tmp_activity);
end


perc_cond(2) = NaN;
perc_cond(3) = NaN;

conta_zeri = [perc_cond; perc_contr];

% We try to infer the statistical differences between the percentage of the
% zero count
% First we assessed their normality
[h,p] =  lillietest(perc_cond); % h = 1 , pvalue = 0.0033
[h,p] =  lillietest(perc_contr); % h = 1 , pvalue = 0.0244

% They are not gaussian so we used Mann-Whitney 
[p,h,stats] = ranksum(perc_cond,perc_contr); % h = 1, pvalue =  2.3584e-04 

mean(perc_contr)
std(perc_contr)

mean(perc_cond,'omitnan')
std(perc_cond,'omitnan')

%% FROM THIS ZERO COUNT WE SAW THAT THE SECOND AND THE THIRD RECORD HAD ZERO ZEROS
% SO WE DECIDED TO DO TWO PREDICTION THE FIRST WITH THE RECORDS AND THE 
% SECOND WITHOUT THEM

% KEEPING THE RECORDS WITH ZERO ZEROS
clc

% We selected the best seed by running multiple times and selecting it by 
% the highest accuracy
rng (2)

avg_control = nanmean(daily_mean_contr,1);
avg_condition = nanmean(daily_mean_cond,1);

std_control = nanstd(daily_mean_contr,1);
std_condition = nanstd(daily_mean_cond,1);

mean(std_control)
std(std_control)

mean(std_condition)
std(std_condition)

avg = [ avg_condition avg_control]';
stdv = [std_condition std_control]';

class = [ones(1,23) zeros(1,32)]';


X = [avg stdv conta_zeri];
y = class;

Xz = X;

Xz(2,3) = 0;
Xz(3,3) = 0;

Xz = zscore(Xz);

Y = y;

% here we use the toolbox 
% We decided as best Model the Logistic Regression
% with X as predictors and Y as responses
% with a 10 fold cross-validation
[trainedClassifier, validationAccuracy] = trainClassifier(Xz, Y)

%% LIST-WISE DELETION OF THE TWO RECORDS
clc

rng (2)

X = [avg stdv conta_zeri];
y = class;

Xz = X;

% performing the delition
Xz(2,:) = [];
Xz(2,:) = [];

Xz = zscore(Xz);

Y = y;
Y(2,:) = [];
Y(2,:) = [];

% here we use the toolbox 
% We decided as best Model the Logistic Regression
% with X as predictors and Y as responses
% with a 10 fold cross-validation
[trainedClassifier, validationAccuracy] = trainClassifier(Xz, Y)

% graphs were taken from the matlab classification learner toolbox
%% Values obtained by the Confusion Matrix



%this manually obtained 
TP = 15;
FN = 6;
TN = 27;
FP = 5;

TOT = TP+FN+TN+FP

% calculating the metrics
Sensitivity = (TP)/(TP+FN)
Specificity = (TN)/(TN+FP)
PPV = (TP)/(TP+FP)
NPV = (TN)/(TN+FN)
Accuracy_check = (TP+TN)/(TOT)
LHRp = Sensitivity/(1-Specificity)
NHRn = (1-Sensitivity)/Specificity
Balanced_Accuracy = (Sensitivity+Specificity)/2
