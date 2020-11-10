function [peaks_ind,peaks_amp]=QRS_FP_FN(ECG,fs,subjectid,subjectdir)
% The detection of QRS complexes was performed using a modified version of the Pan-Tompkins algorithm (Pan and Tompkins, 1985),
% which consists of submitting the annotations resulting from the original implementation to the following additional 
% routines. Additional steps to correct for both false negative and positive QRS detection were applied. 
% False negatives/positives were first detected by probing the time differences between two consecutive
% annotations for values higher/lower than 1.5/0.4 times the median time difference across all annotations. 
% In the case of a false negative, an R peak was manually added in the middle of the two adjacent R peak occurrences; 
% regarding the false positives, the annotation yielding the lowest ECG amplitude between the two adjacent occurrences 
% was deemed spurious and removed. Next, the corrected annotations were subjected to an alignment step, 
% by shifting them to match the occurrence of the R peak. Finally, the modified algorithm was applied to the 
% time-inverted ECG signal, performing a back-search for the R peaks, in order to overcome the poor detection
% at both ends of the ECG signal. The final set is obtained by taking the union of the two sets of detected 
% R peaks. An almost flawless QRS detector is of the utmost importance, as some of the ECG features that are
% modulated by the respiratory cycle, and thus used to compute a surrogate of the respiratory signal,
% are intimately related with changes within each heartbeat.

%% Author: unknown
% No information was found on who initially wrote this code. 

%% Modified: Reyhaneh Bakhtiari
% Modifications: Visually inspecting the detected QRS and possibility to add or remove R peak is added to the code.
% University of Alberta
% email : reyhaneh@gmail.com

%%
% clear FN FP ind_FN amp_FN min_i_FP min_FP ind_FP
% clear qrs_i_FP_mod qrs_amp_FP_mod qrs_i_FP_FN_mod tmp_qrs_amp_FP_FN ind_FP_sort qrs_amp_FP_FN_mod
% clear qrs_amp_raw qrs_i_raw
% fs=250; %sampling frequency
max_hr=100; % Thr for high Heart rate
min_hr=60; % Thr for low Heart rate
% avg_hr=(max_hr+min_hr)/2;
% avg_interval=60*fs/avg_hr;
avg_interval=250; % use this thr if detected HR is not in range
for i_run=1:2
    if i_run==1;
        ecg(:,i_run)=double(ECG);
        color_code='c*';
    else
        ecg(:,i_run)=double(ECG(end:-1:1));
        color_code='b+';
    end
    [qrs_amp_raw{i_run},qrs_i_raw{i_run},delay,ecg_m,ecg_h]=pan_tompkin(ecg(:,i_run),fs,1,subjectid);
    med_qrs_i_raw(i_run,:)=median(diff(qrs_i_raw{i_run})); % median of time diff across all peaks
    
    med_beat(i_run)=60*fs/med_qrs_i_raw(i_run,:);
    if med_beat(i_run)<=max_hr && med_beat(i_run)>=min_hr
        disp(['Median heart beat rate: ',num2str(med_beat(i_run)),' beats/min, med_qrs_r_raw is to: ', num2str(med_qrs_i_raw(i_run,:))])
        FN_Thr(i_run,:)=med_qrs_i_raw(i_run,:)*1.5;
        FP_Thr(i_run,:)=med_qrs_i_raw(i_run,:)*.4;
        
    else
        disp(['@@@ Warning!!! Atypical heart beat rate: ',num2str(med_beat(i_run)),' beat/min',' beat/min, med_qrs_r_raw is set to: ', num2str(avg_interval)])
        FN_Thr(i_run,:)=avg_interval*1.5;
        FP_Thr(i_run,:)=avg_interval*.4;
    end
    
    % find false negative and false positive peaks based on the threshold
    FN{i_run}=find(diff(qrs_i_raw{i_run})>FN_Thr(i_run,:));
    FP{i_run}=find(diff(qrs_i_raw{i_run})<FP_Thr(i_run,:));
   
    % FN adding
    % manually add R_peak for missing False negative peaks in the middle of
    % two adjacent R_peak occurness
    ind_FN_notaligned{i_run}=round((qrs_i_raw{i_run}(FN{i_run}+1)+qrs_i_raw{i_run}(FN{i_run}))/2);
    amp_FN_notaligned{i_run}=((qrs_amp_raw{i_run}(FN{i_run}+1)+qrs_amp_raw{i_run}(FN{i_run}))/2);
    %% This needs work!!!
    % The corrected annotations should be subjected to an alignment step by shifting them to match
    % the occurence of the R_peak
    %%%
%     tic
    ind_FN_tmp=[];
    amp_FN_tmp=[];
    for i_FN=1:length(ind_FN_notaligned{i_run})
        % find a Rpeak in the moving averaged filtered signal, in interval
        % +/-100 datapoints around the mannually added points, maximum peak
        % in +/-50 will be selected as a R_peak if available otherwise the
        % mean value will be used.
        my_lim=[ind_FN_notaligned{i_run}(i_FN)-100:ind_FN_notaligned{i_run}(i_FN)+100];
        [pks,locs] = findpeaks(ecg_m(my_lim),'MINPEAKDISTANCE',round(0.2*fs));
        ind_close_peaks=find(locs>50 & locs<150);
        if size(ind_close_peaks)>0
            [rpeak,ind_rpeak]=max(pks(ind_close_peaks));
            my_Rpeak=my_lim(locs(ind_close_peaks(ind_rpeak)));
        else
            my_Rpeak=ind_FN_notaligned{i_run}(i_FN);
        end
        %% locate the corresponding peak in the filtered signal
        if my_Rpeak-round(0.150*fs)>= 1 && my_Rpeak<= length(ecg_h)
            [y_i,x_i] = max(ecg_h(my_Rpeak-round(0.150*fs):my_Rpeak));
        else
            if my_Rpeak-round(0.150*fs)< 1
                [y_i,x_i] = max(ecg_h(1:my_Rpeak));
                ser_back = 1;
            elseif my_Rpeak>= length(ecg_h)
                [y_i,x_i] = max(ecg_h(my_Rpeak-round(0.150*fs):end));
            end
        end
        ind_FN_tmp(i_FN)=my_Rpeak-round(0.150*fs)+x_i;
        amp_FN_tmp(i_FN)=y_i;
        
    end
    ind_FN{i_run}=ind_FN_tmp;
    amp_FN{i_run}=amp_FN_tmp;
    
    % FP removing
    % annotation yielding the lowest ECG amp between two adjacent occureness was deemed spurius and removed  
    for i_FP=1:length(FP{i_run})
        [min_FP,min_i_FP]=min([qrs_amp_raw{i_run}(FP{i_run}(i_FP)),qrs_amp_raw{i_run}(FP{i_run}(i_FP)+1)]);
        ind_FP{i_run}(i_FP)=FP{i_run}(i_FP)+min_i_FP-1;
    end
    if length(FP{i_run})==0
        ind_FP{i_run}(i_FP)=[];
    end
    
    qrs_i_FP_mod{i_run}=qrs_i_raw{i_run};
    qrs_amp_FP_mod{i_run}=qrs_amp_raw{i_run};
    qrs_i_FP_mod{i_run}(ind_FP{i_run})=[];
    qrs_amp_FP_mod{i_run}(ind_FP{i_run})=[];
    
    % FN & FP adjustment
    [qrs_i_FP_FN_mod{i_run},ind_FP_sort{i_run}]=sort([qrs_i_FP_mod{i_run},ind_FN{i_run}]);
    tmp_qrs_amp_FP_FN{i_run}=[qrs_amp_FP_mod{i_run},amp_FN{i_run}];
    qrs_amp_FP_FN_mod{i_run}=tmp_qrs_amp_FP_FN{i_run}(ind_FP_sort{i_run});
    
    my_fig{i_run}=gcf;
    az(i_run,1)=subplot(311); az(i_run,2)=subplot(312);
    hold on,scatter(qrs_i_FP_FN_mod{i_run},qrs_amp_FP_FN_mod{i_run},color_code);
    az(i_run,3)=subplot(313);
    hold on,scatter(qrs_i_FP_FN_mod{i_run},10*qrs_amp_FP_FN_mod{i_run},color_code);
    linkaxes(az(i_run,:),'x');
    zoom on;
    disp(['  Number of  detected R-peaks in run ', num2str(i_run),': ',num2str(length(tmp_qrs_amp_FP_FN{i_run})), ' (initially: ',num2str(length(qrs_i_raw{i_run})),')'])

end

%% taking the union of two sets of R-peaks
uninverted_qrs_i_FP_FN_mod=size(ecg,1)+1-qrs_i_FP_FN_mod{2}(end:-1:1);
uninverted_qrs_amp_FP_FN_mod=qrs_amp_FP_FN_mod{2}(end:-1:1);

if exist('my_fig{1}.Number')
    figure(my_fig{1}.Number)
else
    figure(my_fig{1})
end

ax(1) = subplot(311); ax(2) = subplot(312);
hold on,scatter(uninverted_qrs_i_FP_FN_mod,uninverted_qrs_amp_FP_FN_mod,color_code);
ax(3)= subplot(313);
hold on,scatter(uninverted_qrs_i_FP_FN_mod,uninverted_qrs_amp_FP_FN_mod,color_code);
linkaxes(ax,'x');

[union_qrs_i,ia,ib] = union(qrs_i_FP_FN_mod{1},uninverted_qrs_i_FP_FN_mod);
[union_tmp,Ind_tmp]=sort([qrs_i_FP_FN_mod{1}(ia),uninverted_qrs_i_FP_FN_mod(ib)]);
cat_qrs_amp=[qrs_amp_FP_FN_mod{1}(ia),uninverted_qrs_amp_FP_FN_mod(ib)];
union_qrs_amp=cat_qrs_amp(Ind_tmp);
% find the peaks very close due to numerical calculation
close_peaks=find(abs(diff(union_qrs_i))<2);
if length(close_peaks)>0
    disp('Near peaks (due to calculation) are found at the following positions, the top row is selected')
    disp(num2str([union_qrs_i(close_peaks);union_qrs_i(close_peaks+1)]))
    union_qrs_i(close_peaks+1)=[];
    union_qrs_amp(close_peaks+1)=[];
end

if exist('my_fig{1}.Number')
    figure(my_fig{1}.Number)
else
    figure(my_fig{1})
end
ax(1) = subplot(311);
hold on;scatter(union_qrs_i,union_qrs_amp+.05,'filled','kd');
ax(2) = subplot(312);
hold on; scatter(union_qrs_i,union_qrs_amp/10,'filled','kd');
ax(3) = subplot(313);
hold on; scatter(union_qrs_i,union_qrs_amp+.05,'filled','kd');
linkaxes(ax,'x')
%% Final step: manual screening removing of peaks
last_FP=find(diff(union_qrs_i)<FP_Thr(1)*1.1);
if length(last_FP)>0
    if last_FP(1)==1
        last_FP(1)=[];
    end
    if last_FP(end)==length(union_qrs_i)
        last_FP(end)=[];
    end
end
disp(['There are ',num2str(length(last_FP)),' close points that should be checked further']);

%% Automatic FP rejection phase 1:
% a single isolated FP is found (two close points). Between the two close points the one that is furthest to the middle point is rejected 
FP_automatic_rej_PH1=[];
FP_remain_PH1=[];
exception=0;
for i_CFP=1:length(last_FP)
    if i_CFP<length(last_FP) & i_CFP>1
        if last_FP(i_CFP)+1==last_FP(i_CFP+1) || last_FP(i_CFP)-1==last_FP(i_CFP-1)
%             disp(['check these points: ', num2str(union_qrs_i(last_FP(i_CFP))), '(',num2str(last_FP(i_CFP)), '), ', num2str(union_qrs_i(last_FP(i_CFP)+1)),'(',num2str(last_FP(i_CFP)+1),')']);
            FP_remain_PH1=[FP_remain_PH1 last_FP(i_CFP)];
            exception=1;
        end
    elseif i_CFP==1 & length(last_FP)>1 %Reyhaneh
        if last_FP(i_CFP)+1==last_FP(i_CFP+1)
%             disp(['check these points: ', num2str(union_qrs_i(last_FP(i_CFP))), '(',num2str(last_FP(i_CFP)), '), ', num2str(union_qrs_i(last_FP(i_CFP)+1)),'(',num2str(last_FP(i_CFP)+1),')']);
            FP_remain_PH1=[FP_remain_PH1 last_FP(i_CFP)];
            exception=1;
        end
    elseif i_CFP==length(last_FP) & length(last_FP)>1 %Reyhaneh
        if last_FP(i_CFP)-1==last_FP(i_CFP-1)
%             disp(['check these points: ', num2str(union_qrs_i(last_FP(i_CFP))), '(',num2str(last_FP(i_CFP)), '), ', num2str(union_qrs_i(last_FP(i_CFP)+1)),'(',num2str(last_FP(i_CFP)+1),')']);
            FP_remain_PH1=[FP_remain_PH1 last_FP(i_CFP)];
            exception=1;
        end
    end
    if  exception==0
        %     disp(['check these points: ', num2str(union_qrs_i(last_FP(i_CFP))), '(',num2str(last_FP(i_CFP)), '), ', num2str(union_qrs_i(last_FP(i_CFP)+1)),'(',num2str(last_FP(i_CFP)+1),')']);
        mean_closeness_i=abs((union_qrs_i(last_FP(i_CFP))-union_qrs_i(last_FP(i_CFP)-1))- (union_qrs_i(last_FP(i_CFP)+2)-union_qrs_i(last_FP(i_CFP))));
        mean_closeness_iP1=abs((union_qrs_i(last_FP(i_CFP)+1)-union_qrs_i(last_FP(i_CFP)-1))- (union_qrs_i(last_FP(i_CFP)+2)-union_qrs_i(last_FP(i_CFP)+1)));
        [max_val,max_ind]=max([mean_closeness_i,mean_closeness_iP1]);
        FP_automatic_rej_PH1=[FP_automatic_rej_PH1, (last_FP(i_CFP)+max_ind-1)];
    end
    exception=0;
end

%% Automatic FP rejection phase 2:
% only 2 consequent FP are found, (three close points) the middle one is kept, the other two are rejected
FP_automatic_rej_PH2=[];
FP_remain_PH2=[];
exception=0;
skip=0;
jump_node=find(diff(FP_remain_PH1)~=1);
diff_jump_node=diff([0,jump_node]);
for i_CFP=1:length(diff_jump_node)
    if diff_jump_node(i_CFP)==2
        FP_automatic_rej_PH2=[FP_automatic_rej_PH2, FP_remain_PH1(jump_node(i_CFP))-1, FP_remain_PH1(jump_node(i_CFP))+1];
    else
        FP_remain_PH2=[FP_remain_PH2 FP_remain_PH1(jump_node(i_CFP))-(diff_jump_node(i_CFP)-1:-1:0)];
    end
end

%% Manual FP rejection:
FP_manual_rej=[];
jump_node_manual=[0,find(diff([FP_remain_PH2,inf])~=1)];

for i_CFP=2:1
% for i_CFP=2:length([jump_node_manual])
%     satisfied='N';
%     while ~strcmpi(satisfied,'Y')
        disp(['check these points: ', ...
            num2str(union_qrs_i(FP_remain_PH2(jump_node_manual(i_CFP-1)+1))),' : ',...
            num2str(union_qrs_i(FP_remain_PH2(jump_node_manual(i_CFP))+1)),' ('...
            num2str([FP_remain_PH2(jump_node_manual(i_CFP-1)+1):FP_remain_PH2(jump_node_manual(i_CFP))+1]),')'])
        if exist('my_fig{1}.Number')
            figure(my_fig{1}.Number)
        else
            figure(my_fig{1})
        end
        az(i_run,1)=subplot(311);az(i_run,2)=subplot(312);az(i_run,3)=subplot(313);
        linkaxes(az(i_run,:),'x');
        xlim([union_qrs_i(FP_remain_PH2(jump_node_manual(i_CFP-1)+1)-3),union_qrs_i(FP_remain_PH2(jump_node_manual(i_CFP))+1+3)])
        FP_points_tmp=input('index of the point to be removed (value in parentes) ');
%         satisfied=input(['Happy with selected point at ',num2str(union_qrs_i(FP_points)),'(',num2str(FP_points),')? Y/N'],'s');
%     end
    FP_manual_rej=[FP_manual_rej, FP_points_tmp];
end

%%
union_qrs_i_final=union_qrs_i;
union_qrs_amp_final=union_qrs_amp;

union_qrs_i_final(sort([FP_automatic_rej_PH1, FP_automatic_rej_PH2,FP_manual_rej]))=[];
union_qrs_amp_final(sort([FP_automatic_rej_PH1, FP_automatic_rej_PH2,FP_manual_rej]))=[];

if exist('my_fig{1}.Number')
    figure(my_fig{1}.Number)
else
    figure(my_fig{1})
end

subplot(311);
hold on;scatter(union_qrs_i_final,union_qrs_amp_final+.05,'filled','mo');
subplot(312);
hold on;scatter(union_qrs_i_final,union_qrs_amp_final/10,'filled','mo');
subplot(313);
hold on;scatter(union_qrs_i_final,union_qrs_amp_final+.05,'filled','mo');

%% Final visual inspection
xlim('auto')
added_vec=1:length(union_qrs_i_final);
subplot(313);
line(repmat(union_qrs_i_final,[2 1]),repmat([-5000; 5000],size(union_qrs_i_final))+[-added_vec;added_vec],'LineWidth',1,'LineStyle','-','Color','g');
% subplot(312);
% line(repmat(union_qrs_i_final,[2 1]),repmat([-0000; 0000],size(union_qrs_i_final))+0.0001*[-added_vec;added_vec],'LineWidth',1,'LineStyle','-','Color','g');
disp('********************************************')
disp(['Visual Inspection to Remove Incorrect R-peaks, in Figure(', num2str(my_fig{1}),'), the 3rd subplot, green train'])
Final_screening=sort(input('Enter index of other points to be removed:  '));
union_qrs_i_final(Final_screening)=[];
union_qrs_amp_final(Final_screening)=[];

disp('********************************************')
disp(['Visual Inspection to add Incorrectly removed R-peaks, in Figure(', num2str(my_fig{1}),'), the 1st subplot'])
Final_screening_added=input('Enter time stamp and amplitude of qrs point (eg. [325,.53; ...] ');
if ~isempty(Final_screening_added)
    union_qrs_i_final=[union_qrs_i_final,Final_screening_added(:,1)'];
    union_qrs_amp_final=[union_qrs_amp_final,Final_screening_added(:,2)'];
end
added_vec=1:length(union_qrs_i_final);
subplot(313);
line(repmat(union_qrs_i_final,[2 1]),repmat([-5000; 5000],size(union_qrs_i_final))+[-added_vec;added_vec],'LineWidth',1,'LineStyle','--','Color','r');

satisfied=input(' Happy with Rpeaks? If not, terminate program, and re-run it by pressing Ctrl^C ','s');


% save R-peak
Pulse_width_THR=avg_interval;
peaks_ind=union_qrs_i_final;
peaks_amp=union_qrs_amp_final;
save([subjectdir,'R_peaks.mat'],'peaks_ind','peaks_amp','Pulse_width_THR','Final_screening');
disp(['  Number of finaly detected R-peaks: ',num2str(length(peaks_ind))])
disp(['Median R-peaks interval: ',num2str(60*fs/median(diff(peaks_ind))),'(',num2str(median(diff(peaks_ind))),')',' beats/min(bins)'])

