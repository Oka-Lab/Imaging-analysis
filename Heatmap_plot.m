
order = [1 2 3];
% order ; to look at all cells responding to a particular stimulus

stim_times = [ 
              56 80; ...
              56 80; ...
              56 80; ... 
              ];
          
folderlist = {
              };
suffix = '_extractedData.mat'; % '_extractedData_withSpikes.mat' or '_extractedData.mat'

threshold = 2.5; %threshold for being counted as an active cell

F_scaled_allmice = [];
responses_allmice = [];
stim_preferences = []; % counts % of cells that prefer a stimulus, for each mouse
stim_responsive = []; % counts the NUMBER of cells that responded to a stimulus, for each mouse
mean_responses = []; % computes the mean population activity (relative to baseline) for each stim/mouse

for mouse = 1:length(folderlist)
    fileList = dir([folderlist{mouse} '*' suffix]);
    
    F_scaled_all = [];
    responses = [];
    stimName = {};
    stimDuration = [];
    
    count=0;
    for i=order
        count=count+1;       
        load([fileList(i).folder filesep fileList(i).name]);
                
        use_signal = F_dff_saved; % F_full_saved or F_dff_saved
        F_scaled = use_signal;
        
        baseline = mean(use_signal(:,1:20),2);
        F_scaled = use_signal - baseline;
        F_scaled = bsxfun(@times,F_scaled,1./std(use_signal(:,1:20),[],2));

        F_scaled_all = [F_scaled_all F_scaled];
        responses(:,count) = mean(F_scaled(:,stim_times(count,1):stim_times(count,2)),2);
        [~,stimName{count}] = fileparts(fileList(i).name);
        stimName{count} = strrep(stimName{count},['_nr' suffix(1:end-4)],'');
        stimDuration(count) = size(F_scaled,2);
    end
    
    F_scaled_allmice = [F_scaled_allmice; F_scaled_all];
    responses_allmice = [responses_allmice; responses];
    
    [peak_response,preferred_stim] = max(responses,[],2);
    active_cells = peak_response>=threshold;
    
    stim_preferences(mouse,:) = hist(preferred_stim(active_cells),1:max(preferred_stim));
    stim_preferences(mouse,:) = stim_preferences(mouse,:)/sum(stim_preferences(mouse,:))*100;
    
    for i = 1:size(responses,2)
        stim_responsive(mouse,i) = sum(responses(:,i)>=threshold);
        mean_responses(mouse,i) = mean(responses(:,i));
    end
end


[peak_response,preferred_stim] = max(responses_allmice,[],2);
active_cells = peak_response>=threshold;

stimName = strrep(stimName','_',' ');


F_sorted = [];
for i= 1:size(responses_allmice,2)
    F_temp = [responses_allmice(preferred_stim==i & active_cells,i) ...
              F_scaled_allmice(preferred_stim==i & active_cells,:)];
    F_temp = flipud(sortrows(F_temp,1));
    F_sorted = [F_sorted; F_temp(:,2:end)];
end

figure(1);clf;
imagesc(F_sorted);
hold on;
for i=1:size(responses,2)
    plot(sum(stimDuration(1:i))*[1 1]+.5,get(gca,'ylim'),'w');
     plot(sum(stimDuration(1:i-1))*[1 1] + stim_times(i,1),get(gca,'ylim'),'r');
     plot(sum(stimDuration(1:i-1))*[1 1] + stim_times(i,2),get(gca,'ylim'),'r');
   % plot(get(gca,'xlim'),sum(preferred_stim<=i & active_cells)*[1 1]+.5,'w');
end
%set(gca,'xtick',cumsum([0 stimDuration(1:end-1)]) + stimDuration/2,...
   %'xticklabel',stimName,'xticklabelrotation',45);
colorbar;
box off;


figure(2);clf;
stim1 = 1;
stim2 = 2;
stim3 = 3;

stim_A = sum(responses_allmice(:,stim1)>=threshold);
stim_A_only = sum(responses_allmice(:,stim1)>=threshold & responses_allmice(:,stim2)<threshold);
stim_B = sum(responses_allmice(:,stim2)>=threshold);
stim_B_only = sum(responses_allmice(:,stim1)<threshold & responses_allmice(:,stim2)>=threshold);
stim_C = sum(responses_allmice(:,stim3)>=threshold);
stim_D = stim_A + stim_B + stim_C;
%stim_F = sum(responses_allmice(:,stim4)>=threshold);
both_stimsAB = sum(responses_allmice(:,stim1)>=threshold & responses_allmice(:,stim2)>=threshold);
both_stimsAC = sum(responses_allmice(:,stim1)>=threshold & responses_allmice(:,stim3)>=threshold);
both_stimsBC = sum(responses_allmice(:,stim2)>=threshold & responses_allmice(:,stim3)>=threshold);
both_stimsABC = sum(responses_allmice(:,stim1)>=threshold & responses_allmice(:,stim2)>=threshold & responses_allmice(:,stim3)>=threshold);
neither_stim = sum(responses_allmice(:,stim1)<threshold & responses_allmice(:,stim2)<threshold);


% need to download venn first:

 venn([stim_A stim_B stim_C], [both_stimsAB both_stimsAC both_stimsBC both_stimsABC], 'FaceColor',{'r','g','c'},'FaceAlpha',{1,0.6,0.6},'EdgeColor','black');
 legend({[stimName{stim1}], [stimName{stim2}], [stimName{stim3}]});

 box off;
 axis off;
axis equal;