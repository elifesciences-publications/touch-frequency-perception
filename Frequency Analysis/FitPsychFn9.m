function  [DATA,bandNo,mydata] = FitPsychFn9(exptdir,rawdata,Opts,stimParams,paramcols)

% exptdir='C:\Google Drive\NeuRA-shakers\Dynamic Load Playwave\Experiments\2013-06-05\30Hz_0Hz_0Ph_0Ph_Max50Hz_Min25Hz_alternatingshakers\Data\';
% fname='2013-06-05_13-17-06_sarah_forearm_close';

%%%%%%% input variables %%%%%
%
%resptype = '12'; % responses were 'first' or 'second' in the interval
%resptype = 'LR'; % responses were 'left' or 'right'
labelcols = [5 6]; %columns of the data file that contain stimulus labels (1,2,3,..)
% freqcols = 1:4; %columns of the data file that contain stimulus frequencies
respcol = 7; %column of the data file that contains the responses
%
%%%%%%%%%%%%%%

% datadir = strcat(exptdir,'Data',filesep);
id = strfind(exptdir,filesep);
DATA.figdir = strcat(exptdir(1:id(end-1)),'Figures',filesep);
mkdir(DATA.figdir);

% check for missing responses, remove them and print a warning
missing = sum(rawdata(:,respcol)==0);
if missing > 0
    warning('%d missing responses excluded from analysis.',missing);
    rawdata = rawdata(rawdata(:,respcol)~=0,:);
end

% if strcmp(stands,'separate')
%     stimParams = stimParams(stimParams(:,12)==1,:);
% end
stimlab = rawdata(:,labelcols);
% stimfreqRaw = rawdata(:,freqcols);
% stimfreq = [max(stimfreqRaw(:,1),stimfreqRaw(:,2)) max(stimfreqRaw(:,3),stimfreqRaw(:,4))];

freq_t = stimParams(:,1:2);
amp_t = stimParams(:,3:4);
amp_t(freq_t==0) = 0;
freq_t(amp_t==0) = 0;
stimParams(:,1:2) = freq_t;
stimParams(:,3:4) = amp_t;
amp_t = []; freq_t = [];
if paramcols == 1 || paramcols == 2    
    t = stimParams(:,1:2);
    for i = 1:numel(stimParams(:,1))
        if sum(stimParams(i,1:2)==0) == 1      
            freq_t(i,1) = t(i,t(i,:)~=0);
        else
            freq_t(i,1) =  t(i,paramcols);
        end 
    end
    stimParams(:,paramcols) = freq_t(:,1);
elseif paramcols == 3 || paramcols == 4   
    t = stimParams(:,3:4);
    for i = 1:numel(stimParams(:,1))
        if sum(stimParams(i,3:4)==0) == 1
            amp_t(i,1) = t(i,t(i,:)~=0);
        else
            amp_t(i,1) =  t(i,(paramcols-2));
        end
    end
    stimParams(:,paramcols) = amp_t(:,1);
end



stimfreq = zeros(numel(stimlab(:,1)),2);
for i = 1:numel(stimlab(:,1))
    for j = 1:2
        stimfreq(i,j) = stimParams(stimlab(i,j),paramcols);
    end
end

bandlab = stimParams(:,12);
if strcmp(Opts.Stands,'combine')
    stndlabel = find(stimParams(:,11)); %label for the standard stimulus
    for i = 2:numel(stndlabel)
        stimlab(stimlab==stndlabel(i)) = stndlabel(1);
    end
    stndlabel = stndlabel(1);
    bandNo = 1;     
else
    stndlabel = find(stimParams(:,11)); %label for the standard stimulus
    bandNo = numel(unique(bandlab));  
end



for j = 1:bandNo  
    %determine which is the standard stimulus and what value of
    %the comparison stimulus was used  
    ids = logical(sum(stimlab==stndlabel(j),2));
    rawdata_temp = rawdata(ids,:);
    stimlab_temp = stimlab(ids,:);
    stimfreq_temp = stimfreq(ids,:);
    stndi = zeros(length(rawdata_temp),1);
    stndf = zeros(length(rawdata_temp),1);
    compf = zeros(length(rawdata_temp),1);
    for r=1:size(rawdata_temp,1)
        stndmatch = find(stimlab_temp(r,:)==stndlabel(j));
        %if comparison = standard, alternate which one is the standard
        if length(stndmatch) > 1
            stndi(r) = mod(r,2)+1;
        elseif length(stndmatch)==1
            stndi(r) = stndmatch;
        end
        stndf(r) = stimfreq_temp(r,stndi(r));
        %comparison is the one that isn't the standard
        compf(r) = stimfreq_temp(r,mod(stndi(r),2)+1);
    end
    
    % code response vector where:
    % 0 means standard was faster
    % 1 means comparison was faster
    if strcmp(Opts.resptype,'LR')
        resp = rawdata_temp(:,respcol) == 2;
    elseif strcmp(Opts.resptype,'12')
        resp = rawdata_temp(:,respcol) ~= stndi;
    end
    resp = resp(rawdata_temp(:,respcol) ~= 0);

    %count responses by comparison frequency
    DATA.compfreq{j} = unique(compf);
    nCompFaster = zeros(length(DATA.compfreq{j}),1);
    nTrials = zeros(length(DATA.compfreq{j}),1);
    for f = 1:length(DATA.compfreq{j})
        thisFreq = DATA.compfreq{j}(f);
        compResp = resp(compf==thisFreq);
        nCompFaster(f) = sum(compResp);
        nTrials(f) = length(compResp);
    end
    
    if sum(nCompFaster) == 0
        error('Aborting. Comparison was NEVER perceived faster')
    elseif sum(nTrials - nCompFaster) == 0
        error('Aborting. Comparison was ALWAYS perceived faster')
    end

    pCompFaster = nCompFaster./nTrials;
    DATA.pCompFaster_adj{j} = pCompFaster;

    %adjust for floor and ceiling responses
    compfreqC = DATA.compfreq{j} - mean(DATA.compfreq{j});
    if corr(DATA.compfreq{j},pCompFaster) < 0
        compfreqC = -compfreqC;
    end

    fl = pCompFaster == 0;
    while sum(fl) > 1
        DATA.pCompFaster_adj{j}( compfreqC.*fl == min(compfreqC.*fl) ) = NaN;
        fl( compfreqC.*fl == min(compfreqC.*fl) ) = 0;
    end

    cl = pCompFaster == 1;
    while sum(cl) > 1
        DATA.pCompFaster_adj{j}( compfreqC.*cl == max(compfreqC.*cl) ) = NaN;
        cl( compfreqC.*cl == max(compfreqC.*cl) ) = 0;
    end

    DATA.pCompFaster_adj{j}(fl) = 1./(nTrials(fl).*2);
    DATA.pCompFaster_adj{j}(cl) = (nTrials(cl).*2-1)./(nTrials(cl).*2);

    %logit transform the data
    logitTransCompFaster = log(DATA.pCompFaster_adj{j}./(1-DATA.pCompFaster_adj{j}));

    %data
    mydata{j} = [DATA.compfreq{j} nCompFaster nTrials pCompFaster DATA.pCompFaster_adj{j} ...
        logitTransCompFaster];
    disp([DATA.compfreq{j} nCompFaster nTrials pCompFaster DATA.pCompFaster_adj{j} ...
        logitTransCompFaster])

    %fit linear regression, calculate PSE & slope
    linfit = polyfit(DATA.compfreq{j}(isnan(logitTransCompFaster)==0),...
        logitTransCompFaster(isnan(logitTransCompFaster)==0),1);
    DATA.PSE(j) = -linfit(2)/linfit(1)
    slope = linfit(1);

    thresh25 = (log(0.25/(1-0.25))-linfit(2)) / linfit(1);
    thresh75 = (log(0.75/(1-0.75))-linfit(2)) / linfit(1);

    DATA.WeberLower(j) = (DATA.PSE(j)-thresh25)/DATA.PSE(j)
    DATA.WeberHigher(j) = (thresh75-DATA.PSE(j))/DATA.PSE(j)

    %points to make the line for the plot
    DATA.freqFine{j} = min(DATA.compfreq{j}) ...
        : (max(DATA.compfreq{j})-min(DATA.compfreq{j}))/100 : max(DATA.compfreq{j});
    logitTransFine = polyval(linfit,DATA.freqFine{j});
    DATA.curveFine{j} = exp(logitTransFine)./(1+exp(logitTransFine));
    DATA.stndf(j) = stndf(1);
    
%     h = figure('position', [500 500 800 500]); % create new figure 
    
    %plot logit transformed data
%     subplot(1,2,1)
    
% %     plot(DATA.compfreq{j},DATA.pCompFaster_adj{j},'k.','markersize',10);
% 
%     set(gca, 'fontsize',12);
%     axis([min(DATA.compfreq{j})-std(DATA.compfreq{j})/3 max(DATA.compfreq{j})+std(DATA.compfreq{j})/3 -0.1 1.1])
%     %axis([min(DATA.compfreq{j}) max(DATA.compfreq{j}) 0 1]);
%     hold on;
%     %plot fitted function
%     plot(DATA.freqFine{j}{j},DATA.curveFine{j},'r-','linewidth',1);
% 
%     xlabel('Comparison Frequency (Hz)');
%     ylabel('Proportion Comparison Faster');
%     title(strcat('Standard frequency = ',num2str(stndf(1)),...
%         'Hz; PSE = ',num2str(PSE),'; WeberLower = ',num2str(WeberLower),...
%         '; WeberHigher = ',num2str(WeberHigher)));
% 
% %     subplot(1,2,2)
% %     axis off
% %     text(-0.1,0.5,num2str(mydata))
% 
%     %save plot to .tiff file in /Figures folder
% %     print(h, '-dtiffnocompression', strcat(figdir,fname,'.tif'));
%     hold off;

end
