%% This version of code was used to analyze OXT neuron activity and behavior in response to chemosensory cues
%% SECTION 1 is for aggregating raw data into structured array
% SKIP TO SECTION 2 if dataarray is already generated.
% email me for dataarray to test code (too large to upload)

clf; close all; clear all;

celltype = sprintf('OXTPO');
experiment = 'all';  %allcues
condition = 'water';
nchannel = 1;
minframes = 2800; 
%short = 1hr or less
%medium = 1hr or more
%long = 3hr or more

%specify here which pump channels would correspond to which odors
controlodor = 'water';
firstodor = 'kin';
secondodor = 'nonkin';
thirdodor = 'adult';
fourthodor = 'para';

if strcmp(condition, 'water')
    pump = 1;
    stim = 0;
elseif strcmp(condition, 'optovin')
    pump = 1;
    stim = 1;
end

disk = 15;


pathname =  sprintf('%s%d%s/%s/%s/', '/Volumes/CW-', disk, '/2P_socialodor/oxyGC6s/aggregated/', condition, experiment)
cd(sprintf('%s/%s/%s/', '/Users/weel/CW Lab Dropbox/Caroline Wee/Carolinebackup/OXTai/epsfiles/liquidlens_oxyGC6s/', celltype, 'social',experiment));

files = dir(pathname);

FISH = 0;

conditionarray = {};
conditionarray{1} = 'water'
conditionarray{2} = 'optovin'
%% AGGREGATE DATA

% SET UP DATA ARRAY
catstruct = struct('info',[],  'stim', [], 'bout', [], 'calcium', [], 'tail', []);
dataarray = struct('water', catstruct, 'optovin', catstruct);

for c = 1:length(conditionarray)
    
    condition_tmp = char(conditionarray{c});
    
    %stimtriggered (account for either pump or teensy or both stims)
    stimstruct = struct('lastpointangle', [], 'cumangle', [], 'power', [], 'calciumpeaks', [], 'calciumsum', [], 'intensities', [], 'location', []);
    dataarray.(condition_tmp).stim = struct('pump', stimstruct, 'teensy', stimstruct);
    
    stimstruct = struct('control', [], 'odor1', [], 'odor2', [], 'odor3', [], 'odor4', [], 'fishIDcontrol', [], 'fishIDodor1', [], 'fishIDodor2', [], 'fishIDodor3', [], 'fishIDodor4', [], 'neuronIDcontrol', [], 'neuronIDodor1', [], 'neuronIDodor2', [], 'neuronIDodor3', [], 'neuronIDodor4', []);
    dataarray.(condition_tmp).stim.pump = setfield(dataarray.(condition_tmp).stim.pump, 'stimtriggeredtail', stimstruct);
    dataarray.(condition_tmp).stim.pump = setfield(dataarray.(condition_tmp).stim.pump, 'stimtriggeredcalcium', stimstruct);
    dataarray.(condition_tmp).stim.teensy = setfield(dataarray.(condition_tmp).stim.teensy, 'stimtriggeredtail', stimstruct);
    dataarray.(condition_tmp).stim.teensy = setfield(dataarray.(condition_tmp).stim.teensy, 'stimtriggeredcalcium', stimstruct);
    
    %bout-triggered
    struct_bout = struct( 'bouttriggeredcalcium', [], 'bouttriggeredtail', [], 'fishIDcalcium', [], 'fishIDtail', [], 'boutID', [], 'neuronID', [], 'stimID', [], 'neuronsperfish', []);
    dataarray.(condition_tmp).bout = setfield(dataarray.(condition_tmp).bout, 'poststim', struct_bout);
    dataarray.(condition_tmp).bout = setfield(dataarray.(condition_tmp).bout, 'nonstim', struct_bout);
    dataarray.(condition_tmp).bout = setfield(dataarray.(condition_tmp).bout, 'firstpoststim', struct_bout);
    
    
    dataarray.(condition_tmp).bout.poststim = struct('pump', struct_bout, 'teensy', struct_bout);
    dataarray.(condition_tmp).bout.nonstim = struct('pump', struct_bout, 'teensy', struct_bout);
    dataarray.(condition_tmp).bout.firstpoststim = struct('pump', struct_bout, 'teensy', struct_bout);
    
    struct_kinematics = struct('boutduration', [], 'peaklastpointangle', [], 'peakcumangle',[], 'power', [], 'maxvelocity', [], 'meanvelocity', [], 'symmetry', [], 'fishID', [], 'activityindex', []);
    dataarray.(condition_tmp).bout.poststim.pump = setfield(dataarray.(condition_tmp).bout.poststim.pump, 'kinematics', struct_kinematics);
    dataarray.(condition_tmp).bout.poststim.teensy = setfield(dataarray.(condition_tmp).bout.poststim.teensy, 'kinematics', struct_kinematics);
    dataarray.(condition_tmp).bout.nonstim.pump = setfield(dataarray.(condition_tmp).bout.nonstim.pump, 'kinematics', struct_kinematics);
    dataarray.(condition_tmp).bout.nonstim.teensy = setfield(dataarray.(condition_tmp).bout.nonstim.teensy, 'kinematics', struct_kinematics);
    dataarray.(condition_tmp).bout.firstpoststim.pump = setfield(dataarray.(condition_tmp).bout.firstpoststim.pump, 'kinematics', struct_kinematics);
    dataarray.(condition_tmp).bout.firstpoststim.teensy = setfield(dataarray.(condition_tmp).bout.firstpoststim.teensy, 'kinematics', struct_kinematics);
    
    
    %calciumtriggered
    
    %others
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'fishID', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'ROIsizes', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'zindices', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'ycoordinates', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'xcoordinates', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'dfoverf', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'dfoverf_smoothed', []);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_tailreg', []);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_motorstimreg', []);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_motorsponreg', []);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_stimreg', []);
    
end

count = 0;

for i =1:length(files)
    
    name = files(i).name;
    findstring = strfind(name, '.mat');
    findstring2 = strfind(name, '._');
    if isempty(findstring) || ~isempty(findstring2)
        continue;
    else
        filename = sprintf('%s/%s', pathname, name);
        
        %try
        load(filename)
         FISH = FISH+1;
        %catch
        %    continue
        %end
       
    end
    
    FISH

    exp_binary = [pump stim];
    
    for e = 1:length(exp_binary)
        
        if exp_binary(e) == 0
            break
        elseif exp_binary(e) == 1
            
            if e ==1
                stimtype = 'pump';
            elseif e == 2
                stimtype = 'teensy';
            end
            
            
            % BOUT TRIGGERED ANALYSES
            
            %to aggregate stim protocol/intensities etc.
            nstims = length(calciumtaildata.(condition).stim.(stimtype).stimtimestamps);

            try
                stimprotocol = calciumtaildata.(condition).stim.(stimtype).pumpprotocol;
                stimprotocolcopy = stimprotocol;
            catch
                stimintensities = calciumtaildata.(condition).stim.(stimtype).intensities;
                stimports = calciumtaildata.(condition).stim.(stimtype).stimports;
            end

            %if stimprotocol > 8; only use first 6 stims
            if length(stimprotocol)>8
                stimprotocolcopy = stimprotocol;
                stimprotocol = {};
                stimprotocol = stimprotocolcopy(1:6);
                nstims = length(stimprotocol);
            end

            % To aggregate tail and calcium parameters
            % Focus on post-stim or non post-stim bout

            poststimbouts = calciumtaildata.(condition).bout.poststimbouts.(stimtype);
            nonstimbouts = calciumtaildata.(condition).bout.nonstimbouts.(stimtype);
            idx_poststim = find(poststimbouts>0);

            %to obtain very first post stim bout

            if isempty(idx_poststim)
                idx_firstpoststim = [];
            else
                [C,ia,ic] = unique(poststimbouts);
                idx_firstpoststim = ia;
            end

            try
                idx_expstart =  find(calciumtaildata.(condition).bout.poststimbouts.teensy>0);
            catch
                idx_expstart = 1;
            end

            expstart = idx_expstart(1);
            expend = idx_expstart(end);

            idx_nonstim = find(nonstimbouts>0);
            idx_nonstim = idx_nonstim(idx_nonstim>=expstart & idx_nonstim<expend);

            stimcats = {};
            stimcats{1} = 'poststim';
            stimcats{2} = 'nonstim';
            stimcats{3} = 'firstpoststim';

            for c = 1:length(stimcats)

                stimcat = char(stimcats{c});

                stimID = {};
                if strcmp(stimcat, 'poststim')
                    idxcat = idx_poststim;
                    stimnums = poststimbouts(idxcat);
                elseif  strcmp(stimcat, 'firstpoststim')
                    idxcat = idx_firstpoststim;
                    stimnums = poststimbouts(idxcat);
                    stimnums = stimnums(stimnums>0);
                elseif   strcmp(stimcat, 'nonstim')
                    idxcat = idx_nonstim;
                    stimnums = nonstimbouts(idxcat);
                end

                try  %if pump

                    %if length(stimprotocolcopy)>8
                    % stimnums = stimnums(stimnums<length(stimprotocolcopy)/2); %less than 6 if it is 12stims total
                    % for s = 1:length(stimnums)
                    % stimID{s,1} =stimprotocol{stimnums(s)};
                    % end
                    % else
                    for s = 1:length(stimnums)
                        stimID{s,1} =stimprotocol{stimnums(s)};
                    end
                    %end

                catch %if teensy

                    stimID{s,1} =[];
                    %                         for s = 1:length(stimnums)
                    %                             if mod(stimnums(s), 2) ==1
                    %                                 stimID{s,1} =stimports(1);
                    %                             elseif mod(stimnums(s), 2) == 0
                    %                                 stimID{s,1} = stimports(2);
                    %                             end
                    %                         end

                end


                for s = 1:length(stimID)
                    dataarray.(condition).bout.(stimcat).(stimtype).stimID{end+1,1} = stimID{s};
                end

                % each cell array is a different bout parameter (also save fish ID)
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.boutduration(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.boutdurations(idxcat);  %bout durations
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.peaklastpointangle(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.peaklastpointangles(idxcat); %peak last point angle
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.maxvelocity(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.maxvelocity(idxcat);%maxvelocity
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.meanvelocity(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.meanvelocity(idxcat);%maxvelocity
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.peakcumangle(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.peakcumangles(idxcat);%peak cumulative angles
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.power(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.power(idxcat);%power
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.symmetry(end+1:end+length(idxcat),1) = calciumtaildata.(condition).bout.symmetry(idxcat);%symmetry
                dataarray.(condition).bout.(stimcat).(stimtype).kinematics.fishID(end+1:end+length(idxcat),1) = FISH*ones(1, length(idxcat));

                % to save tail and calcium traces around bout (bout-triggered)
                %before I start assigning, find where last bout count stopped
                lastbout = length(dataarray.(condition).bout.(stimcat).(stimtype).fishIDtail);

                if FISH == 1
                    lastneuron= 0;
                else
                    lastneuron = sum(dataarray.(condition).bout.(stimcat).(stimtype).neuronsperfish(1:FISH-1));
                end
                %will be zero if empty
                
                %fish ID saves fish identity
                %boutID saves bout identity
                %neuronID saves neuron identity
               
                boutcalcium = calciumtaildata.(condition).bout.bouttriggeredcalcium;
                bouttail = calciumtaildata.(condition).bout.bouttriggeredtail;
                
                extra =  find(idxcat > size(boutcalcium, 1)); %because I do not save bout traces if they end past the end of experiment
                if isempty(extra)
                else
                    idxcat(extra) = nan;
                    idxcat = idxcat(~isnan(idxcat));
                end
                
                dataarray.(condition).bout.(stimcat).(stimtype).bouttriggeredtail(end+1:end+length(idxcat), 1:size(bouttail,2)) = bouttail(idxcat, :);
                dataarray.(condition).bout.(stimcat).(stimtype).fishIDtail(end+1:end+length(idxcat),1) = FISH*ones(length(idxcat),1);
                
                calciumarray4 = boutcalcium(idxcat, :,:);
                
                timewindow2 = calciumtaildata.fishinfo.timewindow2;
                if strcmp(stimtype, 'pump')
                    tracelength = timewindow2*2; %this makes it convenient so timewindow = frames 
                    tracelength2 = 15;
                elseif strcmp(stimtype, 'teensy')
                    tracelength = 15;
                end

                %for each calcium trace, which tail trace is associated with it
                boutID = []; neuronID = [];
                for b = 1: size(calciumarray4, 1)
                    boutID(end+1:end+size(calciumarray4, 2),1) = (lastbout+b)*ones(size(calciumarray4, 2),1);
                    neuronID(end+1:end+size(calciumarray4, 2)) = lastneuron+1: lastneuron + size(calciumarray4,2);
                end
                
                dataarray.(condition).bout.(stimcat).(stimtype).neuronsperfish(end+1) = size(calciumarray4,2);
                
                %now reshape.
                
                calciumarraynew = [];
                for row = 1:size(calciumarray4,1)
                    for column = 1:size(calciumarray4,2)
                        calciumarraynew(end+1,1:size(calciumarray4,3)) = calciumarray4(row, column, :);
                    end
                end
                try
                    %if different number of zplanes -- downsample
                    if size(calciumarraynew,2) > tracelength2 || size(calciumarraynew,2) < tracelength2  %this means there are 2 planes imaged instead of 3;
                        calciumarraynew2 =  resample(calciumarraynew', tracelength2, size(calciumarraynew,2));
                        calciumarraynew = calciumarraynew2';
                    end
                catch
                end
                
                %calciumarray4 = reshape(calciumarray4, [size(calciumarray4,1)*size(calciumarray4,2), size(calciumarray4,3)]);
                dataarray.(condition).bout.(stimcat).(stimtype).bouttriggeredcalcium(end+1:end+size(calciumarraynew,1), 1:size(calciumarraynew,2)) = calciumarraynew;
                dataarray.(condition).bout.(stimcat).(stimtype).fishIDcalcium(end+1:end+size(calciumarraynew,1),1) = FISH*ones(size(calciumarraynew,1),1);
                dataarray.(condition).bout.(stimcat).(stimtype).boutID(end+1:end+length(boutID), 1) = boutID;
                dataarray.(condition).bout.(stimcat).(stimtype).neuronID(end+1:end+length(neuronID), 1) = neuronID;
                
            end
            
            % STIM TRIGGERED ANALYSES
            
            %to aggregate tail angle and mean dfover f data for correlation analysis.
            %Keep track of stim number.
            
            dataarray.(condition).stim.(stimtype).lastpointangle(FISH, 1:nstims) = calciumtaildata.(condition).stim.(stimtype).peaklastpointangle(1:nstims);
            dataarray.(condition).stim.(stimtype).cumangle(FISH, 1:nstims) = calciumtaildata.(condition).stim.(stimtype).peakcumangle(1:nstims);
            dataarray.(condition).stim.(stimtype).power(FISH, 1:nstims) = calciumtaildata.(condition).stim.(stimtype).power(1:nstims);
            dataarray.(condition).stim.(stimtype).calciumpeaks(FISH, 1:nstims) = mean(calciumtaildata.(condition).stim.(stimtype).peakcalcium(1:nstims),2);
            dataarray.(condition).stim.(stimtype).calciumsum(FISH, 1:nstims) = mean(calciumtaildata.(condition).stim.(stimtype).sumcalcium(1:nstims),2);
            
            
            try
                dataarray.(condition).stim.(stimtype).protocol(FISH, 1:nstims) = stimprotocol;
            catch
                dataarray.(condition).stim.(stimtype).protocol(FISH, 1:nstims) = stimports;
                dataarray.(condition).stim.(stimtype).intensities(FISH, 1:nstims) = stimintensities;
            end
            
            %To aggregate calcium and tail arrays surrounding stim by intensity
            %Fish identity saved as FishID
            c2 = calciumtaildata.(condition).stim.(stimtype).stimtriggeredtail;
            calcium = calciumtaildata.(condition).stim.(stimtype).stimtriggeredcalcium;
            
            calciumtraces_stim = [];
            
            %ODOR STIM 1
            try idx_odor = find(strcmp(stimprotocol, firstodor)==1);
                
                tailarray1 = c2(idx_odor,:);
                calciumarray1 = calcium(idx_odor, :, :);
                
                
                %for each stimulus-triggered trace, which neurons are associated with it
                neuronID = [];
                for s = 1: size(idx_odor, 1)
                    neuronID(end+1:end+size(calciumarray1, 2)) = lastneuron+1: lastneuron + size(calciumarray1,2);
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor1(end+1:end+length(neuronID), 1) = neuronID;
                
                calciumtraces_stim = [];
                for row = 1:size(calciumarray1,1)
                    for column = 1:size(calciumarray1,2)
                        calciumtraces_stim(end+1,1:size(calciumarray1,3)) = calciumarray1(row, column, :);
                    end
                end
                
                %if different number of zplanes -- downsample
                if size(calciumtraces_stim,2) >tracelength || size(calciumtraces_stim,2) <tracelength %this means there are 2 planes imaged instead of 3;
                    calciumtraces_stim2 =  resample(calciumtraces_stim', tracelength, size(calciumtraces_stim,2));
                    calciumtraces_stim = calciumtraces_stim2';
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.odor1(end+1:end+size(tailarray1,1), 1:length(tailarray1)) = tailarray1;
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.odor1(end+1:end+size(calciumtraces_stim,1), 1:size(calciumtraces_stim,2))= calciumtraces_stim;
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDodor1(end+1:end+size(tailarray1,1),1) = FISH*ones(size(tailarray1,1),1);
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.fishIDodor1(end+1:end+size(calciumtraces_stim,1), 1) = FISH*ones(size(calciumtraces_stim,1),1);
                
            catch
            end
            
            %ODOR STIM 2
            try idx_odor = find(strcmp(stimprotocol, secondodor)==1);
                
                tailarray1 = c2(idx_odor,:);
                calciumarray1 = calcium(idx_odor, :, :);
                
                
                %for each stimulus-triggered trace, which neurons are associated with it
                neuronID = [];
                for s = 1: size(idx_odor, 1)
                    neuronID(end+1:end+size(calciumarray1, 2)) = lastneuron+1: lastneuron + size(calciumarray1,2);
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor2(end+1:end+length(neuronID), 1) = neuronID;
                
                calciumtraces_stim = [];
                for row = 1:size(calciumarray1,1)
                    for column = 1:size(calciumarray1,2)
                        calciumtraces_stim(end+1,1:size(calciumarray1,3)) = calciumarray1(row, column, :);
                    end
                end
                
                %if different number of zplanes -- downsample
                if size(calciumtraces_stim,2) >tracelength || size(calciumtraces_stim,2) < tracelength %this means there are 2 planes imaged instead of 3;
                    calciumtraces_stim2 =  resample(calciumtraces_stim', tracelength, size(calciumtraces_stim,2));
                    calciumtraces_stim = calciumtraces_stim2';
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.odor2(end+1:end+size(tailarray1,1), 1:length(tailarray1)) = tailarray1;
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.odor2(end+1:end+size(calciumtraces_stim,1), 1:size(calciumtraces_stim,2))= calciumtraces_stim;
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDodor2(end+1:end+size(tailarray1,1),1) = FISH*ones(size(tailarray1,1),1);
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.fishIDodor2(end+1:end+size(calciumtraces_stim,1), 1) = FISH*ones(size(calciumtraces_stim,1),1);
                
            catch
            end
            
            %ODOR STIM 3
            try idx_odor = find(strcmp(stimprotocol, thirdodor)==1);
                
                tailarray1 = c2(idx_odor,:);
                calciumarray1 = calcium(idx_odor, :, :);
                
                
                %for each stimulus-triggered trace, which neurons are associated with it
                neuronID = [];
                for s = 1: size(idx_odor, 1)
                    neuronID(end+1:end+size(calciumarray1, 2)) = lastneuron+1: lastneuron + size(calciumarray1,2);
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor3(end+1:end+length(neuronID), 1) = neuronID;
                
                calciumtraces_stim = [];
                for row = 1:size(calciumarray1,1)
                    for column = 1:size(calciumarray1,2)
                        calciumtraces_stim(end+1,1:size(calciumarray1,3)) = calciumarray1(row, column, :);
                    end
                end
                
                %if different number of zplanes -- downsample
                if size(calciumtraces_stim,2) >tracelength || size(calciumtraces_stim,2) <tracelength %this means there are 2 planes imaged instead of 3;
                    calciumtraces_stim2 =  resample(calciumtraces_stim', tracelength, size(calciumtraces_stim,2));
                    calciumtraces_stim = calciumtraces_stim2';
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.odor3(end+1:end+size(tailarray1,1), 1:length(tailarray1)) = tailarray1;
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.odor3(end+1:end+size(calciumtraces_stim,1), 1:size(calciumtraces_stim,2))= calciumtraces_stim;
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDodor3(end+1:end+size(tailarray1,1),1) = FISH*ones(size(tailarray1,1),1);
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.fishIDodor3(end+1:end+size(calciumtraces_stim,1), 1) = FISH*ones(size(calciumtraces_stim,1),1);
                
            catch
            end
            
            %ODOR STIM 4
            
            try idx_odor = find(strcmp(stimprotocol, fourthodor)==1);
                
                tailarray1 = c2(idx_odor,:);
                calciumarray1 = calcium(idx_odor, :, :);
                
                %for each stimulus-triggered trace, which neurons are associated with it
                neuronID = [];
                for s = 1: size(idx_odor, 1)
                    neuronID(end+1:end+size(calciumarray1, 2)) = lastneuron+1: lastneuron + size(calciumarray1,2);
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor4(end+1:end+length(neuronID), 1) = neuronID;
                
                calciumtraces_stim = [];
                for row = 1:size(calciumarray1,1)
                    for column = 1:size(calciumarray1,2)
                        calciumtraces_stim(end+1,1:size(calciumarray1,3)) = calciumarray1(row, column, :);
                    end
                end
                
                %if different number of zplanes -- downsample
                if size(calciumtraces_stim,2) >tracelength || size(calciumtraces_stim,2) <tracelength %this means there are 2 planes imaged instead of 3;
                    calciumtraces_stim2 =  resample(calciumtraces_stim', tracelength, size(calciumtraces_stim,2));
                    calciumtraces_stim = calciumtraces_stim2';
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.odor4(end+1:end+size(tailarray1,1), 1:length(tailarray1)) = tailarray1;
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.odor4(end+1:end+size(calciumtraces_stim,1), 1:size(calciumtraces_stim,2))= calciumtraces_stim;
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDodor4(end+1:end+size(tailarray1,1),1) = FISH*ones(size(tailarray1,1),1);
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.fishIDodor4(end+1:end+size(calciumtraces_stim,1), 1) = FISH*ones(size(calciumtraces_stim,1),1);
                
            catch
            end
            
            %CONTROL STIM
            try idx_control = find(strcmp(stimprotocol, controlodor)==1);
                
                tailarray1 = c2(idx_control,:);
                calciumarray1 = calcium(idx_control, :, :);
                
                
                %for each stimulus-triggered trace, which neurons are associated with it
                neuronID = [];
                for s = 1: size(idx_control, 1)
                    neuronID(end+1:end+size(calciumarray1, 2)) = lastneuron+1: lastneuron + size(calciumarray1,2);
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDcontrol(end+1:end+length(neuronID), 1) = neuronID;
                
                calciumtraces_stim = [];
                for row = 1:size(calciumarray1,1)
                    for column = 1:size(calciumarray1,2)
                        calciumtraces_stim(end+1,1:size(calciumarray1,3)) = calciumarray1(row, column, :);
                    end
                end
                
                %if different number of zplanes -- downsample
                if size(calciumtraces_stim,2) >tracelength || size(calciumtraces_stim,2) <tracelength%this means there are 2 planes imaged instead of 3;
                    calciumtraces_stim2 =  resample(calciumtraces_stim', tracelength, size(calciumtraces_stim,2));
                    calciumtraces_stim = calciumtraces_stim2';
                end
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.control(end+1:end+size(tailarray1,1), 1:length(tailarray1)) = tailarray1;
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.control(end+1:end+size(calciumtraces_stim,1), 1:size(calciumtraces_stim,2))= calciumtraces_stim;
                
                dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDcontrol(end+1:end+size(tailarray1,1),1) = FISH*ones(size(tailarray1,1),1);
                dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.fishIDcontrol(end+1:end+size(calciumtraces_stim,1), 1) = FISH*ones(size(calciumtraces_stim,1),1);
                
            catch
            end
            
        end
        
    end
    
    % dfoverf and others
    
    %to aggregate fishID of each neuron
    nneurons = length(calciumtaildata.(condition).spatial.ROIsize);
    nneurons = calciumtaildata.fishinfo.number_of_neurons;
    dataarray.(condition).calcium.fishID(end+1:end+nneurons, 1) = FISH*ones(nneurons,1);
    
    idxnotnan  = find(~isnan(calciumtaildata.(condition).spatial.APpositions));
    % to aggregate ROI sizes
    dataarray.(condition).calcium.ROIsizes(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.ROIsize(1:nneurons);
    
    %to aggregate indices of z-plane numbers
    dataarray.(condition).calcium.zindices(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.DVpositions(1:nneurons);
    
    %to aggregate coordinates
    dataarray.(condition).calcium.ycoordinates(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.APpositions(idxnotnan);
    dataarray.(condition).calcium.xcoordinates(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.LRpositions(idxnotnan);
    
    %XCORR
    dataarray.(condition).calcium.xcorr_tailreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_tailreg;
    dataarray.(condition).calcium.xcorr_motorstimreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_motorstimreg;
    dataarray.(condition).calcium.xcorr_motorsponreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_motorsponreg;
    dataarray.(condition).calcium.xcorr_stimreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_stimreg;

    % aggregate regressors according to name

    stimreg_name = calciumtaildata.(condition).calcium.xcorr_stimregname;

    for c = 1:length(stimreg_name)
        
        try
            dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))(count+1:count+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg(:,c);
            
        catch
            dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))(count+1:count+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg(:,c);
        end
        
    end
 
    %make the zeros nans
    for c = 1:length(stimreg_name)
        dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))(dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))==0) = nan;
    end
    

        count = length(dataarray.(condition).calcium.xcorr_stimreg.water);
        
        %to aggregate dfoverf
        dfoverf_consolidated = calciumtaildata.(condition).calcium.dfoverf(idxnotnan, :);
        
        if nchannel == 2
            dfoverf_RED_consolidated = calciumtaildata.(condition).calcium.dfoverf_RED(idxnotnan, :);
        end
        
        dfoverf_smoothed = calciumtaildata.(condition).calcium.dfoverf_smoothed(idxnotnan,:);
        
        
        if size(dfoverf_smoothed,2) > minframes || size(dfoverf_smoothed,2) < minframes %this is miminum number of frames
            dfoverf_smoothed2 = resample(dfoverf_smoothed', minframes, size(dfoverf_smoothed, 2));
            dfoverf_smoothed = dfoverf_smoothed2';
        end
        
        dataarray.(condition).calcium.dfoverf_smoothed(end+1:end+nneurons,:) = dfoverf_smoothed;
        
        
        if size(dfoverf_consolidated,2) > minframes || size(dfoverf_consolidated,2) < minframes %this is miminum number of frames
            dfoverf_consolidated2 = resample(dfoverf_consolidated', minframes, size(dfoverf_consolidated, 2));
            dfoverf_consolidated = dfoverf_consolidated2';
        end
        
        dataarray.(condition).calcium.dfoverf(end+1:end+nneurons,:) = dfoverf_consolidated;
        
        if nchannel == 2
            
            if size(dfoverf_RED_consolidated,2) > minframes || size(dfoverf_RED_consolidated,2) < minframes %this is miminum number of frames
                dfoverf_RED_consolidated2 = resample(dfoverf_RED_consolidated', 1770, size(dfoverf_RED_consolidated, 2));
                dfoverf_RED_consolidated = dfoverf_RED_consolidated2';
            end
            
            dataarray.(condition).calcium.dfoverf_RED(end+1:end+nneurons,:) = dfoverf_RED_consolidated;
        end
end

clear calciumtaildata

% clearvars -except dataarray celltype exps folderpath_aggregated exptype disk nchannel
%
f = 15;
lw = 2;
zstart = 1;
zend = 3;

% replace zeroes with nans

stimregarray = fieldnames(dataarray.(condition).calcium.xcorr_stimreg);

for c = 1:length(stimregarray)
    
    cat = stimregarray{c}
    dataarray.(condition).calcium.xcorr_stimreg.(cat)(dataarray.(condition).calcium.xcorr_stimreg.(cat) ==0) = nan;
    
    if length(dataarray.(condition).calcium.xcorr_stimreg.(cat)) < length(dataarray.(condition).calcium.xcorr_stimreg.water)
    dataarray.(condition).calcium.xcorr_stimreg.(cat)(end+1:length(dataarray.(condition).calcium.xcorr_stimreg.water)) = nan;
    
    end
end

%% SECTION 2: THIS IS WHERE ANALYSIS STARTS
%requires dataarray

% BOUT KINEMATICS PANEL
%per fish
stimcat = 'poststim'
stimtype = 'pump'
stimchannel = 'kin'
ymax = 0.4;
kinematicsstats = [];
meankinematics = {};

fnames = fieldnames(dataarray.(condition).bout.(stimcat).(stimtype).kinematics);
parameters =  [3 4 5 6];

figure(1)
clf;

m = 1

for i = 1:length(parameters)

    parameter = char(fnames(parameters(i)));

    %extract tail bouts
    tail = dataarray.(condition).bout.(stimcat).(stimtype).kinematics.(parameter);
    stimID = dataarray.(condition).bout.(stimcat).(stimtype).stimID;
    fishID = dataarray.(condition).bout.(stimcat).(stimtype).kinematics.fishID;

    if strcmp(parameter,'symmetry')
        for j = 1:length(tail)
            if tail(j)<0.5
                tail(j) = 1-DMSO_tail(j);
            end
        end
    end

    %per fish
    fishmeans_control = [];
    fishmeans_treatment = [];

    for fish = 1:max(fishID)
        fishmeans_control(fish) = nanmean(tail(fishID==fish & strcmp(stimID, 'water')));
        fishmeans_treatment(fish) = nanmean(tail(fishID==fish & strcmp(stimID, stimchannel)));
    end

    [p,h] = signrank(fishmeans_control, fishmeans_treatment);

    kinematicsstats(i,1) = p;
    kinematicsstats(i,2) = nanmean(fishmeans_control);
    kinematicsstats(i,3) = nanmean(fishmeans_treatment);
    meankinematics{i}(:,1) = fishmeans_control;
    meankinematics{i}(:,2) = fishmeans_treatment;

    figure(1)
    subplot(length(parameters),2,i*2-1)
    scatter(ones(1, length(fishmeans_control)), fishmeans_control, 'bo', 'filled');
    hold on;
    scatter(2*ones(1, length(fishmeans_treatment)), fishmeans_treatment, 'ro', 'filled');

    for fish = 1:max(fishID)
        plot([1 2], [fishmeans_control(fish), fishmeans_treatment(fish)], 'k');
    end
    ylabel(strcat(parameter, ' p=', num2str(p)));
    box off

    if strcmp(parameter,'symmetry')
        axis([0 3 0.5 1]);
    elseif strcmp(parameter, 'boutduration')
        axis([0 3 0 600])
    elseif strcmp(parameter, 'power')
        axis([0 3 0 40000])
    elseif strcmp(parameter, 'maxvelocity')
        axis([0 3 0 50]);
    elseif strcmp(parameter, 'meanvelocity')
        axis([0 3 0 5]);
    elseif strcmp(parameter, 'peakcumangle')
        axis([0 3 0 250]);
    elseif strcmp(parameter, 'frequency')
        axis([0 3 0 50]);
    end
    
end

name = strcat('meankinematicsperfish-', stimchannel);
print(name, '-depsc2');

%output is kinematicsstats and meankinematics


%% Distribution of bouts 
stimcat = 'poststim'
stimtype = 'pump'
stimchannel = 'kin'
ymax = 0.4;

fnames = fieldnames(dataarray.(condition).bout.(stimcat).(stimtype).kinematics);
parameters =  [1 3 4 5 6];

figure(2);
clf;

tsnearray = [];

for i = 1:length(parameters)

    parameter = char(fnames(parameters(i)));

    %poststim

    tail =  dataarray.(condition).bout.(stimcat).(stimtype).kinematics.(parameter);
    stimID = dataarray.(condition).bout.(stimcat).(stimtype).stimID;
    control_tail = tail(find(strcmp(stimID, 'water')));
    odor_tail = tail((find(strcmp(stimID, stimchannel))));

    figure(2)
    subplot(length(parameters), 2, i*2-1)

    if strcmp(parameter, 'boutduration')
        dx =40;
    elseif strcmp(parameter, 'peakcumangle')
        dx =15;
    elseif strcmp(parameter, 'power')
        dx =500;
    elseif strcmp(parameter, 'symmetry')
        dx = 0.5;
    elseif strcmp(parameter, 'maxvelocity')
        dx = 2.5;
    elseif strcmp(parameter, 'meanvelocity')
        dx = 0.25;

    end

    xvalues = min(control_tail):dx:max(control_tail);
    nelements  = hist(control_tail,xvalues);
    
    try
        
    histogram(control_tail, xvalues, 'Normalization', 'probability');

    h=findobj(gca, 'Type', 'patch');
    set(h, 'Facecolor', 'none', 'EdgeColor', 'k');
    box off;

    if strcmp(parameter, 'boutduration')
        axis([0 900 0 ymax])
    elseif strcmp(parameter, 'peakcumangle')
        axis([0 360 0 ymax])
    elseif strcmp(parameter, 'power')
        axis([0 10000 0 ymax])
    elseif strcmp(parameter, 'maxvelocity')
        axis([0 80 0 ymax])
    elseif strcmp(parameter, 'meanvelocity')
        axis([0 5 0 ymax]);
    end

    ylabel(parameter);

    figurename = strcat('kinematicshistogram-control-', stimcat, stimtype);
    print(figurename, '-depsc2');

    figure(3)
    subplot(length(parameters), 2, i*2-1)
    xvalues = min(odor_tail):dx:max(odor_tail);
    nelements  = hist(odor_tail,xvalues);
    histogram(odor_tail, xvalues, 'Normalization', 'probability');

    if strcmp(parameter, 'boutduration')
        axis([0 900 0 ymax])
    elseif strcmp(parameter, 'peakcumangle')
        axis([0 360 0 ymax])
    elseif strcmp(parameter, 'power')
        axis([0 10000 0 ymax])
    elseif strcmp(parameter, 'maxvelocity')
        axis([0 80 0 ymax])
    elseif strcmp(parameter, 'meanvelocity')
        axis([0 5 0 ymax]);
    end

    h=findobj(gca, 'Type', 'patch');
    set(h, 'Facecolor', 'none', 'EdgeColor', 'r');
    box off;
    if strcmp(stimcat, 'poststim') || strcmp(stimcat, 'firstpoststim')
    [p,h] = ranksum(control_tail, odor_tail,'tail', 'right');
    elseif strcmp(stimcat, 'nonstim')
    [p,h] = ranksum(control_tail, odor_tail, 'tail', 'right'); 
    end
    
    n1 = length(control_tail);
    n2 = length(odor_tail);
    ylabel(strcat('p=', num2str(p), "n=", num2str(n1), '/', num2str(n2)));
    mean1 = mean(control_tail);
    mean2 = mean(odor_tail);
    err1 = std(control_tail)/sqrt(n1);
    err2 = std(odor_tail)/sqrt(n2);

    xlabel(strcat(parameter, '',num2str(mean1), '+-',  num2str(err1), '/', num2str(mean2), '+-',  num2str(err2)));

    figurename = strcat('kinematicshistogram-odor-', stimcat, stimtype, stimchannel);
    print(figurename, '-depsc2');

    catch
    end


    tsnearray = horzcat(tsnearray, vertcat(control_tail, odor_tail));
end

%% STIM TRIGGERED AVERAGES

close all;
stimtype = 'pump'

if strcmp(condition, 'optovin')
    ymin = -0.8;
    ymax =1.8;
else
    ymin = -0.4;
    ymax = 0.5;
end

dataarray.(condition).stim.(stimtype).tailkinematics = struct;

stimchannel = {};

if strcmp(condition, 'optovin')
    stimchannel{1} = 'control';
    stimchannel{2} = 'odor1';
else
    stimchannel{1} = 'control';
    stimchannel{2} = 'odor1';
    stimchannel{3} = 'odor2';
    stimchannel{4} = 'odor3';
end


calciumsum_pre = [];
calciumsum_post = [];

calciummax_pre = [];
calciummax_post = [];

calciumsum_prepump = [];
calciumsum_postpump = [];

trace_byneuron = {};


if strcmp(condition, 'water')
    window = 89; %this is the range in which you average calcium responses
    windowstart = 0; %start quantifying after stimulus has been presented
    statswindow = 60; %60 seconds
elseif strcmp(condition, 'optovin')
    window = 5;
    windowstart = 0;
end

if strcmp(stimtype, 'teensy')
    window = 5;
end

figure(4);
clf;

map = [0 0 0; 1 140/255 0; 30/255 144/255 1 ;199/255 21/255 133/255; 50/255 205/255 50/255; 0 0 1]

statstable = {};
abs_tailtraces = struct;
abs_tailtraces.tails = struct;
abs_tailtraces.fish = struct;


for s = 1:length(stimchannel)

    stimchan = char(stimchannel{s});
    calcium_raw = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.(stimchan);
    ypos = dataarray.(condition).calcium.ycoordinates;
    xpos = dataarray.(condition).calcium.xcoordinates;

    figure(6)
    scatter(xpos, ypos);

    dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan) = struct;

    abs_tail = abs(dataarray.(condition).stim.(stimtype).stimtriggeredtail.(stimchan));
    tail = dataarray.(condition).stim.(stimtype).stimtriggeredtail.(stimchan);

    frameunit = mean(diff(calciumtaildata.(condition).tail.tailtime))*1000;
    timeframe = round(window*1000/frameunit); %postcue

    % extract kinematics
    for j = 1:size(abs_tail, 1)

        %STIM TYPE NEEDS TO BE PUMP
        if strcmp(stimtype, 'pump')

            tmp = abs_tail(j, size(abs_tail,2)/2:size(abs_tail,2)/2+timeframe);
            tmp_pre = abs_tail(j, 1:size(abs_tail,2)/2);

            %find peaks using non-smoothed trace (for peak angle, latencies)
            [pks, locs] = findpeaks(tmp, 'MinPeakHeight', 25,'MinPeakDistance', 25, 'MinPeakProminence', 25);

            %find peaks using smoothed trace (for counting bouts)
            tmp1 = smooth(tmp, 5);
            [pks1, locs1] = findpeaks(tmp1, 'MinPeakHeight',30,'MinPeakDistance', 25, 'MinPeakProminence', 30);

            [pks_pre, locs_pre] = findpeaks(tmp_pre, 'MinPeakHeight', 25 ,'MinPeakDistance', 25, 'MinPeakProminence', 25);

            tmp_pre1 = smooth(tmp_pre, 5);
            [pks_pre1, locs_pre1] = findpeaks(tmp_pre1, 'MinPeakHeight', 30 ,'MinPeakDistance', 25, 'MinPeakProminence', 30);

            %use not smoothed trace
            %peak angle
            try
                meanpks = mean(pks);
                peakangles( j,1) = meanpks;
            catch
                peakangles( j,1) = 0;
            end

            try
                meanpks_pre = mean(pks_pre);
                peakangles_pre(j,1) = meanpks_pre;
            catch
                peakangles_pre(j,1) = 0;
            end

            %uses smoothed trace, bouts per second

            boutthreshold = 50;
            maxboutthreshold = 100;
            nbouts(j,1) = length(pks1) / window;
            nbouts_large(j,1) = length(pks1(pks1>= maxboutthreshold)) / window;
            nbouts_medium(j,1) = length(pks1(pks1>=boutthreshold & pks1<maxboutthreshold))/ window;
            nbouts_small(j,1) = length(pks1(pks1<boutthreshold))/window;
            nbouts_large2(j,1) = length(pks1(pks1>=boutthreshold))/window;

            try
                nbouts_pre(j,1) = length(pks_pre1)/window;
                nbouts_pre_large(j,1) = length(pks_pre1(pks_pre1>=maxboutthreshold)) / window;
                nbouts_pre_medium(j,1) = length(pks_pre1(pks_pre1>=boutthreshold & pks_pre1<maxboutthreshold))/ window;
                nbouts_pre_small(j,1) = length(pks_pre1(pks_pre1<boutthreshold)) / window;
                nbouts_pre_large2(j,1) = length(pks_pre1(pks_pre1>=boutthreshold)) / window;
            catch
                nbouts_pre(j,1) = 0;
                nbouts_pre_large(j,1) = 0;
                nbouts_pre_medium(j,1) = 0;
                nbouts_pre_small(j,1) = 0;
                nbouts_pre_large2(j,1) = 0;
            end

            %using normal trace
            sumangle = sum(tmp); %up to end of bouts
            sumangles(j,1) = sumangle;

            try
                sumangle_pre = sum(tmp_pre);
                sumangles_pre(j,1) = sumangle_pre;
            catch
                sumangle_pre = 0;
                sumangles_pre(j,1) = sumangle_pre;
            end

            %velocity
            try
                maxvelocity = max(abs(diff(tmp(1:locs(2)),1,2)), [], 2)/frameunit; %for first poststim bout
                meanvelocity = nanmean(abs(diff(tmp(1:locs(2)),1,2)), 2)/frameunit;%for first poststim bout
                maxvelocities(j, 1) = maxvelocity;
                meanvelocities( j, 1) = meanvelocity;
                
            catch
                maxvelocity = 0;
                meanvelocity = 0;
                maxvelocities(j, 1) = maxvelocity;
                meanvelocities( j, 1) = meanvelocity;
            end
            
            try
                maxvelocity_pre = max(abs(diff(tmp_pre,1,2)), [], 2)/frameunit;
                meanvelocity_pre = nanmean(abs(diff(tmp_pre,1,2)), 2)/frameunit;
                maxvelocities_pre(j, 1) =maxvelocity_pre;
                meanvelocities_pre( j, 1) = meanvelocity_pre;
            catch
                maxvelocity_pre = 0;
                meanvelocity_pre = 0;
                maxvelocities_pre(j, 1) =maxvelocity_pre;
                meanvelocities_pre( j, 1) = meanvelocity_pre;
            end

        end

        %pre
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).peakanglepre = peakangles_pre;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).sumanglepre = sumangles_pre;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).maxvelocitypre = maxvelocities_pre;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).meanvelocitypre = meanvelocities_pre;

        %post
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).peakangle = peakangles;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).sumangle = sumangles;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).maxvelocity = maxvelocities;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).meanvelocity = meanvelocities;

        %number of bouts

        %pre
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqpre = nbouts_pre;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqprealllarge = nbouts_pre_large2; %large2 if you include medium
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqprevlarge = nbouts_pre_large;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqpremedium = nbouts_pre_medium;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqpresmall = nbouts_pre_small;

        %post
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreq = nbouts;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqalllarge = nbouts_large2; %large2 if you include medium
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqvlarge = nbouts_large;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqmedium = nbouts_medium;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).boutfreqsmall = nbouts_small;

    end

    if strcmp(stimchan, 'control')
        neuronID = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDcontrol;
    elseif strcmp(stimchan, 'odor1')
        neuronID = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor1;
    elseif strcmp(stimchan, 'odor2')
        neuronID = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor2;
    elseif strcmp(stimchan, 'odor3')
        neuronID = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor3;
    elseif strcmp(stimchan, 'odor4')
        neuronID = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDodor4;
    end

    neuronIDcontrol = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.neuronIDcontrol;
    %same for control or stim
    fishID = dataarray.(condition).stim.(stimtype).stimtriggeredtail.fishIDcontrol;

    totaltimepoints = size(calcium_raw,2);
    totaltailpoints = size(abs_tail,2);
    baselineend = size(calcium_raw,2)/2;

    calcium_raw(calcium_raw==0) = nan;
    calcium = calcium_raw -nanmean(calcium_raw(:, 1:baselineend), 2);

    abs_tailtraces.tails.(stimchan) = [];
    abs_tailtraces.fish.(stimchan) = [];

    for stim = 1:size(abs_tail,1)
        if nanmean(abs_tail(stim,:))==0
        else
            abs_tailtraces.tails.(stimchan)(end+1, 1:size(abs_tail, 2)) = abs_tail(stim, 1:size(abs_tail, 2))
            abs_tailtraces.fish.(stimchan)(end+1) = fishID(stim);
        end
    end

    %this is to define range over which to quantify response

    stimstart = 1;
    if strcmp(condition, 'optovin')
        stimstart = size(calcium,2)/6
    end

    %this is just to average activity for each neuron during each stim.
    %Not sorted by neuron
    for n = 1:length(neuronID)

        calciumsum_post(n,s) = nansum(calcium(n, baselineend+stimstart+windowstart:baselineend+stimstart + windowstart+ window))/window; %average response/s
        calciumsum_pre(n,s) = nansum(calcium(n, baselineend+stimstart-window:baselineend+stimstart-1))/window; %average response/s
        calciummax_post(n,s) = max(calcium(n, baselineend+stimstart:baselineend+stimstart+window)); %maximum pre stim
        calciummax_pre(n,s) = max(calcium(n, baselineend+stimstart-window:baselineend+stimstart-1)); %maximum post stim

        if strcmp(condition, 'optovin')
            calciumsum_prepump(n,s) =  nansum(calcium(n, baselineend-window-1:baselineend+window))/window; %average response/s
            %this is currently same as calciumsum_pre
            calciumsum_postpump(n,s) =  nansum(calcium(n, baselineend+stimstart-window:baselineend+stimstart-1))/window*3; %average response/s
        end

    end

    calciumsum_pre(calciumsum_pre ==0) = nan;
    calciumsum_post(calciumsum_post==0) = nan;
    calciumsum_prepump(calciumsum_prepump==0) = nan;
    calciumsum_postpump(calciumsum_postpump==0) = nan;
      
    % now average both calcium activity and average activity by neuron
    for n = 1:max(neuronIDcontrol)
        idx_n = find(neuronID == n);
        
        calciumsum_post_byneuron(n,s) = nanmean(calciumsum_post(idx_n,s));
        calciumsum_pre_byneuron(n,s) = nanmean(calciumsum_pre(idx_n,s));
        calciummax_post_byneuron(n,s) = nanmean(calciummax_post(idx_n,s));
        calciummax_pre_byneuron(n,s) = nanmean(calciummax_pre(idx_n,s));

        trace_byneuron{s}(n,1:size(calcium,2)) = nanmean(calcium(idx_n,:),1);

        if strcmp(condition, 'optovin')
            calciumsum_prepump_byneuron(n,s)= nanmean(calciumsum_prepump(idx_n,s));
            calciumsum_postpump_byneuron(n,s)= nanmean(calciumsum_postpump(idx_n,s));
        end

    end

    calcium = trace_byneuron{s};
    tmp = calcium(:,1);
    nneurons = length(tmp(~isnan(tmp)));

    calcium_norm = trace_byneuron{s}-trace_byneuron{1};

    %specific to OXT neurons
    anteriorthreshold = 175;

    populations = {'PO', 'PT'};

    for j = 1:length(populations)

        population = populations{j};

        if j==1
            calcium = calcium(ypos <= anteriorthreshold, :);
            calcium_norm = calcium_norm(ypos <= anteriorthreshold,:);
            calcium_forstats = calciumsum_post_byneuron(ypos <= anteriorthreshold,:);
            nneurons = size(calcium_norm, 1);
        else
            calcium = calcium(ypos > anteriorthreshold, :);
            calcium_norm = calcium_norm(ypos > anteriorthreshold, :);
            calcium_forstats = calciumsum_post_byneuron(ypos > anteriorthreshold,:);
            nneurons = size(calcium_norm, 1);
        end

        writematrix(calcium, strcat(stimchan, '-', population,'-stimtriggeredaverage.csv'));

        if j == 1
            figure(4);
        else
            figure(5);
        end
        subplot(2,2,1);
        plot(1:size(calcium,2), nanmean(calcium,1), 'Color', map(s, :));
        hold on;
        %errorbars
        plot(1:size(calcium,2), nanmean(calcium,1)+nanstd(calcium,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);
        plot(1:size(calcium,2), nanmean(calcium,1)-nanstd(calcium,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);

        line([size(calcium,2)/2 size(calcium,2)/2], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);

        if strcmp(condition, 'optovin')
            line([size(calcium,2)/2+size(calcium,2)/6  size(calcium,2)/2 + size(calcium,2)/6], [-10 10], 'LineStyle', '--', 'Color', 'r', 'LineWidth', lw);
        end
        set(gca, 'FontSize', f)
        box off;
        axis([1 size(calcium,2) ymin ymax]);

        subplot(2,2,2); %zoomed in version
        plot(1:size(calcium,2), nanmean(calcium,1), 'Color', map(s, :));
        hold on;
        %errorbars
        plot(1:size(calcium,2), nanmean(calcium,1)+nanstd(calcium,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);
        plot(1:size(calcium,2), nanmean(calcium,1)-nanstd(calcium,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);

        line([size(calcium,2)/2 size(calcium,2)/2], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);

        if strcmp(condition, 'optovin')
            line([size(calcium,2)/2+size(calcium,2)/6  size(calcium,2)/2 + size(calcium,2)/6], [-10 10], 'LineStyle', '--', 'Color', 'r', 'LineWidth', lw);
        end
        set(gca, 'FontSize', f)
        box off;

        %zoomed in axis
        axis([size(calcium,2)/2+20 size(calcium,2)/2+40 ymin ymax]);

        %tail
        subplot(2,2,3);
        plot(1:size(abs_tail,2), nanmean(abs_tail,1),  'Color', map(s, :));
        hold on;
        axis([1 size(abs_tail,2) 0 40]);
        line([size(abs_tail,2)/2 size(abs_tail,2)/2], [-100 100], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);
        set(gca, 'FontSize', f)
        box off;

        subplot(2,2,4);
        if strcmp(condition, 'water')
            if s==1
            else
                plot(1:size(calcium_norm,2), nanmean(calcium_norm,1), 'Color', map(s, :));
                hold on;
                %errorbars
                plot(1:size(calcium_norm,2), nanmean(calcium_norm,1)+nanstd(calcium_norm,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);
                plot(1:size(calcium_norm,2), nanmean(calcium_norm,1)-nanstd(calcium_norm,1)/sqrt(nneurons), 'Color', [0.5 0.5 0.5]);

                line([size(calcium,2)/2 size(calcium,2)/2], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);

                if strcmp(condition, 'optovin')
                    line([size(calcium,2)/2+size(calcium,2)/6  size(calcium,2)/2 + size(calcium,2)/6], [-10 10], 'LineStyle', '--', 'Color', 'r', 'LineWidth', lw);
                end
                set(gca, 'FontSize', f)
                box off;
                axis([1 size(calcium,2) ymin ymax]);
            end

        else

            plot(1:size(abs_tail,2), nanmean(abs_tail,1),  'Color', map(s, :));
            hold on;
            %axis([1 size(abs_tail,2) 0 40]);

            axis([size(abs_tail,2)/2+2*size(abs_tail,2)/18 size(abs_tail,2)/2+4*size(abs_tail,2)/18 0 40]);

            line([size(abs_tail,2)/2 size(abs_tail,2)/2], [-100 100], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);
            set(gca, 'FontSize', f)
            box off;

        end

        statstable{1,j} = calcium_forstats;

        figure(4)
        print(strcat('stimtriggeredaverage_', 'PO_', condition, stimtype), '-depsc2')
        figure(5)
        print(strcat('stimtriggeredaverage_', 'PT_', condition, stimtype), '-depsc2')
    end

end

% stats

pvals = [];
for i = 1:2
    [p, tbl, stats] = kruskalwallis(statstable{1,i}(:,:));
    c = multcompare(stats, 'Estimate','row')
    pvals(:,i) = c(:, 6);
end

%% plot kinematics per stim (using stim-triggered)
stimtype = 'pump'

%probably should use teensy setting above
stimchannel = {};
stimchannel{1} = 'control';
stimchannel{2} = 'odor2';

databyfish = {};

figure(5)
clf;
figure(6)
clf;

fnames = fieldnames(dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan));

% same for control and odor (experiment-specific)

count = 1;

for i = 1:8 %length(fnames)

    parameter = char(fnames(i));

    array = [];
    arraybyfish = [];

    for s = 1:length(stimchannel)
        stimchan = char(stimchannel{s});

        stimname = strcat('fishID', stimchan);
        fishID = dataarray.(condition).stim.(stimtype).stimtriggeredtail.(stimname);


        array1 = dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).(parameter);
        array(1:length(array1),s) = array1;


        for fish = 1:max(fishID)
            arraybyfish(fish, s) = nanmean(array1(find(fishID==fish)));
        end

        figure(5)
        subplot(3,4, count)
        Barandscatter_color(arraybyfish, [0 3]);
        ylabel(parameter);

    end

    [p2,h2] = signrank(arraybyfish(:,1), arraybyfish(:,2), 'tail', 'right')
    xlabel(num2str(p2));

    %ylabel(strcat('n=', num2str(length(arraybyfish))));


    databyfish{i}=arraybyfish;
    count = count+1;

end

figure(5)
print(strcat('tailkinematicsperstimperfish', condition, stimtype), '-depsc2')

count = 1;
for i = 9:length(fnames)

    parameter = char(fnames(i));

    array = [];
    arraybyfish = [];

    for s = 1:length(stimchannel)
        stimchan = char(stimchannel{s});
        array1 = dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).(parameter);
        array(1:length(array1),s) = array1;


        for fish = 1:max(fishID)
            arraybyfish(fish, s) = nanmean(array1(find(fishID==fish)));
        end

        figure(6)
        subplot(3,5, count)
        Barandscatter_color(arraybyfish, [0 3]);
        ylabel(parameter);

    end

    [p2,h2] = signrank(arraybyfish(:,1), arraybyfish(:,2), 'tail', 'right')
    xlabel(num2str(p2));

    %ylabel(strcat('n=', num2str(length(arraybyfish))));


    databyfish{i}=arraybyfish;
    count = count+1;

end

figure(6)
print(strcat('boutfreqperstimperfish', condition, stimtype), '-depsc2')


%% calcium sum analysis per neuron
figure(2);
clf;

for s = 1:length(stimchannel)

    figure(2)

    subplot(2,length(stimchannel),s);
    hist(calciumsum_post_byneuron(:,s)-calciumsum_pre_byneuron(:,s), [-0.6:0.05:0.6]);
    axis([-0.6 0.6 0 300])
    box off;

end

print(strcat('sumcalcium_histograms', condition, stimtype), '-depsc2')

[sorted, idx] = sort(calciumsum_post_byneuron(:,1));
calciumsum_sorted = [];

for n = 1:length(calciumsum_post_byneuron)
    calciumsum_sorted(n,:) = calciumsum_post_byneuron(idx(n), :);
end

lims = [-0.3 0.3]
figure(3);
imagesc(calciumsum_post_byneuron, lims);
box off;


%% dfoverf concat for each neuron
% concatenate trace_byneuron

dfoverf_concat = [];
m = 0;
for s = 1:length(trace_byneuron)

    array = trace_byneuron{s};

    dfoverf_concat(1:size(array,1), m+1:m+size(array,2)) = array;
    m = m+ size(array,2);
end

dfoverf = dfoverf_concat;

%% compare all signals to 1) baseline 2) water
%NOTE: need to run stim triggered average again with stim type 'pump'!!!!


neurontype = 'all'

%specific to OXT neurons
idx_PT = find(ypos >anteriorthreshold);

nstimchannels = length(fieldnames(dataarray.(condition).calcium.xcorr_stimreg));
ROIsizes =  dataarray.(condition).calcium.ROIsizes; %no particular reason why ROI sizes, just to get right length
neuronvector = 1:length(ROIsizes);


stimcatarray = {};
stimcatarray{1} = 'water'
stimcatarray{2} = firstodor;
stimcatarray{3} = secondodor;
stimcatarray{4} = thirdodor;

stimchannel = {};

if strcmp(celltype, 'OXTPO')
    stimchannel{1} = 'control';
    stimchannel{2} = 'odor1';
    stimchannel{3} = 'odor2';
    stimchannel{4} = 'odor3';
end

xcorr_stimreg = dataarray.(condition).calcium.xcorr_stimreg;
xcorr_motorstimreg = dataarray.(condition).calcium.xcorr_motorstimreg;

%2 x SEM is threshold
diffthreshold = 0.05;

% are we doing analysis based on xcorr or sum calcium?
useregressor = 0;

corr_thresholds = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.1];
diff_thresholds = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.05];

fractions = {};
ratios = {};
fractionspercorrcoeff = {};

for d = 1:length(diff_thresholds)

    activity = {}; %this is to determine if neurons are activated or suppressed BY STIMULUS (e.g.UV + Optovin or water flow or odor flow) = use regressor
    binary = {}; %this is to determine if neurons are activated or suppressed by ODOR, RELATIVE TO WATER - use difference
    %selectivity = {}; %only for Optovin condition. This determines if neurons are only activated by UV, or already activated by water flow

    diffthreshold = diff_thresholds(d);
    corrthreshold = corr_thresholds(d);

    for s = 1:nstimchannels

        cat = stimcatarray{s};

        for n = 1:length(neuronvector)

            idx_n = neuronvector(n);
            % relative to baseline. Get a choice to use regressor
            % or mean calcium

            if useregressor == 0

                %default
                activity{1}(n,s) = 0;
                activity{2}(n,s) = 0;

                if calciumsum_post_byneuron(idx_n,s) < calciumsum_pre_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,s)-calciumsum_pre_byneuron(idx_n,s))>= diffthreshold
                    activity{1}(n,s) = 1; %suppressed2
                elseif calciumsum_post_byneuron(idx_n,s) > calciumsum_pre_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,s)-calciumsum_pre_byneuron(idx_n,s))>= diffthreshold
                    activity{2}(n,s) = 1; %activated2
                elseif  isnan(calciumsum_post_byneuron(idx_n,s))
                    activity{1}(n,s) = nan;
                    activity{2}(n,s) = nan;
                end

            elseif useregressor == 1 %use correlation coefficient


                activity{1}(n,s) = 0;
                activity{2}(n,s) = 0;

                % activity per odor
                if xcorr_stimreg.(cat)(idx_n,1) < -corrthreshold %&& xcorr_stimreg_pvalues.(cat)(idx_n,1)< pcutoff
                    activity{1}(n,s) = 1; %suppressed2
                    activity{2}(n,s) = 0;
                elseif xcorr_stimreg.(cat)(idx_n,1) > corrthreshold %&& xcorr_stimreg_pvalues.(cat)(idx_n,1)< pcutoff
                    activity{2}(n,s) = 1; %activated2
                    activity{1}(n,s) = 0;
                elseif  isnan(xcorr_stimreg.(cat)(idx_n,1))
                    activity{1}(n,s) = nan;
                    activity{2}(n,s) = nan;

                end
            end


            % odor relative to water - use calcium sum difference ONLY

            if calciumsum_post_byneuron(idx_n,1) > calciumsum_post_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,1)-calciumsum_post_byneuron(idx_n,s))>= diffthreshold
                binary{1}(n,s) = 1; %suppressed
            elseif calciumsum_post_byneuron(idx_n,1) < calciumsum_post_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,1)-calciumsum_post_byneuron(idx_n,s))>= diffthreshold
                binary{2}(n,s) = 1; %activated
            elseif  isnan(calciumsum_post_byneuron(idx_n,s))
                binary{1}(n,s) = nan;
                binary{2}(n,s) = nan;
            else
                binary{1}(n,s) = 0;
                binary{2}(n,s) = 0;
            end

        end

        %activated/suppressed relative to water
        suppressed = binary{1}(:,s);
        activated = binary{2}(:,s);

        %activated/suppressed by odor
        suppressed2 = activity{1}(:,s);
        activated2 = activity{2}(:,s);

        %suppressed/activated relative to WATER
        fractionsuppressed = length(suppressed(suppressed==1))/length(suppressed(~isnan(suppressed)));
        fractionactivated = length(activated(activated==1))/length(activated(~isnan(activated)));

        %suppressed/activated by stimulus (i.e. water or odor)
        fractionsuppressed2 = length(suppressed2(suppressed2==1))/length(suppressed2(~isnan(suppressed2)));
        fractionactivated2 = length(activated2(activated2==1))/length(activated2(~isnan(activated2)));

        fraction_nochange =  length(find(suppressed==0 & activated==0))/length(suppressed);

        fraction_nochange2 =  length(find(suppressed2==0 & activated2==0))/length(suppressed2);


        fractions{s}(d,1) = fractionsuppressed; %odor relative to water
        fractions{s}(d,2) = fractionactivated;
        fractions{s}(d,3) = fraction_nochange;
        fractions{s}(d,4) = fractionsuppressed2; %relative to baseline
        fractions{s}(d,5) = fractionactivated2;
        fractions{s}(d,6) = fraction_nochange2;

    end

end

row = length(diff_thresholds);
figure(5);
clf
colormap(jet)

try
    b = bar([fractions{2}(row,1:3); fractions{3}(row, 1:3); fractions{4}(row, 1:3)], 'stacked');
catch
    b = bar([fractions{2}(row,1:3); fractions{3}(row, 1:3)], 'stacked');
end

b(1).FaceColor = 'b';
b(2).FaceColor = 'r';
b(3).FaceColor = 'k';
box off
print('fractions', '-depsc2')


figure(6)
clf;
map3 = [0 0 1; 1 0 0; 0 0 0];

if strcmp(condition, 'optovin')
    addnum = 6;
else
    addnum = 0;
end

for s = 2:length(stimchannel)
    subplot(3,4,s-1)
    for i = 1:3
        column = fractions{s}(:,i+addnum);
        plot(diff_thresholds(1:row-1), column(1:row-1), 'color', map3(i,:))
        hold on;
        box off;
        axis([0 diff_thresholds(end-1) 0 1]);
    end
end

print('fractions-allthresholds', '-depsc2')

%% statistics

close all;
binary_suppressed(binary_suppressed == -1) = 1;

for i = 1:length(binary_suppressed)

    if binary_suppressed(i,2) == 1 | binary_suppressed(i,3) == 1 |binary_suppressed(i,4) == 1
        binary_suppressed(i,1) = 1;
    end
end

suppressedbyodor = calciumsum_post_byneuron.*binary_suppressed;
suppressedbyodor(suppressedbyodor ==0) = nan;


for i = 1:length(binary_activated)

    if binary_activated(i,2) == 1 | binary_activated(i,3) == 1 | binary_activated(i,4) == 1
        binary_activated(i,1) = 1;
    end
end

activatedbyodor = calciumsum_post_byneuron.*binary_activated;
activatedbyodor(activatedbyodor ==0) = nan;

pvals2 = []

[p, tbl, stats] = kruskalwallis(suppressedbyodor(:,:));
c = multcompare(stats, 'Estimate','row')
pvals2(:,1) = c(:, 6);


[p, tbl, stats] = kruskalwallis(activatedbyodor(:,:));
c = multcompare(stats, 'Estimate','row')
pvals2(:,2) = c(:, 6);

samplesize = []

for i = 1:4

    column_suppressed = suppressedbyodor(:,i);
    samplesize(i, 1) = length(column_suppressed(~isnan(column_suppressed)));

    column_activated = activatedbyodor(:,i);
    samplesize(i, 2) = length(column_activated(~isnan(column_activated)));

end

%close all;

%% quantify overlap between odor responses (with varying intensities)

%first combine suppressed and enhanced neurons
activity_combined = [];
activity_suppressed = activity{1,1};
activity_suppressed(activity_suppressed == 1) = -1;
activity_activated =  activity{1,2};
activity_combined = activity_suppressed + activity_activated;

if strcmp(celltype, 'OXTPO')
    orderarray = [nan 1 2; 1 2 3; 1 2 4] %non-kin, then adult
elseif nstimchannels == 3
    orderarray = [1 2 3] %change if more odors
elseif nstimchannels == 4
    orderarray = [1 2 3 4]
elseif nstimchannels ==2
    orderarray = [1 2]
end

index = 1;
sortedxcorrarray = {};
sortedcalciumarray = {};
catnamesarray = {};

for m = 1:size(orderarray,1)

    catnames = {};
    order = orderarray(m, :);
    sortedxcorrarray{m} = [];
    sortedcalciumarray{m} = [];
    catnamesarray{m} = {};
    activity_odor = [];
    xcorr_odor = [];
    calciumsum_odor = [];

    if m == 1
        shortestcolumn = xcorr_stimreg.kin;
    elseif  m == 2
        shortestcolumn = xcorr_stimreg.nonkin;
    elseif m == 3
        shortestcolumn = xcorr_stimreg.adult;
    end

    neuronvector = find(~isnan(shortestcolumn));


    index = 1;
    for s = order(~isnan(order))

        cat = stimcatarray{s};

        %activity currently determined by whether regressor = 1 or 0
        %and the threshold

        activity_odor(1:length(neuronvector),index) = activity_combined(neuronvector,s); %combine activity for each odor including water

        %xcorr
        column1 = xcorr_stimreg.(cat);
        xcorr_odor(1:length(neuronvector),index) = column1(neuronvector,1);

        %calciumsum
        column2 = calciumsum_post_byneuron(:,s)-calciumsum_pre_byneuron(:,s);

        calciumsum_odor(1:length(neuronvector),index) = column2(neuronvector,1);

        calciumsum_odor_norm = calciumsum_odor-calciumsum_odor(:,1);


        catnames{1,index} = cat;
        index = index + 1;

    end

    startnum = 1;
    idxold = [];
    sortedxcorr = [];
    sortedcalcium = [];
    sortedcalcium_norm = [];

    %which to sort by

    if useregressor == 1
        sortarray = xcorr_odor;
    else
        sortarray = calciumsum_odor;
    end

    for s = 1:length(order)-1 %the last column is already sorted based on previous column
        % this is the column by which to sort by
        column = activity_odor(:,s);

        idxk = find(column==1);
        idxk2 = find(column==-1);
        idxk3 = find(column==0);

        %do not want to sort neurons that have already been sorted
        idx_unique = setdiff(idxk, idxold);
        idx_unique2 = setdiff(idxk2, idxold);
        idx_unique3 = setdiff(idxk3, idxold);

        [sorted, sorted_idx] = sortrows(sortarray(idx_unique, :), s, 'descend');
        [sorted2, sorted_idx2] = sortrows(sortarray(idx_unique2, :), s, 'descend');
        [sorted3, sorted_idx3] = sortrows(sortarray(idx_unique3, :), s, 'descend');

        if useregressor == 1

            sortedxcorr(startnum:startnum+length(idx_unique)-1, 1:size(activity_odor,2)) = sorted;
            sortedxcorr(startnum+length(idx_unique)+1:startnum+length(idx_unique)+length(idx_unique2),1:size(activity_odor,2)) = sorted2;

            tmparray = calciumsum_odor(idx_unique, :);
            tmparray2 = calciumsum_odor(idx_unique2, :);
            tmparray3 = calciumsum_odor(idx_unique3, :);

            sortedcalcium(startnum:startnum+length(idx_unique)-1, 1:size(activity_odor,2)) = tmparray(sorted_idx, :);
            sortedcalcium(startnum+length(idx_unique):startnum+length(idx_unique)+length(idx_unique2)-1,1:size(activity_odor,2)) = tmparray2(sorted_idx2, :)

            % but you are sorting ALL the columns ultimately
            if s == length(order)-1
                sortedxcorr(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = sorted3;
            end

            if s == length(order)-1
               sortedcalcium(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = tmparray3(sorted_idx3, :)
            end


        else
            sortedcalcium(startnum:startnum+length(idx_unique)-1, 1:size(activity_odor,2)) = sorted;
            sortedcalcium(startnum+length(idx_unique):startnum+length(idx_unique)+length(idx_unique2)-1,1:size(activity_odor,2)) = sorted2;

            tmparray = xcorr_odor(idx_unique, :);
            tmparray2 = xcorr_odor(idx_unique2, :);
            tmparray3 = xcorr_odor(idx_unique3, :);

            sortedxcorr(startnum:startnum+length(idx_unique)-1, 1:size(activity_odor,2)) = tmparray(sorted_idx, :);
            sortedxcorr(startnum+length(idx_unique)+1:startnum+length(idx_unique)+length(idx_unique2),1:size(activity_odor,2)) = tmparray2(sorted_idx2, :);

            % but you are sorting ALL the columns ultimately
            if s == length(order)-1
                sortedxcorr(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = tmparray3(sorted_idx3, :)
            end

            if s == length(order)-1
                sortedcalcium(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = sorted3;
            end

        end


        %idxold = vertcat(idxold, idx_unique, idx_unique2, idx_unique3);
        idxold = vertcat(idxold, idx_unique, idx_unique2); %idx_unique3 is the remaining
        startnum = length(idxold) + 1;
    end

    catnamesarray{m} = catnames;
    sortedxcorrarray{m} = sortedxcorr;
    sortedcalciumarray{m} = sortedcalcium;

    %calcium sum normalized to water
    sortedcalciumarray_norm{m} = sortrows(calciumsum_odor_norm(:, 2:end));

end

%you can use the sorted xcorr or sum calcium array from now on

% plot correlation maps

%to generate a color map red to blue
if strcmp(celltype,'OB')
    corrlimits = [-0.6 0.6];
else
    corrlimits = [-0.3 0.3];
    corrlimits2 = [-0.3 0.3];
end

cmapred = [];
cmapblue = [];
nc_red =  45;
nc_blue = 45;
for c = 1:nc_red
    cmapred(c,1) = (c-1)*1/nc_red;
    cmapred(c,2) = (c-1)*1/nc_red;
    cmapred(c,3) = 1;
end

for c = 1:nc_blue
    cmapblue(c,1) = 1;
    cmapblue(c,2) = 1-c*1/nc_blue;
    cmapblue(c,3) = 1-c*1/nc_blue;
end

cmapwhite = ones(10,3);
cmap_fine = vertcat(cmapred, cmapwhite, cmapblue)

cmap2 = [];
cmap2(11,:) = [1 0 0];
cmap2(10,:) = [1 0.3 0.3];
cmap2(9,:) = [1 0.6 0.6];
cmap2(8,:) = [1 0.9 0.9];
cmap2(7,:) = [1 1 1];
cmap2(6,:) = [1 1 1];
cmap2(5,:) = [1 1 1];
cmap2(4,:) = [0.9 0.9 1];
cmap2(3,:) = [0.6 0.6 1];
cmap2(2,:) = [0.3 0.3 1];
cmap2(1,:) = [0 0 1];

figure(15)
clf
figure(16)
clf
figure(17)
clf

for i = 1:length(sortedxcorrarray)

    figure(15)
    subplot(1,3,i)
    colormap(cmap_fine)

    imagesc(sortedxcorrarray{i}, corrlimits)
    xlabel(catnamesarray{i});
    ylabel(num2str(length(sortedxcorrarray{i})))
    box off;
    print(strcat('clusteractivity_xcorr_', neurontype), '-deps')

    figure(16)
    subplot(1,3,i)
    colormap(cmap_fine)
    imagesc(sortedcalciumarray{i}, corrlimits2)
    xlabel(catnamesarray{i});
    ylabel(num2str(length(sortedcalciumarray{i})))
    box off;
    print(strcat('clusteractivity_calciumsum_', neurontype), '-deps')


    figure(17)
    subplot(1,3,i)
    colormap(cmap_fine)
    imagesc(sortedcalciumarray_norm{i}, corrlimits2)
    xlabel(catnamesarray{i});
    ylabel(num2str(length(sortedcalciumarray{i})))
    box off;
    print(strcat('clusteractivity_calciumsum_normalizedtowater', neurontype), '-deps')

end

%% quantify overlap
%cross compare active neurons

crosscomparisons.(neurontype).overlap = struct;
crosscomparisons.(neurontype).exclusive = struct;
crosscomparisons.(neurontype).nneurons = size(activity_odor,1);
neurontype = 'all'
stimlist = [1: size(activity_odor, 2)];

for s = 1:size(activity_odor, 2)

    cat = stimcatarray{s};
    column = activity_odor(:,s);
    otherstims = setdiff(stimlist,s);
   
    for k = otherstims
        othercat = stimcatarray{k};
        othercolumn = activity_odor(:,k);

        %A intersect B
        idx_intersect = find(column == 1 & othercolumn == 1);
        catname = strcat(cat, 'and', othercat);
        crosscomparisons.(neurontype).overlap.(catname) = idx_intersect;
        %A not B

        idx_exclusive = find(column == 1 & othercolumn ~= 1);
        catname = strcat(cat, 'not', othercat);
        crosscomparisons.(neurontype).exclusive.(catname) = idx_exclusive;

    end
end

%% Unsupervised k means clustering of xcorr values
% use sorted values
% kmeans
%
% k = 9;
% clusterarrays = {};
% idxes = {};
%
% for i = 1:length(sortedxcorrarray)
% clusterarray = sortedxcorrarray{i};
% idx = kmeans(clusterarray, k, 'Display', 'iter', 'OnlinePhase', 'on', 'Replicates', 5);
% clusterarrays{i} = clusterarray;
% idxes{i} = idx;
% end

%% to save k means clusters if necessary
c = clock;
time = strcat(num2str(c(4)), num2str(c(5)));
save(sprintf('%s-%s-%s-%s%s', 'analyzeddata', celltype, date, time, '.mat'));
%% process k means clusters
%to generate a color map red to blue

if strcmp(celltype, 'OB')
    corrlimits = [-0.6 0.6];
else
    corrlimits = [-0.6 0.6];
end

corrlimits2 = [-0.3 0.3]

cmapred = [];
cmapblue = [];
nc_red =  45;
nc_blue = 45;
for c = 1:nc_red
    cmapred(c,1) = (c-1)*1/nc_red;
    cmapred(c,2) = (c-1)*1/nc_red;
    cmapred(c,3) = 1;
end

for c = 1:nc_blue
    cmapblue(c,1) = 1;
    cmapblue(c,2) = 1-c*1/nc_blue;
    cmapblue(c,3) = 1-c*1/nc_blue;
end

cmapwhite = ones(10,3);
cmap_fine = vertcat(cmapred, cmapwhite, cmapblue)

dfoverf = dfoverf_concat;
rearrangearray = [];

%rearrange = [5 1 8 7 6 4 2 3] %GABAoxy
if strcmp(celltype, 'OB')
    rearrange = 1:k %OB
    rearrange = [10 9 3 1 12 6 2 4 8 7 11 5];
elseif strcmp(celltype, 'GABAoxy')
    rearrange = 1:k;
    rearrange = [10 6 1 11 4 12 7 5 2 3 8 9];
elseif strcmp(celltype, 'GABAonly')
    rearrange = 1:k;
    rearrange = [3 6 7 9 12 10 4 5 1 2 8 11];
elseif strcmp(celltype, 'OXTPO')
    rearrangearray(1,:) = [1 9 6 5 8 2 3 4 7];
    rearrangearray(2,:) = [2 6 8 7 9 1 3 4 5];
    rearrangearray(3,:) = [1 9 4 6 7 2 8 5 3];
end

fishID = dataarray.(condition).calcium.fishID;

sorted_diffs = {};
sorted_calciumsums = {};
sortedclusters = {};
clustersizesarray = {};
sorted_neuronID = {};
sorted_idxes = {};

for i = 1:length(sortedxcorrarray)
    rearrange = rearrangearray(i,:);
    clusterarray = clusterarrays{i};
    calciumsumarray = sortedcalciumarray{i};
    idx = idxes{i};

    sortedcalciumsum = [];
    sortedcluster = [];
    sorteddfoverf = [];
    clustersizes = [];
    clustermeans = [];
    sorted_neuronID = [] ;
    sorted_idx = [] ;
    sorted_fishID = [];


    for s = rearrange %needs to be a single row not column
        % this is the column by which to sort by

        %but you are sorting ALL the columns ultimately
        sortedcluster(length(sortedcluster)+1: length(sortedcluster) + length(find((idx==s))), 1:size(clusterarray,2)) = clusterarray(idx==s, :);
        sortedcalciumsum(length(sortedcalciumsum)+1: length(sortedcalciumsum) + length(find((idx==s))), 1:size(calciumsumarray,2)) = calciumsumarray(idx==s, :);
        clustersizes(s,1) =  length(find((idx==s)));
        sorted_neuronID(length(sorted_neuronID)+1: length(sorted_neuronID) + length(find((idx==s))), 1) = find(idx==s);
        sorted_idx(length(sorted_idx)+1: length(sorted_idx) + clustersizes(s), 1) = repmat(s, [clustersizes(s), 1]);

        sorteddfoverf(size(sorteddfoverf,1)+1: size(sorteddfoverf,1) + length(find((idx==s))), 1:size(dfoverf,2)) = dfoverf(idx==s, :);
        clustermeans(s, 1:size(dfoverf,2)) = mean(dfoverf(idx==s, :),1);
        sorted_fishID(length(sorted_fishID)+1: length(sorted_fishID) + length(find((idx==s))), 1) = fishID(find(idx==s));
    end

    sortedclusters{i} = sortedcluster;
    sortedcalciumsums{i} = sortedcalciumsum;
    clustersizesarray{i} = clustersizes;
    sorted_neuronIDs{i} = sorted_neuronID;
    sorted_idxes{i} = sorted_idx;
    sorted_diffs{i} = sortedcalciumsum-sortedcalciumsum(:,1);

end

figure(16)
clf;
figure(17)
clf;
figure(18);
clf;

for i = 1:length(sortedxcorrarray)
    figure(15+i);

    catnames = catnamesarray{i};
    try
        type = strcat(catnames(1), catnames(2), catnames(3))
    catch
        type = strcat(catnames(1), catnames(2))
    end


    subplot(1,3,1);
    colormap(cmap_fine)
    sortedcluster = sortedclusters{i};

    imagesc(sortedcluster, corrlimits2)
    xlabel(type)
    box off;

    subplot(1,3,2);
    colormap(cmap_fine)
    calciumsumcluster = sortedcalciumsums{i};
    imagesc(calciumsumcluster, corrlimits2)
    box off;

    subplot(1,3,3);
    colormap(cmap_fine)
    calciumdiffcluster = sorted_diffs{i};

    if i ==1
        imagesc(calciumdiffcluster(:, 2), corrlimits2)
    else
        imagesc(calciumdiffcluster(:, end-1:end), corrlimits2)
    end

    box off;
    colorbar

    colorarray = hsv(k);
    count = 1;
    nneurons = length(sortedcluster);
    clustersizes = clustersizesarray{i};

    rearrange = rearrangearray(i,:);

    for s = rearrange
        color = colorarray(s, :);
        text(size(sortedcluster,2)+0.5, count, strcat(num2str(s), ': n = ', num2str(clustersizes(s)), ' (', num2str(clustersizes(s)/nneurons*100), '%)'), 'Color', color)
        count = count+ clustersizes(s);
        hold on;
        ylim([0 nneurons]);
        set(gca, 'YDir', 'reverse')
    end


    print(strcat('kmeansclustersorted_', char(type)), '-depsc2', '-painters')

end

dataarray.(condition).calcium.kmeansxcorr_idx = idxes;
dataarray.(condition).calcium.kmeansxcorr_clustersizes = clustersizesarray;
dataarray.(condition).calcium.kmeansxcorr_sorted = sortedclusters;
dataarray.(condition).calcium.kmeansxcorr_sortedidx = sorted_idxes;
dataarray.(condition).calcium.kmeansxcorr_dfoverf = sorteddfoverf;
dataarray.(condition).calcium.kmeansxcorr_fishIDs = sorted_fishID;
dataarray.(condition).calcium.kmeansxcorr_neuronIDs = sorted_neuronIDs;

figure(17)
clf;
subplot(1,2,1);
clims = [-0.5 2];
imagesc(sorteddfoverf, clims)
box off;
count = 1;
n = k;
for s = rearrange
    subplot(1,2,1);
    color = colorarray(s, :);
    text(size(sorteddfoverf,2) + 300, count, num2str(s), 'Color', color)
    count = count+ clustersizes(s);
    hold on;
    ylim([0 nneurons]);
    set(gca, 'YDir', 'reverse')
    colorbar;

    subplot(1,2,2)
    plot(1:size(sorteddfoverf,2), clustermeans(s,:)+n, 'Color', color)
    hold on;
    text(size(sorteddfoverf,2) + 50, n, num2str(s), 'Color', color)
    box off;
    n = n-1;
    xlim([0  size(sorteddfoverf,2)])
    ylim([0 k+1])
end

print(strcat('kmeansclusterraster_', neurontype), '-depsc2','-painters')

%by fish
dataarray.(condition).calcium.kmeansxcorr_idx = idx;
dataarray.(condition).calcium.kmeansxcorr_clustersizes = clustersizes;
dataarray.(condition).calcium.kmeansxcorr_sorted = sortedcluster;
dataarray.(condition).calcium.kmeansxcorr_sortedidx = sorted_idx;
dataarray.(condition).calcium.kmeansxcorr_dfoverf = sorteddfoverf;
dataarray.(condition).calcium.kmeansxcorr_fishIDs = sorted_fishID;
dataarray.(condition).calcium.kmeansxcorr_neuronIDs = sorted_neuronID;
%save(sprintf('%s-%s%s', 'analyzeddata', celltype, '.mat'));

%% split clusters per fish
clustersizes_perfish = [];
sortedidx_perfish = [];

figure(18)
clf

for fish = 1:max(fishID)
    sortedcluster_perfish = sortedcluster(sorted_fishID == fish, :);
    sortedidx_perfish = sorted_idx(sorted_fishID == fish, :);

    figure(18)
    colormap(cmap_fine)
    subplot(ceil(max(fishID)/2), 2, fish)
    imagesc(sortedcluster_perfish, corrlimits)

    count = 1;
    for s = rearrange
        color = colorarray(s, :);
        clustersizes_perfish(s, fish) = length(find(sortedidx_perfish == s));
        text(size(sortedcluster,2)+0.7, count, strcat(num2str(s)), 'Color', color);
        count = count+ clustersizes_perfish(s, fish)+2;
        hold on;
        %ylim([0 nneurons]);
        set(gca, 'YDir', 'reverse')
    end


end

dataarray.(condition).calcium.kmeansxcorr_clustersizesperfish = clustersizes_perfish;

print(strcat('kmeansclusterraster_perfish', neurontype), '-depsc2')
%% GMM
k = 1:10;
nK = numel(k);
Sigma = {'diagonal','full'};
nSigma = numel(Sigma);
SharedCovariance = {true,false};
SCtext = {'true','false'};
nSC = numel(SharedCovariance);
RegularizationValue = 0.01;
options = statset('MaxIter',10000);

m = cell(nK,nSigma,nSC);
aic = zeros(nK,nSigma,nSC);
bic = zeros(nK,nSigma,nSC);
converged = false(nK,nSigma,nSC);

% Fit all models
for m = 1:nSC
    for j = 1:nSigma
        for i = 1:nK
            gm{i,j,m} = fitgmdist(clusterarray,k(i),...
                'CovType',Sigma{j},...
                'SharedCov',SharedCovariance{m},...
                'Options', options)
            aic(i,j,m) = gm{i,j,m}.AIC;
            bic(i,j,m) = gm{i,j,m}.BIC;
            converged(i,j,m) = gm{i,j,m}.Converged;
        end
    end
end

allConverge = (sum(converged(:)) == nK*nSigma*nSC)
figure;
bar(reshape(aic,nK,nSigma*nSC));
title('AIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
xlabel('$k$','Interpreter','Latex');
ylabel('AIC');
legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
    'Full-unshared'});
box off;

figure;
bar(reshape(bic,nK,nSigma*nSC));
title('BIC For Various $k$ and $\Sigma$ Choices','Interpreter','latex');
xlabel('$c$','Interpreter','Latex');
ylabel('BIC');
legend({'Diagonal-shared','Full-shared','Diagonal-unshared',...
    'Full-unshared'});
box off;

%% plot stim triggered averaged for suppressed vs activated neurons

%set threshold here
%can also be a difference threshold
diff_threshold = 0.05;
corr_threshold = 0.2;
useregressor = 0;
figure(7)
clf

stimtriggereddata = struct;

figure(8)
clf

figure(9)
clf

figure(10)
clf

if strcmp(condition, 'optovin')
    ymin = -1;
    ymax =4;
else
    ymin = -0.7;
    ymax = 0.7;
end

for s = 1:length(stimchannel)

    stimchan = char(stimchannel{s});

    if strcmp(condition, 'water')
        calcium_raw = trace_byneuron{s}-trace_byneuron{1}; %normalized
    elseif strcmp(condition, 'optovin')
        calcium_raw = trace_byneuron{s}; %not normalized
    end

    totaltimepoints = size(calcium_raw,2);
    baselineend = size(calcium_raw,2)/2;

    if s == 1
        suppressed = binary{1}(:,2);
        activated = binary{2}(:,2);
    else
        suppressed = binary{1}(:,s);
        activated = binary{2}(:,s);
    end

    %relative to water

    calciumraw_suppressed = calcium_raw(find(suppressed==1),:);
    calciumraw_activated = calcium_raw(find(activated==1),:);
    calciumraw_nochange = calcium_raw(find(suppressed==0 & activated==0),:);

    calciumraw_suppressed(calciumraw_suppressed==0) = nan;
    calcium_suppressed = calciumraw_suppressed -nanmean(calciumraw_suppressed(:, 1:baselineend), 2);

    calciumraw_activated(calciumraw_activated==0) = nan;
    calcium_activated = calciumraw_activated -nanmean(calciumraw_activated(:, 1:baselineend), 2);

    calciumraw_nochange(calciumraw_nochange==0) = nan;
    calcium_nochange = calciumraw_nochange-nanmean(calciumraw_nochange(:, 1:baselineend), 2);

    %relative to baseline
    suppressed2 = activity{1}(:,1);
    activated2 = activity{2}(:,1);

    %this is to use regressor instead of difference threshold
    if useregressor ==1
        activated2 = zeros(length(activated), 1);
        suppressed2 = zeros(length(suppressed), 1);
        xcorr_stimreg = dataarray.(condition).calcium.xcorr_stimreg;
        xcorr_motorstimreg = dataarray.(condition).calcium.xcorr_motorstimreg;
        xcorr_motorsponreg = dataarray.(condition).calcium.xcorr_motorsponreg;

        map2 = [204/255 0 102/255; 1  102/255 178/255; 1 153/255 51/255; 0 0 0];

        slidingstats = [];
        slidingstats2 = [];
        corr_thresholds = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6];
        corr_thresholds2 = [0 -0.05 -0.1 -0.15 -0.2 -0.25 -0.3 -0.35 -0.4 -0.45 -0.5 -0.55 -0.6];

        for k = 1:length(corr_thresholds)
            idx_activated = find(xcorr_stimreg>=corr_thresholds(k));
            slidingstats(k,1) =  length(idx_activated)/nneurons*100;
            idx_activated = find(xcorr_motorstimreg>=corr_thresholds(k));
            slidingstats(k,2) =  length(idx_activated)/nneurons*100;
            idx_activated = find(xcorr_motorsponreg>=corr_thresholds(k));
            slidingstats(k,3) =  length(idx_activated)/nneurons*100;

            idx_suppressed = find(xcorr_stimreg<=corr_thresholds2(k));
            slidingstats2(k,1) =  length(idx_suppressed)/nneurons*100;
            idx_suppressed = find(xcorr_motorstimreg<=corr_thresholds2(k));
            slidingstats2(k,2) =  length(idx_suppressed)/nneurons*100;
            idx_suppressed = find(xcorr_motorsponreg<=corr_thresholds2(k));
            slidingstats2(k,3) =  length(idx_suppressed)/nneurons*100;

        end

        figure(10)
        clf;
        for ss = 1:3
            figure(3);
            subplot(1,2,2);
            plot(slidingstats(:,ss), 'Color', map2(ss,:))
            axis([1 13 0 100]);
            set(gca, 'XTick', [2 4 6 8 10 12]);
            set(gca, 'XTickLabel', [0.05 0.15 0.25 0.35 0.45 0.55]);
            hold on;
            box off;

            subplot(1,2,1);
            plot(slidingstats2(:,ss), 'Color', map2(ss,:))
            axis([1 13 0 100]);
            set(gca, 'Xdir', 'reverse')
            set(gca, 'YAxisLocation', 'right')
            set(gca, 'XTick', [2 4 6 8 10 12]);
            set(gca, 'XTickLabel', [-0.05 -0.15 -0.25 -0.35 -0.45 -0.55]);
            hold on;
            box off;
        end

        print('slidingstats', '-depsc2')
    end

    %%%%%% for a specific corr coeff threshold %%%%%%%

    activated2 = zeros(length(activated), 1);
    suppressed2 = zeros(length(suppressed), 1);

    if useregressor == 1
        idx_activated = find(xcorr_stimreg>=corr_threshold  | xcorr_motorstimreg>=corr_threshold);
        activated2(idx_activated) = 1;
        idx_suppressed = find(xcorr_stimreg<=-corr_threshold  | xcorr_motorstimreg<=-corr_threshold);
        suppressed2(idx_suppressed) = 1;

    else

        activated2 = zeros(length(activated), 1);
        suppressed2 = zeros(length(suppressed), 1);

        for n = 1:length(calciumsum_post_byneuron)
            if calciumsum_post_byneuron(n,1) < calciumsum_pre_byneuron(n,1) && abs(calciumsum_post_byneuron(n,1)-calciumsum_pre_byneuron(n,1))>= corr_threshold
                suppressed2(n) = 1; %suppressed2
            elseif calciumsum_post_byneuron(n,1) > calciumsum_pre_byneuron(n,1) && abs(calciumsum_post_byneuron(n,1)-calciumsum_pre_byneuron(n,1))>= corr_threshold
                activated2(n) = 1; %activated2
            else
                suppressed2(n) = 0;
                activated2(n) = 0;
            end
        end
    end


    if length(activated2)<length(activated)
        diff = length(activated)-length(activated2);
        activated2(end+1:end+diff) = 0;
        activity{2}(end+1:end+diff,:) = 0;
    elseif length(suppressed2)<length(suppressed)
        diff = length(suppressed)-length(suppressed2);
        suppressed2(end+1:end+diff)=0;
        activity{1}(end+1:end+diff,:) = 0;
    end


    calciumraw_suppressed2 = calcium_raw(find(suppressed2==1),:);
    calciumraw_activated2 = calcium_raw(find(activated2==1),:);

    calciumraw_suppressed2(calciumraw_suppressed2==0) = nan;
    calcium_suppressed2 = calciumraw_suppressed2 -nanmean(calciumraw_suppressed2(:, 1:baselineend), 2);

    calciumraw_activated2(calciumraw_activated2==0) = nan;
    calcium_activated2 = calciumraw_activated2 -nanmean(calciumraw_activated2(:, 1:baselineend), 2);

    figure(7);
    title('All OXT neurons');
    subplot(2,2,1);
    plot(1:size(calcium_suppressed,2), nanmean(calcium_suppressed,1), 'Color', map(s, :));
    hold on;
    %errorbars
    plot(1:size(calcium_suppressed,2), nanmean(calcium_suppressed,1)+nanstd(calcium_suppressed,1)/sqrt(size(calcium_suppressed,1)), 'Color', [0.5 0.5 0.5]);
    plot(1:size(calcium_suppressed,2), nanmean(calcium_suppressed,1)-nanstd(calcium_suppressed,1)/sqrt(size(calcium_suppressed,1)), 'Color', [0.5 0.5 0.5]);

    line([size(calcium_suppressed,2)/2 size(calcium_suppressed,2)/2], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);

    if strcmp(condition, 'optovin')
        line([size(calcium_suppressed,2)/2+size(calcium_suppressed,2)/6  size(calcium_suppressed,2)/2 + size(calcium_suppressed,2)/6], [-10 10], 'LineStyle', '--', 'Color', 'r', 'LineWidth', lw);
    end
    set(gca, 'FontSize', f)
    box off;
    axis([1 size(calcium_suppressed,2) ymin ymax]);
    ylabel('Suppressedbyodor');


    subplot(2,2,2);
    plot(1:size(calcium_activated,2), nanmean(calcium_activated,1), 'Color', map(s, :));
    hold on;
    %errorbars
    plot(1:size(calcium_activated,2), nanmean(calcium_activated,1)+nanstd(calcium_activated,1)/sqrt(size(calcium_activated,1)), 'Color', [0.5 0.5 0.5]);
    plot(1:size(calcium_activated,2), nanmean(calcium_activated,1)-nanstd(calcium_activated,1)/sqrt(size(calcium_activated,1)), 'Color', [0.5 0.5 0.5]);

    line([size(calcium_activated,2)/2 size(calcium_activated,2)/2], [-10 10], 'LineStyle', '--', 'Color', 'k', 'LineWidth', lw);

    if strcmp(condition, 'optovin')
        line([size(calcium_activated,2)/2+size(calcium_activated,2)/6  size(calcium_activated,2)/2 + size(calcium_activated,2)/6], [-10 10], 'LineStyle', '--', 'Color', 'r', 'LineWidth', lw);
    end
    set(gca, 'FontSize', f)
    box off;
    axis([1 size(calcium_activated,2) ymin ymax]);
    ylabel('Enhancedbyodor');

    stimtriggereddata.(stimchan)= struct;
    stimtriggereddata.(stimchan).activated = calcium_activated;
    stimtriggereddata.(stimchan).suppressed = calcium_suppressed;

    % scatterplots

    g1 = calciumsum_post_byneuron(:,1);
    g2 = calciumsum_post_byneuron(:,s);

    figure(11)
    subplot(3,4,s)

    if strcmp(condition, 'optovin')
        scatter(g2(activated2==1), g1(activated2==1), 'k.');
        hold on;
        scatter(g2(suppressed3==1),g1(suppressed3==1), 'b.');
        scatter(g2(activated3==1),g1(activated3==1), 'r.');
    else
        scatter(g2, g1, 'k.');
        hold on;
        scatter(g2(suppressed==1),g1(suppressed==1), 'b.');
        scatter(g2(activated==1),g1(activated==1), 'r.');
    end

    line([-0.5 1],[-0.5 1], 'color', 'k')
    axis([-0.5 1 -0.5 1])

end

figure(7)
print(strcat('stimtriggeredaverage-classified', condition, stimtype, 'threshold', num2str(corr_threshold*100)), '-depsc2')

figure(8)
print(strcat('stimtriggeredaverage-classified2', condition, stimtype, 'threshold', num2str(corr_threshold*100)), '-depsc2')

figure(9)
print(strcat('stimtriggeredaverage-classified-zoomin', condition, stimtype, 'threshold', num2str(corr_threshold*100)), '-depsc2')

figure(10)
print(strcat('calciumsum_scatterplot', condition, stimtype), '-depsc2')


%% quantify overlap between odor responses

%if responds to kin, does it also respond to other odors? If so which
%direction?
binary_combined = [];
binary_suppressed = binary{1,1};
binary_suppressed(binary_suppressed == 1) = -1;
binary_activated =  binary{1,2};
binary_combined = binary_suppressed + binary_activated;

% first, split binary into 2 groups
idx_odor3 = find(~isnan(binary_combined(:,3)));
idx_odor4 = find(~isnan(binary_combined(:,4)));

binary_odor3 = [];
binary_odor4 = [];
binary_odor3(:,1) = binary_combined(idx_odor3,2);
binary_odor3(:,2) = binary_combined(idx_odor3,3);

binary_odor4(:,1) = binary_combined(idx_odor4,2);
binary_odor4(:,2) = binary_combined(idx_odor4,4);

binary_perodor = {};
binary_perodor{1,1} = binary_odor3;
binary_perodor{1,2} = binary_odor4;

sortedbinary_perodor = {};

for s = 1:length(binary_perodor)

    startnum = 1;
    idxold = [];
    sortedbinary = [];

    for k = 1:2

        column = binary_perodor{1,s}(:,k);
        idxk = find(column==-1);
        idxk2 = find(column==0);
        idxk3 = find(column==1);
        %do not want to sort fish that have already been sorted
        idx_unique = setdiff(idxk, idxold);
        idx_unique2 = setdiff(idxk2, idxold);
        idx_unique3 = setdiff(idxk3, idxold);

        sortedbinary(startnum:startnum+length(idx_unique)-1,1:size(binary_perodor{1,s},2)) =sort(binary_perodor{1,s}(idx_unique,:));
        sortedbinary(startnum+length(idx_unique):startnum+length(idx_unique)+length(idx_unique2)-1,1:size(binary_perodor{1,s},2)) = sort(binary_perodor{1,s}(idx_unique2,:));
        sortedbinary(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(binary_perodor{1,s},2)) = sort(binary_perodor{1,s}(idx_unique3,:));

        idxold = vertcat(idxold, idx_unique, idx_unique2, idx_unique3);
        startnum = length(idxold) + 1;
    end

    sortedbinary_perodor{1,s} = sortedbinary;
end

% plot correlation maps
cmap2(11,:) = [1 0 0];
cmap2(10,:) = [1 0.3 0.3];
cmap2(9,:) = [1 0.6 0.6];
cmap2(8,:) = [1 0.9 0.9];
cmap2(7,:) = [1 1 1];
cmap2(6,:) = [1 1 1];
cmap2(5,:) = [1 1 1];
cmap2(4,:) = [0.9 0.9 1];
cmap2(3,:) = [0.6 0.6 1];
cmap2(2,:) = [0.3 0.3 1];
cmap2(1,:) = [0 0 1];

range = [-1 1]
colormap(cmap2)
figure(15)
subplot(1,3,1);
imagesc(sortedbinary_perodor{1,1}, range)
box off;
%print('clusterbinaryodor3', '-depsc2', '-painters')

hold on;
colormap(cmap2)
subplot(1,3,2)
imagesc(sortedbinary_perodor{1,2}, range)
box off;
print('clusterbinaryodor3and4', '-depsc2', '-painters')


%% overlap between kin and other cues

% -/- -/+
% +/- +/+  (suppressed or activated or none)

odor3overlap = zeros(2,2);
odor4overlap = zeros(2,2);

odor3overlap(1,1) = length(find(binary_odor3(:,1)==-1 & binary_odor3(:,2)==-1))/length(find(binary_odor3(:,1)==-1));
odor3overlap(1,2) = length(find(binary_odor3(:,1)==-1 & binary_odor3(:,2)==1))/length(find(binary_odor3(:,1)==-1));
odor3overlap(2,1) = length(find(binary_odor3(:,1)==1 & binary_odor3(:,2)==-1))/length(find(binary_odor3(:,1)==1));
odor3overlap(2,2) = length(find(binary_odor3(:,1)==1 & binary_odor3(:,2)==1))/length(find(binary_odor3(:,1)==1));

odor4overlap(1,1) = length(find(binary_odor4(:,1)==-1 & binary_odor4(:,2)==-1))/length(find(binary_odor4(:,1)==-1));
odor4overlap(1,2) = length(find(binary_odor4(:,1)==-1 & binary_odor4(:,2)==1))/length(find(binary_odor4(:,1)==-1));
odor4overlap(2,1) = length(find(binary_odor4(:,1)==1 & binary_odor4(:,2)==-1))/length(find(binary_odor4(:,1)==1));
odor4overlap(2,2) = length(find(binary_odor4(:,1)==1 & binary_odor4(:,2)==1))/length(find(binary_odor4(:,1)==1));

%% compare maximum signal post-stim to pre-stim
maxdiffs = [];
maxdiffs_sum = [];
maxdiffs_activated = [];
maxdiffs_activated2 = [];
sumdiffs = [];
sumdiffs_activated = [];
sumdiffs_activated2 = [];

if strcmp(condition, 'optovin')
    for s = 1:length(stimchannel)

        maxdiffs(:,s) = calciummax_post_byneuron(:,s) - calciummax_pre_byneuron(:,s);
        maxdiffs_sum(:,s) = calciumsum_post_byneuron(:,s) - calciumsum_pre_byneuron(:,s);
    end

    usesum = 0;

    if strcmp(condition, 'optovin')

        for s = 1:length(stimchannel)

            if usesum == 1
                column = maxdiffs_sum(:,s);
            else
            end
            column = maxdiffs(:,s);
            maxdiffs_activated(:,s) = column(suppressed3==1); %these are activated neurons that are suppressed
            maxdiffs_activated2(:,s) = column(activated3==1);%these are activated neurons that are activated
        end
    end


    thetathreshold = 0;
    idx_enhanced = find(maxdiffs_activated(:,2)>maxdiffs_activated(:,1) & abs(maxdiffs_activated(:,2)-maxdiffs_activated(:,1))>=thetathreshold);
    idx_suppressed = find(maxdiffs_activated(:,2)<maxdiffs_activated(:,1) & abs(maxdiffs_activated(:,2)-maxdiffs_activated(:,1))>=thetathreshold);

    idx_enhanced2 = find(maxdiffs_activated2(:,2)>maxdiffs_activated2(:,1) & abs(maxdiffs_activated2(:,2)-maxdiffs_activated2(:,1))>=thetathreshold);
    idx_suppressed2 = find(maxdiffs_activated2(:,2)<maxdiffs_activated2(:,1) & abs(maxdiffs_activated2(:,2)-maxdiffs_activated2(:,1))>=thetathreshold);


    figure(15)
    clf;
    subplot(3,4,1);
    scatter(maxdiffs_activated(:,2), maxdiffs_activated(:,1), 'k.');
    hold on;
    %scatter(maxdiffs_activated(idx_enhanced,2), maxdiffs_activated(idx_enhanced,1),'r.');
    %scatter(maxdiffs_activated(idx_suppressed,2), maxdiffs_activated(idx_suppressed,1), 'b.');
    axis([0 5 0 5]);
    line([0 10],[0 10], 'Color', 'k');

    subplot(3,4,2);
    scatter(maxdiffs_activated2(:,2), maxdiffs_activated2(:,1), 'k.');
    hold on;
    %scatter(maxdiffs_activated2(idx_enhanced2,2), maxdiffs_activated2(idx_enhanced2,1), 'r.');
    %scatter(maxdiffs_activated2(idx_suppressed2,2), maxdiffs_activated2(idx_suppressed2,1), 'b.');
    axis([0 5 0 5]);
    line([0 10],[0 10], 'Color', 'k');

    [p,h] = signrank(maxdiffs_activated(:,1),maxdiffs_activated(:,2), 'tail', 'right');
    [p2,h2] = signrank(maxdiffs_activated2(:,1), maxdiffs_activated2(:,2),'tail', 'left');

    print(strcat('maxdiffs', num2str(corr_threshold*100)), '-depsc2');
end
%% SPATIAL DISTRIBUTION
YCOORDINATES = dataarray.(condition).calcium.ycoordinates;
XCOORDINATES = dataarray.(condition).calcium.xcoordinates;
fishIDs = dataarray.(condition).calcium.fishID;
zindices = dataarray.(condition).calcium.zindices;
xmin = 0;
xmax = 240;
ymax = 240;
markersize = 25;

map = [0.4 0 1;  0.5 0.5 1; 102/255 204/255 0; 204/255 0 102/255; 1  102/255 178/255; 1 153/255 51/255; 0 0 0];
map = [0 0 0; 1 140/255 0; 30/255 144/255 1 ;199/255 21/255 133/255; 50/255 205/255 50/255; 0 0 1];

%maxneurons =500; %number of neurons
%randseq = randperm(length(suppressed));

%normalize X and Y coordinates

YXCOORDINATES =[];
fisharray = unique(fishIDs);
YXCOORDINATES(:,1) = YCOORDINATES;
x_width = xmax;
y_width = ymax;

%only normalize x not y
for i = 1:length(fisharray)

    fish = fisharray(i);
    Xpos = XCOORDINATES(fishIDs == fish);
    Ypos = YCOORDINATES(fishIDs == fish);
    norm_Xpos = x_width*(Xpos-min(Xpos))/(max(Xpos)-min(Xpos));
    %norm_Ypos = y_width*(Ypos-min(Ypos))/(max(Ypos)-min(Ypos));
    YXCOORDINATES(fishIDs==fish,2) = norm_Xpos;
    %YXCOORDINATES(fishIDs==fish,1) = norm_Ypos;
end

%%
nplot_perfish = 25; %total number to sample per fish
nplot =nplot_perfish*8; %total number to plot (nplot_perfish*min number of fish)

figure(3);
clf;
figure(11)
clf;
figure(12)
clf;

for s = 2:length(stimchannel) %first channel is water

    for c = 1:length(binary)

        figure(10+c)

        idx_binary = binary{c}(:,s); %active or suppressed (all fish)
        idx_binary_subset = idx_binary(~isnan(idx_binary));
        idx_condition = find(idx_binary_subset ==1);%for fish within this condition
        fishIDs_condition = fishIDs(~isnan(idx_binary)); %fishIDs of fish tested on this condition


        if c == 1
            [idx_subset, new_idx] = resampleindices2(idx_condition, fishIDs_condition, nplot_perfish, nplot);
        else
            Lia = ismember(idx_subset, idx_condition);
            new_idx = idx_subset(Lia.*idx_subset > 0);
        end

        color = map(s, :);
        subplot(2,4,s-1);

        scatter(YXCOORDINATES(idx_subset, 2), YXCOORDINATES(idx_subset, 1),  50, [0.5 0.5 0.5], 'filled');
        hold on;
        scatter(YXCOORDINATES(new_idx, 2), YXCOORDINATES(new_idx,1), 50, color,  'filled');
        set(gca, 'Ydir', 'reverse');
        axis([xmin xmax 0 240]);
        ylabel(char(stimchannel{s}));
        box off;

        box off;
        axis off;

    end

end

figure(11)
print(strcat('spatial-suppressed', '-threshold', num2str(corr_threshold*100)),'-depsc2', '-painters')

figure(12)
print(strcat('spatial-activated', '-threshold', num2str(corr_threshold*100)),'-depsc2', '-painters')


%% Histogram of AP position

figure(15)
clf
figure(16)
clf

ymax = 20;

for c = 1:length(binary)
    for s = 2:length(stimchannel) %first channel is water

        column = binary{c}(:,s).*activated2;
        column(column==0) = nan;
        array_AP = column.*YXCOORDINATES(:,1);
        array_AP = array_AP(~isnan(array_AP));
        figure(14+c);
        subplot(4,2,s-1);
        h = histogram(array_AP(:),'Normalization', 'count');
        h.BinWidth = 8;

        xlim([0 200]);
        if s ==2 && strcmp(condition, 'water')
            ylim([0 50]);
        elseif s~=2 && strcmp(condition, 'water')
            ylim([0 20]);
        elseif strcmp(condition, 'optovin')
            ylim([0 ymax]);
        end

        box off;
        hold on;

        line([nanmedian(array_AP(:)) nanmedian(array_AP(:))], [0 100], 'Color', 'k');

    end
end

figure(15)
print(strcat('APposition-odorsuppressed-TRPA1activated', '-threshold', num2str(corr_threshold*100)),'-depsc2','-painters')
figure(16)
print(strcat('APposition-odoractivated-TRPA1activated', '-threshold', num2str(corr_threshold*100)),'-depsc2', '-painters')

%% BOUT TRIGGERED AVERAGE

%probably also want to filter by bout size
%poststim or nonstim

map_Control = [153/255 204/255 1; 51/255 51/255 1; 0 0 1];
map_Odor = [1 204/255 204/255; 1 102/255 102/255; 1 51/255 51/255];

stimcats = {};
stimcats{1} = 'poststim';
stimcats{2} = 'nonstim';
stimtype = 'teensy';
stimchan = 'kin';

if strcmp(stimtype, 'pump')
    ymin = -0.05;
    ymax = 0.05;
else
    ymin = -0.1;
    ymax = 1;
end

figure(5);
clf;

for s = 1:length(stimcats)

    calcium_raw = [];
    Control_calcium = [];
    Odor_calcium = [];

    stimcat = char(stimcats{s});

    calcium = dataarray.(condition).bout.(stimcat).(stimtype).bouttriggeredcalcium;
    calcium(calcium==0) = nan;
    abs_tail = abs(dataarray.(condition).bout.(stimcat).(stimtype).bouttriggeredtail);
    boutID = dataarray.(condition).bout.(stimcat).(stimtype).boutID;

    if s == 1
        stimID = dataarray.(condition).bout.(stimcat).(stimtype).stimID;
    end

    totaltimepoints = size(calcium,2);
    totaltailpoints = size(abs_tail,2);
    baselineend = floor(tracelength2/2);

    idx_control =  find(strcmp(stimID, 'water'));
    idx_odor = find(strcmp(stimID, stimchan));
    idxtrace_control = [];
    idxtrace_odor = [];


    for b = 1:length(idx_control)
        idxtemp = find(boutID == idx_control(b));
        idxtrace_control(end+1:end+length(idxtemp),1) =  idxtemp;
    end

    for b = 1:length(idx_odor)
        idxtemp = find(boutID == idx_odor(b));
        idxtrace_odor(end+1:end+length(idxtemp),1) =  idxtemp;
    end


    calcium = calcium -nanmean(calcium(:, 1:baselineend), 2);

    subplot(2,2,s*2-1);
    plot(nanmean(calcium(idxtrace_control,:),1), 'Color', map_Control(3,:));
    hold on;
    plot(nanmean(calcium(idxtrace_odor,:),1), 'Color', map_Odor(3,:));

    %errorbars
    plot(nanmean(calcium(idxtrace_control,:),1)+nanstd(calcium(idxtrace_control,:),1)/sqrt(length(idxtrace_control)), 'Color', [0.5 0.5 0.5]);
    plot(nanmean(calcium(idxtrace_control,:),1)-nanstd(calcium(idxtrace_control,:),1)/sqrt(length(idxtrace_control)), 'Color', [0.5 0.5 0.5]);
    plot(nanmean(calcium(idxtrace_odor,:),1)+nanstd(calcium(idxtrace_odor,:),1)/sqrt(length(idxtrace_odor)), 'Color', [0.5 0.5 0.5]);
    plot(nanmean(calcium(idxtrace_odor,:),1)-nanstd(calcium(idxtrace_odor,:),1)/sqrt(length(idxtrace_odor)), 'Color', [0.5 0.5 0.5]);
    %axis([1 size(calcium,2) -0.1 0.1]);

    axis([1 size(calcium,2) ymin ymax]);
    set(gca, 'FontSize', f)
    box off;


    subplot(2,2,s*2);
    plot(nanmean(abs_tail(idx_control,:),1), 'Color', map_Control(3,:));
    hold on;
    plot(nanmean(abs_tail(idx_odor,:),1),  'Color', map_Odor(3,:));
    axis([1 size(abs_tail,2) 0 70]);
    set(gca, 'FontSize', f)
    box off;

end

%% Supplementary Figure 2 - plot effects by isolation time or age
fishID = dataarray.(condition).calcium.fishID;

colorarray = prism(max(fishID));

isolationtimes = [90 60 90 150 150 120 30 180 90 90 60 60 240 90 120 210 150 120 150 210 150 180 90 45 90 120 180 120 180];
ages = [9 9 9 9 9 10 10 11 11 11 8 9 9 9 9 9 10 10 11 11 11 8 8 8 8 9 9 9 9]

calciumsum_byfish = [];

for i = 1:max(fishID)

    %this is if you subtract pre-stimulus response
    calciumsum_byfish(i, 1) = mean(calciumsum_post_byneuron(fishID == i ,1)-calciumsum_pre_byneuron(fishID == i ,1));
    calciumsum_byfish(i, 2) = mean(calciumsum_post_byneuron(fishID == i ,2)-calciumsum_pre_byneuron(fishID == i ,2));

    %this is if you don't - currently I am not since it doesn't make much
    %of a diff
    calciumsum_byfish(i, 1) = mean(calciumsum_post_byneuron(fishID == i ,1));
    calciumsum_byfish(i, 2) = mean(calciumsum_post_byneuron(fishID == i ,2));

end

%BY ISOLATION TIME
unit = 4;
figure(10)
clf
subplot(2,1,1);
scatter(isolationtimes, calciumsum_byfish(:, 1), 20, 'k', 'filled');
hold on;
scatter(isolationtimes+unit , calciumsum_byfish(:, 2), 20, 	[255/255 165/255 0], 'filled');

for i = 1:max(fishID)
    plot([isolationtimes(i) isolationtimes(i)+unit], calciumsum_byfish(i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end

xlabel('Isolation time (min)')

subplot(2,1,2);
scatter(isolationtimes, calciumsum_byfish(:, 1)-calciumsum_byfish(:, 2), 20, 'k', 'filled');

xlabel('Isolation time (min)')

%BY AGE
unit2 = 0.2
figure(11)
clf
subplot(2,1,1);
scatter(ages, calciumsum_byfish(:, 1), 20, 'k', 'filled');
hold on;
scatter(ages+unit2 , calciumsum_byfish(:, 2), 20, 	[255/255 165/255 0], 'filled');

for i = 1:max(fishID)
    plot([ages(i) ages(i)+unit2], calciumsum_byfish(i,:), 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
end

xlabel('Age (dpf)')
xlim([7.5 11.5])

subplot(2,1,2);
scatter(ages, calciumsum_byfish(:, 1)-calciumsum_byfish(:, 2), 20, 'k', 'filled');
xlim([7.5 11.5])


%%
save(sprintf('%s-%s%s', 'analyzeddata', celltype, '.mat'));