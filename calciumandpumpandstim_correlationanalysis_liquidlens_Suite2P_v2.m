%% To analyze calcium imaging data from Suite2P-extracted units
%% SECTION 1 is for aggregating raw data into structured array
% SKIP TO SECTION 2 if dataarray is already generated
% Email me for dataarray to test code (too large to upload)

clf; close all; clear all;

celltype = sprintf('OB');
%experiment = 'kinnonkin';
experiment = 'kinnonkinadult';
condition = 'water';
nchannel =1;

minframes = 2800; %6 rounds

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

disk = 21;

pathname =  sprintf('%s%d%s/%s/%s/', '/Volumes/CW-', disk, '/2P_olfactorybulb_new/aggregated/', condition, experiment)
cd('/Users/carolinewee/Dropbox/Carolinebackup/OXTai/epsfiles/olfactory/suite2P/'); %save folder
%pathname =  sprintf('%s%d%s/%s/%s/', '/Volumes/CW-', disk, '/2P_socialodor/GABAoxy/aggregated/', condition, experiment)
%cd(sprintf('%s/%s/', '/Users/carolinewee/Dropbox/Carolinebackup/OXTai/epsfiles/social/GABAoxy/',experiment)); %savefolder

%pathname =  sprintf('%s%d%s/%s/%s/', '/Volumes/CW-', disk, '/2P_socialodor/GABAonly/aggregated/', condition, experiment)
%cd(sprintf('%s/%s/', '/Users/carolinewee/CW Lab Dropbox/Caroline Wee/Carolinebackup/OXTai/epsfiles/social/GABAoxy/',experiment)); %savefolder
%GABA only still saved in GABAoxy

files = dir(pathname);

FISH = 0;

conditionarray = {};
conditionarray{1} = 'water'
conditionarray{2} = 'optovin'

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
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'redcells', []);
    dataarray.(condition_tmp).calcium = setfield(dataarray.(condition_tmp).calcium, 'oxycells', []);
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
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_stimreg', struct);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_stimreg2', struct);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_stimreg_pvalues', struct);
    dataarray.(condition).calcium = setfield(dataarray.(condition).calcium, 'xcorr_stimreg_pvalues2', struct);

end

for i =1:length(files)

    name = files(i).name;
    findstring = strfind(name, '.mat');
    if isempty(findstring)
        continue;
    else
        filename = sprintf('%s/%s', pathname, name);
        load(filename)
        FISH = FISH+1;
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

            %%

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

            %%
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
            
            try %this is to exclude bouts before the experiment begins
                idx_expstart =  find(calciumtaildata.(condition).bout.poststimbouts.teensy>0);
                expstart = idx_expstart(1);
                expend = idx_expstart(end);
            catch

                expstart = 1;
                expend = length(nonstimbouts);

            end
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
            dataarray.(condition).stim.(stimtype).calciumpeaks(FISH, 1:nstims) = mean(calciumtaildata.(condition).stim.(stimtype).peakcalcium(1:nstims, :),2);
            dataarray.(condition).stim.(stimtype).calciumsum(FISH, 1:nstims) = mean(calciumtaildata.(condition).stim.(stimtype).sumcalcium(1:nstims, :), 2);


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

    % to aggregate dfoverf and others

    %to aggregate fishID of each neuron
    nneurons = length(calciumtaildata.(condition).spatial.npix);
    dataarray.(condition).calcium.fishID(end+1:end+nneurons, 1) = FISH*ones(nneurons,1);

    % to aggregate ROI sizes
    dataarray.(condition).calcium.ROIsizes(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.npix(1:nneurons);

    %to aggregate indices of z-plane numbers
    dataarray.(condition).calcium.zindices(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.zpos(1:nneurons);

    %to aggregate coordinates - use the coordinates that are normalized and
    %and centered
    dataarray.(condition).calcium.ycoordinates(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.ypos_norm2;
    dataarray.(condition).calcium.xcoordinates(end+1:end+nneurons, 1) = calciumtaildata.(condition).spatial.xpos_norm2;

    % to aggregate redcells
    try
        dataarray.(condition).calcium.redcells(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.redcell(1:nneurons);
    catch
    end

    % to aggregate oxycells
    try
        dataarray.(condition).calcium.oxycells(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.oxycell(1:nneurons);
    catch
    end

    %XCORR
    %try %if XCORR values have been calculated;
    dataarray.(condition).calcium.xcorr_tailreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_tailreg;
    dataarray.(condition).calcium.xcorr_motorstimreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_motorstimreg;
    dataarray.(condition).calcium.xcorr_motorsponreg(end+1:end+nneurons, 1) = calciumtaildata.(condition).calcium.xcorr_motorsponreg;

    % aggregate regressors according to name
    stimreg_name = calciumtaildata.(condition).calcium.xcorr_stimregname;

    for c = 1:length(stimreg_name)

        try
            dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))(end+1:end+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg(:,c);
        catch

            dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c})) = [];
            dataarray.(condition).calcium.xcorr_stimreg.(char(stimreg_name{c}))(end+1:end+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg(:,c);
            %dataarray.(condition).calcium.xcorr_stimreg2.(char(stimreg_name{c})) = [];
            %dataarray.(condition).calcium.xcorr_stimreg2.(char(stimreg_name{c}))(end+1:end+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg2(:,c);
            dataarray.(condition).calcium.xcorr_stimreg_pvalues.(char(stimreg_name{c})) = [];
            dataarray.(condition).calcium.xcorr_stimreg_pvalues.(char(stimreg_name{c}))(end+1:end+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg_pvalues(:,c);
            %dataarray.(condition).calcium.xcorr_stimreg_pvalues2.(char(stimreg_name{c})) = [];
            %dataarray.(condition).calcium.xcorr_stimreg_pvalues2.(char(stimreg_name{c}))(end+1:end+nneurons, :) = calciumtaildata.(condition).calcium.xcorr_stimreg_pvalues2(:,c);

        end
    end

    %to aggregate dfoverf
    dfoverf_consolidated = calciumtaildata.(condition).calcium.dfoverf;
    dfoverf_smoothed = calciumtaildata.(condition).calcium.dfoverf_smoothed;

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

end

%clear calciumtaildata
%clearvars -except dataarray celltype exps folderpath_aggregated exptype disk nchannel condition

f = 15;
lw = 2;
zstart = 1;
zend = 3;

%% SECTION 2: THIS IS WHERE ANALYSIS STARTS
% Requires dataarray

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

    %fishmeans_control(isnan(fishmeans_control)) = 0;
    %fishmeans_treatment(isnan(fishmeans_treatment)) = 0;

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

% all neurons
close all;
stimtype = 'pump'
neurontype = 'all'; %'all', 'red', 'oxy', 'subpallial', 'vTel', 'PO'

if strcmp(condition, 'optovin')
    ymin = -0.8;
    ymax =1.8;
else
    ymin = -0.3;
    ymax = 0.8;  %0.4 0.8
end

dataarray.(condition).stim.(stimtype).tailkinematics = struct;

stimchannel = {};
stimchannel{1} = 'control';
stimchannel{2} = 'odor1';
stimchannel{3} = 'odor2';
stimchannel{4} = 'odor3';
stimchannel{5} = 'odor4';

calciumsum_pre = [];
calciumsum_post = [];

calciummax_pre = [];
calciummax_post = [];

calciumsum_prepump = [];
calciumsum_postpump = [];

trace_byneuron = {};

calciumsum_post_byneuron = [];
calciumsum_pre_byneuron = [];
calciummax_post_byneuron = [];
calciummax_pre_byneuron = [];
calciumsum_prepump_byneuron = [];
calciumsum_postpump_byneuron = [];

if strcmp(condition, 'water')
    window = 60; %this is the range in which you average calcium responses
elseif strcmp(condition, 'optovin')
    window = 5;
end

if strcmp(stimtype, 'teensy')
    window = 5;
end

map = [0 0 0; 1 140/255 0; 30/255 144/255 1 ;199/255 21/255 133/255; 50/255 205/255 50/255; 0 0 1]

statstable = {};
abs_tailtraces = struct;
abs_tailtraces.tails = struct;
abs_tailtraces.fish = struct;


neurontypearray{1} = 'GABA';
neurontypearray{2} = 'subpallial GABA'; %aSpa
neurontypearray{3} = 'vTel GABA'; %pSpa
neurontypearray{4} = 'PO GABA';
neurontypearray{5} = 'PO nonGABA';
neurontypearray{6} = 'oxy';
neurontypearray{7} = 'all';

nstimchannels = length(fieldnames(dataarray.(condition).calcium.xcorr_stimreg));

% FIRST to scale data across all fish for determining neuron type

xpos_scaled = [];
ypos_scaled = [];

xpos = dataarray.(condition).calcium.xcoordinates;
ypos = dataarray.(condition).calcium.ycoordinates;
zpos = dataarray.(condition).calcium.zindices;
redcells = dataarray.(condition).calcium.redcells;
oxycells = dataarray.(condition).calcium.oxycells;
fishID = dataarray.(condition).calcium.fishID;

if strcmp(celltype,'OB')
    xrange = 320;
    yrange = xrange/2;
elseif strcmp(celltype,'GABAoxy')
    yrange = 320;
    xrange = yrange/2;
elseif strcmp(celltype, 'GABAonly')
    yrange = 320;
    xrange = yrange/3*2;
end

for fish = 1:max(fishID)
    
    xpos_perfish = xpos(fishID == fish);
    ypos_perfish = ypos(fishID == fish);
    
    max_x = max(xpos_perfish);
    max_y = max(ypos_perfish);
    
    new_xposperfish = (xpos_perfish/max_x)*xrange;
    new_yposperfish = (ypos_perfish/max_y)*yrange;
    
    xpos_scaled(end+1: end+length(xpos_perfish), 1) = new_xposperfish;
    ypos_scaled(end+1: end+length(ypos_perfish), 1) = new_yposperfish;
    
end

for type = 1:length(neurontypearray)
    
    
    neurontype = neurontypearray{type};
    
    pause(0.5);
    figure(4);
    clf;
    
    
    for s = 1:nstimchannels

        stimchan = char(stimchannel{s});

        calcium_raw = dataarray.(condition).stim.(stimtype).stimtriggeredcalcium.(stimchan);
        abs_tail = abs(dataarray.(condition).stim.(stimtype).stimtriggeredtail.(stimchan));
        tail = dataarray.(condition).stim.(stimtype).stimtriggeredtail.(stimchan);

        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan) = struct;

        % same for control and odor (only in these experiments)
        frameunit = mean(diff(calciumtaildata.(condition).tail.tailtime))*1000; %why are we using calciumtaildata?
        peakangles = max(abs(tail), [], 2);
        sumangles = trapz(abs(tail), 2)*frameunit;
        maxvelocities = max(abs(diff(tail,1,2)), [], 2)/frameunit;
        meanvelocities = nanmean(abs(diff(tail,1,2)), 2)/frameunit;

        peakangles(peakangles ==0) = nan;
        sumangles(sumangles==0) = nan;
        maxvelocities(maxvelocities==0) = nan;
        meanvelocities(meanvelocities==0) = nan;

        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).peakangle = peakangles;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).sumangle = sumangles;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).maxvelocity = maxvelocities;
        dataarray.(condition).stim.(stimtype).tailkinematics.(stimchan).meanvelocity = meanvelocities;

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

        % same for control or stim
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

            calciumsum_post(n,s) = nansum(calcium(n, baselineend+stimstart:baselineend+stimstart+window))/window; %average response/s
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
        redthreshold = 0.4;  %to determine if cell is red cell

        if strcmp(neurontype, 'GABA')
            redcells = dataarray.(condition).calcium.redcells;
            neuronvector = find(redcells > redthreshold);

        elseif strcmp(neurontype, 'subpallial GABA') && contains(pathname, 'GABAoxy') %set y value thresholds
            neuronvector = find(redcells > redthreshold & ypos_scaled<=100);

        elseif strcmp(neurontype, 'subpallial nonGABA') && contains(pathname, 'GABAoxy') %set y value thresholds
            neuronvector = find(redcells < redthreshold & ypos_scaled<=100);

        elseif strcmp(neurontype, 'vTel GABA') && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells > redthreshold & ypos_scaled<=180 & ypos_scaled>=100)

        elseif strcmp(neurontype, 'vTel nonGABA') && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells < redthreshold & ypos_scaled<=180 & ypos_scaled>=100)

        elseif strcmp(neurontype, 'PO nonGABA')  && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells < redthreshold & ypos_scaled> 180 & oxycells ~=1);

        elseif strcmp(neurontype, 'PO GABA')  && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells > redthreshold & ypos_scaled> 180 & oxycells ~=1);

        elseif strcmp(neurontype, 'oxy')
            oxycells = dataarray.(condition).calcium.oxycells;
            neuronvector = find(oxycells ==1);

        elseif strcmp(neurontype, 'subpallial GABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells > redthreshold & ypos_scaled<=200);

        elseif strcmp(neurontype, 'subpallial nonGABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells < redthreshold & ypos_scaled<=200);

        elseif strcmp(neurontype, 'vTel GABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells > redthreshold & ypos_scaled>=200);

        elseif strcmp(neurontype, 'vTel nonGABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells < redthreshold & ypos_scaled>=200);

        else
            ROIsizes =  dataarray.(condition).calcium.ROIsizes; %no particular reason why ROI sizes, just to get right length
            neuronvector = 1:length(ROIsizes);
        end


        for n = 1: length(neuronvector) %scroll through all desired neurons

            idx_n = find(neuronID == neuronvector(n));

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
        %axis([1 size(calcium,2) ymin ymax]);
        axis([1 177 ymin ymax]);

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
        axis([1 size(abs_tail,2) 0 100]);
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

        statstable{s} = calcium;

    end

    print(strcat('stimtriggeredaverage', condition, stimtype, '_', neurontype), '-depsc2')

end


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

%% compare all signals to:
% 1) baseline (whether activated or suppressed by the odor)
% 2) water (whether activated relative to water
%need to make sure neurontype was "all" in STIMTRIGGEREDAVERAGE section

crosscomparisons = struct; %save it somewhere later?
neurontype = 'all'; %'all', 'red', 'oxy', 'subpallial', 'vTel', 'PO'

stimcatarray = {};
stimcatarray{1} = 'water'
stimcatarray{2} = firstodor;
stimcatarray{3} = secondodor;
stimcatarray{4} = thirdodor;

stimchannel = {};

if strcmp(celltype, 'OB') || strcmp(celltype, 'GABAonly')
    stimchannel{1} = 'control';
    stimchannel{2} = 'odor1';
    stimchannel{3} = 'odor2';
    stimchannel{4} = 'odor3';

else
    stimchannel{1} = 'control';
    stimchannel{2} = 'odor1';
    stimchannel{3} = 'odor2';
end

%NOTE: need to run stim triggered average again with stim type 'pump'!!!!
xcorr_stimreg = dataarray.(condition).calcium.xcorr_stimreg;
%xcorr_stimreg_pvalues = dataarray.(condition).calcium.xcorr_stimreg_pvalues;

%to replace non-signficant values with zeros
% pcutoff = 0.05;         %pvalue cutoff
% for s = 1:nstimchannels
%     cat = stimcatarray{s};
%     xcorr_stimreg.(cat)(xcorr_stimreg_pvalues.(cat)>=pcutoff) = 0;
% end

xcorr_motorstimreg = dataarray.(condition).calcium.xcorr_motorstimreg;
ypos = dataarray.(condition).calcium.ycoordinates;
useregressor = 1; % if zero, use calcium sum

%need to be same length
corr_thresholds = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.1];
diff_thresholds = [0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.05];
%diff_thresholds = [0 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.05];

fractions = {};
% now average both calcium activity and average activity by neuron
redthreshold = 0.4;  %to determine if cell is red cell

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
            %
            if useregressor == 0

                try
                    if calciumsum_post_byneuron(idx_n,s) <= calciumsum_pre_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,s)-calciumsum_pre_byneuron(idx_n,s))>= diffthreshold
                        activity{1}(n,s) = 1; %suppressed2
                    elseif calciumsum_post_byneuron(idx_n,s) >= calciumsum_pre_byneuron(idx_n,s) && abs(calciumsum_post_byneuron(idx_n,s)-calciumsum_pre_byneuron(idx_n,s))>= diffthreshold
                        activity{2}(n,s) = 1; %activated2
                    else
                        activity{1}(n,s) = 0;
                        activity{2}(n,s) = 0;
                    end

                catch %if column is shorter than nneurons in control (water condition)
                    activity{1}(n,s) = nan;
                    activity{2}(n,s) = nan;
                end

            elseif useregressor == 1 %use correlation coefficient
                %
                try
                    % activity per odor
                    if xcorr_stimreg.(cat)(idx_n,1) <= -corrthreshold %&& xcorr_stimreg_pvalues.(cat)(idx_n,1)< pcutoff
                        activity{1}(n,s) = 1; %suppressed2
                        activity{2}(n,s) = 0;
                    elseif xcorr_stimreg.(cat)(idx_n,1) >= corrthreshold %&& xcorr_stimreg_pvalues.(cat)(idx_n,1)< pcutoff
                        activity{2}(n,s) = 1; %activated2
                        activity{1}(n,s) = 0;
                    else
                        activity{1}(n,s) = 0;
                        activity{2}(n,s) = 0;
                    end

                catch
                    activity{1}(n,s) = nan;
                    activity{2}(n,s) = nan;
                end

            end

            % odor relative to water - use difference

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

            %%% odor relative to water
            %
            %                 if xcorr_stimreg(n,1) > xcorr_stimreg(n,s) && abs(xcorr_stimreg(n,1)-xcorr_stimreg(n,s))>= corrthreshold
            %                     binary{1}(n,s) = 1; %suppressed relative to water
            %                 elseif xcorr_stimreg(n,1) < xcorr_stimreg(n,s) && abs(xcorr_stimreg(n,1)-xcorr_stimreg(n,s))>=corrthreshold
            %                     binary{2}(n,s) = 1; %activated relative to water
            %                 else
            %                     binary{1}(n,s) = 0;
            %                     binary{2}(n,s) = 0;
            %                 end
            %
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

% quantify overlap between odor responses

%first combine suppressed and enhanced neurons
activity_combined = [];
activity_suppressed = activity{1,1};
activity_suppressed(activity_suppressed == 1) = -1;
activity_activated =  activity{1,2};
activity_combined = activity_suppressed + activity_activated;

activity_odor = [];
xcorr_odor = [];

if strcmp(celltype, 'OXTPO')
    orderarray = [nan 2 1; 2 3 1; 2 4 1] %non-kin, then adult
elseif nstimchannels == 3
    orderarray = [2 3 1] %change if more odors
elseif nstimchannels == 4
    orderarray = [2 3 4 1]
elseif nstimchannels ==2
    orderarray = [2 1]
end

index = 1;
sortedactivityarray = {};
catnamesarray = {};

for m = 1:size(orderarray,1)

    catnames = {};
    order = orderarray(m, :);
    sortedactivityarray{m} = [];
    catnamesarray{m} = {};
    activity_odor = [];
    xcorr_odor = [];


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

        column = xcorr_stimreg.(cat);

        activity_odor(1:length(neuronvector),index) = activity_combined(neuronvector,s); %combine activity for each odor including water

        % need to also reduce xcorr to subset of cells
        xcorr_odor(1:length(neuronvector),index) = xcorr_stimreg.(cat)(neuronvector,1);

        catnames{1,index} = cat;
        index = index + 1;

    end

    startnum = 1;
    idxold = [];
    sortedactivity = [];


    for s = 1:length(order)-1
        % this is the column by which to sort by
        column = activity_odor(:,s);

        idxk = find(column==1);
        idxk2 = find(column==-1);
        idxk3 = find(column==0);

        %do not want to sort neurons that have already been sorted
        idx_unique = setdiff(idxk, idxold);
        idx_unique2 = setdiff(idxk2, idxold);
        idx_unique3 = setdiff(idxk3, idxold);

        %sort of xcorr values
        %[sorted, idx_sorted] = sort(xcorr_odor(idx_unique,s), 'descend')

        %but you are sorting ALL the columns ultimately
        sortedactivity(startnum:startnum+length(idx_unique)-1, 1:size(activity_odor,2)) = xcorr_odor(idx_unique, :);
        sortedactivity(startnum+length(idx_unique):startnum+length(idx_unique)+length(idx_unique2)-1,1:size(activity_odor,2)) = xcorr_odor(idx_unique2,:);
        %sortedactivity(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = xcorr_odor(idx_unique3,:);

        if s == length(order)-1
            sortedxcorr(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(activity_odor,2)) = xcorr_odor(idx_unique3,:);
        end


        %idxold = vertcat(idxold, idx_unique, idx_unique2, idx_unique3);
        idxold = vertcat(idxold, idx_unique, idx_unique2);
        startnum = length(idxold) + 1;
    end

    catnamesarray{m} = catnames;
    sortedactivityarray{m} = sortedactivity;

end

%% plot correlation maps

%to generate a color map red to blue
if strcmp(celltype,'OB')
    corrlimits = [-0.6 0.6];
else
    corrlimits = [-0.3 0.3];
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
clf;
for i = 1:length(sortedactivityarray)
    subplot(1,3,i)
    colormap(cmap_fine)
    imagesc(sortedactivityarray{i}, corrlimits)
    xlabel(catnames);
    ylabel(num2str(length(sortedactivityarray{i})))
    box off;
    print(strcat('clusteractivity_', neurontype), '-deps')
end

%% quantify overlap
%cross compare active neurons

crosscomparisons.(neurontype).overlap = struct;
crosscomparisons.(neurontype).exclusive = struct;
crosscomparisons.(neurontype).nneurons = size(activity_odor,1);

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

xcorr_stimreg = dataarray.(condition).calcium.xcorr_stimreg;
catnames = stimcatarray(1:length(order));

clusterarray = [];
for i = 1:length(catnames)
    clusterarray(:,i) = xcorr_stimreg.(catnames{i});
end

%% kmeans
k = 12;
idx = kmeans(clusterarray, k, 'Display', 'iter', 'OnlinePhase', 'on', 'Replicates', 5);
%% save clusters at specific time point
c = clock;
time = strcat(num2str(c(4)), num2str(c(5)));
save(sprintf('%s-%s-%s-%s%s', 'analyzeddata', celltype, date, time, '.mat'));
%%
%to generate a color map red to blue

if strcmp(celltype, 'OB')
    corrlimits = [-0.6 0.6];
else
    corrlimits = [-0.6 0.6];
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

sortedcluster = [];
sorteddfoverf = [];
clustersizes = [];
clustermeans = [];
sorted_neuronID = [] ;
sorted_idx = [] ;
sorted_fishID = [];

dfoverf = dfoverf_concat;
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
end

fishID = dataarray.(condition).calcium.fishID;

for s = rearrange
    % this is the column by which to sort by

    %but you are sorting ALL the columns ultimately
    sortedcluster(length(sortedcluster)+1: length(sortedcluster) + length(find((idx==s))), 1:size(clusterarray,2)) = clusterarray(idx==s, :);
    sorteddfoverf(size(sorteddfoverf,1)+1: size(sorteddfoverf,1) + length(find((idx==s))), 1:size(dfoverf,2)) = dfoverf(idx==s, :);
    clustersizes(s,1) =  length(find((idx==s)));
    sorted_neuronID(length(sorted_neuronID)+1: length(sorted_neuronID) + length(find((idx==s))), 1) = find(idx==s);
    sorted_idx(length(sorted_idx)+1: length(sorted_idx) + clustersizes(s), 1) = repmat(s, [clustersizes(s), 1]);

    clustermeans(s, 1:size(dfoverf,2)) = mean(dfoverf(idx==s, :),1);
    sorted_fishID(length(sorted_fishID)+1: length(sorted_fishID) + length(find((idx==s))), 1) = fishID(find(idx==s));
end

figure(16)
clf;
subplot(1,2,1);
colormap(cmap_fine)
imagesc(sortedcluster, corrlimits)
xlabel(catnames);
ylabel(num2str(maxneuron))
box off;
colorbar;

colorarray = hsv(k);
count = 1;
for s = rearrange
    color = colorarray(s, :);
    text(size(sortedcluster,2)+3, count, strcat(num2str(s), ': n = ', num2str(clustersizes(s)), ' (', num2str(clustersizes(s)/nneurons*100), '%)'), 'Color', color)
    count = count+ clustersizes(s);
    hold on;
    ylim([0 nneurons]);
    set(gca, 'YDir', 'reverse')
end


print(strcat('kmeansclustersorted_', neurontype), '-depsc2', '-painters')

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

%% now plot by cell type

if strcmp(celltype, 'GABAoxy') || strcmp(celltype, 'GABAonly')

    xpos = dataarray.(condition).calcium.xcoordinates;
    ypos = dataarray.(condition).calcium.ycoordinates;
    zpos = dataarray.(condition).calcium.zindices;
    redcells = dataarray.(condition).calcium.redcells;
    oxycells = dataarray.(condition).calcium.oxycells;

    % to normalize
    xpos_scaled = [];
    ypos_scaled = [];

    redthreshold = 0.4;

    if strcmp(celltype,'OB')
        xrange = 320;
        yrange = xrange/2;
    elseif strcmp(celltype,'GABAoxy')
        yrange = 320;
        xrange = yrange/2;
    elseif strcmp(celltype, 'GABAonly')
        yrange = 320;
        xrange = yrange/3*2;
    end

    allclustersubsets = {};%for aggregating data
    allclustersizes = {}; %for aggregating data

    for fish = 1:max(fishID)

        xpos_perfish = xpos(fishID == fish);
        ypos_perfish = ypos(fishID == fish);

        max_x = max(xpos_perfish);
        max_y = max(ypos_perfish);

        new_xposperfish = (xpos_perfish/max_x)*xrange;
        new_yposperfish = (ypos_perfish/max_y)*yrange;

        xpos_scaled(end+1: end+length(xpos_perfish), 1) = new_xposperfish;
        ypos_scaled(end+1: end+length(ypos_perfish), 1) = new_yposperfish;

    end

    figure(19)
    clf
    figure(20)
    clf
    neurontypearray = {};

    neurontypearray{1} = 'GABA';
    neurontypearray{2} = 'subpallial GABA';
    neurontypearray{3} = 'vTel GABA';

    if strcmp(celltype, 'GABAoxy')
        neurontypearray{4} = 'PO GABA';
        neurontypearray{5} = 'PO nonGABA';
        neurontypearray{6} = 'oxy';
    end

    clustersizes_percelltype = zeros(k, length(neurontypearray));
    clustersizes_perfish_percelltype = {};
    for fish = 1:max(fishID)
        clustersizes_perfish_percelltype{fish} = zeros(k, length(neurontypearray));
    end


    for type = 1:length(neurontypearray)

        neurontype = neurontypearray{type}

        if strcmp(neurontype, 'GABA')
            redcells = dataarray.(condition).calcium.redcells;
            neuronvector = find(redcells > redthreshold);

        elseif strcmp(neurontype, 'subpallial GABA') && contains(pathname, 'GABAoxy') %set y value thresholds
            neuronvector = find(redcells > redthreshold & ypos_scaled<=100);

        elseif strcmp(neurontype, 'vTel GABA') && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells > redthreshold & ypos_scaled<=180 & ypos_scaled>=100)

        elseif strcmp(neurontype, 'PO nonGABA')  && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells < redthreshold & ypos_scaled> 180 & oxycells ~=1);

        elseif strcmp(neurontype, 'PO GABA')  && contains(pathname, 'GABAoxy')
            neuronvector = find(redcells > redthreshold & ypos_scaled> 180 & oxycells ~=1);

        elseif strcmp(neurontype, 'oxy')
            oxycells = dataarray.(condition).calcium.oxycells;
            neuronvector = find(oxycells ==1);

            %GABAonly - more dorsal, to exclude OB neurons
        elseif strcmp(neurontype, 'subpallial GABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells > redthreshold & ypos_scaled<=200 & ypos_scaled>100 & zpos>1 );

            %GABAonly - more ventral,
        elseif strcmp(neurontype, 'subpallial GABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells > redthreshold & ypos_scaled<=200 & zpos == 1);


        elseif strcmp(neurontype, 'vTel GABA') && contains(pathname, 'GABAonly')
            neuronvector = find(redcells > redthreshold & ypos_scaled>=200);

        else
            ROIsizes =  dataarray.(condition).calcium.ROIsizes; %no particular reason why ROI sizes, just to get right length
            neuronvector = 1:length(ROIsizes);
        end

        cluster_subset = [];
        dfoverf_subset = [];

        for n = 1:nneurons

            if intersect(sorted_neuronID(n), neuronvector)

                cluster_subset(end+1,:) = sortedcluster(n,:);
                dfoverf_subset(end+1,:) = sorteddfoverf(n,:);

                clustersizes_percelltype(sorted_idx(n), type) = clustersizes_percelltype(sorted_idx(n,:), type)+1;
                clustersizes_perfish_percelltype{sorted_fishID(n)}(sorted_idx(n), type) = clustersizes_percelltype(sorted_idx(n,:), type)+1;
            else
            end
        end

        allclustersubsets{type} = cluster_subset;


        figure(19)
        subplot(3,3,type)

        colormap(cmap_fine)
        imagesc(cluster_subset, corrlimits)
        %xlabel(catnames)
        ylabel(num2str(length(neuronvector)))
        title(neurontype)
        box off;

        print('kmeansclustersorted_alltypes', '-depsc2')

        figure(20)
        colormap(parula)
        subplot(3,3, type)
        clims = [-0.5 2];
        imagesc(dfoverf_subset, clims)
        ylabel(num2str(length(neuronvector)))
        title(neurontype)
        box off;

        print('kmeansclusterraster_alltypes', '-depsc2')


    end


    dataarray.(condition).calcium.kmeansxcorr_clustersizespercelltype = clustersizes_percelltype;


else
end

% plot cluster sizes per cell type

figure(50)
clf

  for type = 1:length(neurontypearray)
     tmparray = clustersizes_percelltype(:,type);

     %need to rearrange
     tmparray2 = [];

     for s = rearrange
         tmparray2(end+1) = tmparray(s);
     end

     allclustersizes{type} = tmparray2';

     figure(50)
     subplot(4,3,type)
     bar(tmparray2)
     xlim([0 13])

     box off

  end

  %% spatial distribution for each cluster

  markersize = 25;
  redthreshold = 0.4; %threshold for determining if cell is red. default is 0.65.

  condition = 'water';

  xpos = dataarray.(condition).calcium.xcoordinates;
  ypos = dataarray.(condition).calcium.ycoordinates;
  zpos = dataarray.(condition).calcium.zindices;
  redcells = dataarray.(condition).calcium.redcells;
  oxycells = dataarray.(condition).calcium.oxycells;
  fishID = dataarray.(condition).calcium.fishID;


  if strcmp(celltype,'OB')
      xrange = 320;
      yrange = xrange/2;
      viewxy =  [35 20]
  elseif strcmp(celltype,'GABAoxy')
      yrange = 320;
      xrange = yrange/2;
      viewxy =  [27 13]
  elseif strcmp(celltype, 'GABAonly')
      yrange = 320;
      xrange = yrange/3*2;
      viewxy =  [41 15]
  end

  % to normalize
  xpos_scaled = [];
  ypos_scaled = [];
  close all;

  for fish = 1:max(fishID)

      xpos_perfish = xpos(fishID == fish);
      ypos_perfish = ypos(fishID == fish);

      max_x = max(xpos_perfish);
      max_y = max(ypos_perfish);

      new_xposperfish = (xpos_perfish/max_x)*xrange;
      new_yposperfish = (ypos_perfish/max_y)*yrange;

      xpos_scaled(end+1: end+length(xpos_perfish), 1) = new_xposperfish;
      ypos_scaled(end+1: end+length(ypos_perfish), 1) = new_yposperfish;

  end

  % 3D Plot for each cluster

  markersize = 15;

  colorarray = hsv(max(idx));

  for cluster = 1:max(idx) %we still keep to the original sorted clusters but relabel in illustrator

      color = colorarray(cluster, :);

      neuronids = find(idx==cluster);
      xpos_cluster = xpos_scaled(neuronids);
      ypos_cluster = ypos_scaled(neuronids);
      zpos_cluster = zpos(neuronids);

      neuronids_other = find(idx~= cluster);
      xpos_notcluster = xpos_scaled(neuronids_other);
      ypos_notcluster = ypos_scaled(neuronids_other);
      zpos_notcluster = zpos(neuronids_other);

      figure(35)
      set(gca, 'YDir', 'reverse')
      subplot(1,2,1)
      scatter3(xpos_cluster, ypos_cluster, zpos_cluster,  markersize , color, 'filled');
      hold on;
      %scatter3(xpos_cluster(redcells>=redthreshold), ypos_cluster(redcells>=redthreshold), zpos_cluster(redcells>=redthreshold),  markersize , 'k');

      title('all clusters')
      ylim([0 320]);
      xlim([0 320]);
      view(viewxy)
      %axis off

      %individual clusters
      figure(22 + cluster)
      set(gca, 'YDir', 'reverse')
      subplot(1,2,1)
      scatter3(xpos_notcluster, ypos_notcluster, zpos_notcluster,  markersize , [0.8 0.8 0.8], 'filled');
      hold on;
      scatter3(xpos_cluster, ypos_cluster, zpos_cluster,  markersize , color, 'filled');
      set(gca, 'YDir', 'reverse')
      view(viewxy)
      %axis off

      title(strcat('cluster', num2str(cluster)))
      ylim([0 320]);
      xlim([0 320]);
      figure(22 + cluster)
      print(strcat('spatial3D-kmeans-xcorr-cluster', num2str(cluster)), '-depsc2', '-painters')

  end

  figure(35)
  print('spatial3D-kmeans-xcorr', '-depsc2', '-painters')

  %% GMM
  k = 1:20;
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

  %% quantify overlap between odor responses

  %if responds to kin, does it also respond to other odors? If so which
  %direction?
  %use calcium sum

  %first combine suppressed and enhanced neurons
  binary_combined = [];
  binary_suppressed = binary{1,1};
  binary_suppressed(binary_suppressed == 1) = -1;
  binary_activated =  binary{1,2};
  binary_combined = binary_suppressed + binary_activated;

  binary_odor = [];
  calciumdiff_odor = [];

  shortestcolumn = activity_combined(:,3);%not all fish were imaged with nonkin odor
  maxneuron = min(find(isnan(shortestcolumn)))-1;

  if isempty(maxneuron)
      maxneuron = size(binary_combined,1);
  end

  for s = 2:size(binary_combined,2) %number of categories

      binary_odor(1:maxneuron,s-1) = binary_combined(1:maxneuron,s); %combine binaries for each odor
      calciumdiff_odor(1:maxneuron, s-1) = calciumsum_post_byneuron(neuronvector(1:maxneuron),s) - calciumsum_post_byneuron(neuronvector(1:maxneuron),1); %difference from water
  end

  startnum = 1;
  idxold = [];
  sortedbinary = [];

  for s = 1:size(binary_combined,2)-1

      % this is the column by which to sort by
      column = binary_odor(:,s);

      idxk = find(column== 1);
      idxk2 = find(column==0);
      idxk3 = find(column==-1);

      %do not want to sort neurons that have already been sorted
      idx_unique = setdiff(idxk, idxold);
      idx_unique2 = setdiff(idxk2, idxold);
      idx_unique3 = setdiff(idxk3, idxold);

      %but you are sorting ALL the columns ultimately
      sortedbinary(startnum:startnum+length(idx_unique)-1, 1:size(binary_odor,2)) = calciumdiff_odor(idx_unique, :);
      sortedbinary(startnum+length(idx_unique):startnum+length(idx_unique)+length(idx_unique2)-1,1:size(binary_odor,2)) = calciumdiff_odor(idx_unique2,:);
      sortedbinary(startnum+length(idx_unique)+length(idx_unique2):startnum+length(idx_unique)+length(idx_unique2)+length(idx_unique3)-1, 1:size(binary_odor,2)) = calciumdiff_odor(idx_unique3,:);

      idxold = vertcat(idxold, idx_unique, idx_unique2, idx_unique3);
      startnum = length(idxold) + 1;

  end

  % plot correlation maps
  %to generate a color map red to blue
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
  clf;
  colormap(cmap_fine)
  imagesc(sortedbinary, corrlimits)
  box off;
  print('clusterbinary', '-depsc2', '-painters')

  %% rearrange and also aggregate df over f

  %idx is original clusters
  %rearrange is new order

  dfoverf_bycluster = [];
  renamed_idx = [];
  for i = 1:length(idx)
      renamed_idx(i,1) = find(rearrange==idx(i));

  end

  for i = 1:length(rearrange)
      dfoverf_bycluster(i,:) = mean(sorteddfoverf(find(sorted_idx==rearrange(i)), :), 1);
  end


  %% SPATIAL DISTRIBUTION PER FISH
  % NEED TO CHANGE TO REFLECT NEW WAY OF REPRESENTING XCORR STIMREG
  markersize = 15;
  redthreshold = 0.4; %threshold for determining if cell is red. default is 0.65.

  if strcmp(celltype,'OB')
      corrlimits = [-0.6 0.6];
  else
      corrlimits = [-0.6 0.6];
  end

  range = corrlimits(2)-corrlimits(1);

  stimcatarray = {};
  stimcatarray{1} = 'water';
  stimcatarray{2} = 'kin';
  stimcatarray{3} = 'nonkin';
  stimcatarray{4} = 'adult';
  stimcatarray{5} = 'para';

  condition = 'water';
  xcorr_stimreg = dataarray.(condition).calcium.xcorr_stimreg;
  xcorr_motorstimreg = dataarray.(condition).calcium.xcorr_motorstimreg;

  %to replace non-signficant values with zeros
  % pcutoff = 0.05;         %pvalue cutoff
  % for s = 1:length(stimchannel)
  %     cat = stimcatarray{s};
  %     xcorr_stimreg.(cat)(xcorr_stimreg_pvalues.(cat)>=pcutoff) = 0;
  % end

  xpos = dataarray.(condition).calcium.xcoordinates;
  ypos = dataarray.(condition).calcium.ycoordinates;
  zpos = dataarray.(condition).calcium.zindices;
  redcells = dataarray.(condition).calcium.redcells;
  oxycells = dataarray.(condition).calcium.oxycells;

  %to generate a color map red to blue
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
  colorindices = [];

  figure(12)
  clf;
  figure(13);
  clf;
  figure(14);
  clf;
  figure(15);
  clf;

  for s = 1:size(activity_combined,2)  %number of categories
      cat = stimcatarray{s};
      corrcoeffs = xcorr_stimreg.(cat);

      for i = 1:length(corrcoeffs)
          coeff = corrcoeffs(i);

          colorindex = round((coeff - corrlimits(1))/range*100);
          colorindices(i,s) = colorindex;

          if colorindex <= 0
              color = cmap_fine(1,:);
          elseif colorindex > length(cmap_fine)
              color = cmap_fine(length(cmap_fine), :)
          else
              color = cmap_fine(colorindex,:);
          end

          figure(12)
          set(gca, 'YDir', 'reverse')
          view(viewxy)


          z = zpos(i);
          subplot(1,4,s);

          scatter3(xpos(i), ypos(i), zpos(i), markersize , color, 'filled');
          hold on;
          title(stimcatarray{s});

          if contains(pathname, 'GABAoxy')

              xlim([0 300]);
              ylim([0 300]);
          else
              xlim([0 xrange]);
              ylim([0 xrange]);

          end

          if strcmp(celltype, 'GABAoxy') || strcmp (celltype, 'GABAonly')

              try
                  if redcells(i) >= redthreshold
                      scatter3(xpos(i), ypos(i),zpos(i), markersize , 'MarkerFaceColor', 'none', 'MarkerEdgeColor', [0.5 0.5 0.5]);
                  end
              catch
              end

              try
                  if oxycells(i) == 1
                      scatter3(xpos(i), ypos(i), zpos(i), markersize , 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'y');
                  end
              catch
              end

          end


      end

  end

  isred=zeros(length(redcells),1);
  isred(find(redcells>=redthreshold),1)=1;

  figure(12)
  for s = 1:1:size(activity_combined,2)
      subplot(1,4,s)
      view(viewxy)
      xlim([0 xrange]);
      ylim([0 xrange]);
      set(gca, 'Ydir', 'reverse');
  end

  figure(12)
  print('3Dspatial_xcorr_allstims', '-depsc2', '-painters')


  figure(12)
  print(strcat('spatial', stimcatarray{1}), '-depsc2', '-painters')
  figure(13)
  print(strcat('spatial', stimcatarray{2}),  '-depsc2','-painters')
  figure(14)
  print(strcat('spatial', stimcatarray{3}), '-depsc2', '-painters')
  figure(15)
  print(strcat('spatial', stimcatarray{4}), '-depsc2', '-painters')

  %%
  save(sprintf('%s-%s%s', 'analyzeddata', celltype, '.mat'));