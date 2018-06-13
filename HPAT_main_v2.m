% Copyright 2013, 2014, 2015, 2016 Thomas M. Mosier, Kendra V. Sharp, and 
% David F. Hill
% This file is part of multiple packages, including the GlobalClimateData 
% Downscaling Package, the Hydropower Potential Assessment Tool, and the 
% Conceptual Cryosphere Hydrology Framework.
% 
% The above named packages are free software: you can 
% redistribute it and/or modify it under the terms of the GNU General 
% Public License as published by the Free Software Foundation, either 
% version 3 of the License, or (at your option) any later version.
% 
% The above named packages are distributed in the hope that it will be 
% useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with the Downscaling Package.  If not, see 
% <http://www.gnu.org/licenses/>.


%HPAT main script

%Implements simple degree index hydrologic model and uses output to assess 
%hydropower





%%USER INPUTS:
%Select the name for the region (used for naming output files):
region = {'USGS_14158790'}; %'OR_large_test';

%Select the mode to run the model in:
runType = 'simulate';    %Either 'default' (guess a parameter set), 
                            %'calibrate' (for optimizing parameters), 
                            %'calibrate_resume' (use this if a calibration 
                                %routine was interuptted)
                            %'validate' (apply optimized parameter set and 
                                %compare to observation data)
                            %'simulate' (apply optimized parameter set but
                                %do not compare to observations)
                  
nGage = 1;  %Number of gage files to load (if using 'CCHF_gage-*.txt' compiled gagedate written in 
            %previous run, set nGage to 1); 
             
%Use a shapefile to delineate ice extent?
iceGrid = 'none'; %Can be 'none' or 'same'
            
%Select the time-series elements to run the mdoel for:
startDate = [1981, 10];  %Vector with format '[year, month]' specifying date 
                        %to begin model run (run will include this date).
                        
endDate = [2010, 9];     %Vector with format '[year, month]' specifying date 
                        %to end model run (run will include this date).

monthsSpin = '12 months';	%String defining number of months to 
                            %run model for prior to start date
                        %'startDate'.
                        
timeStep = 'monthly';     %String specifying time resolution.  Can be 
                        %'daily', 'monthly', 'hourly'.
                        
dataRes = 'monthly';    %String specifying resolution of input time-series 
                        %data 
                        
printHydro = 'asc';     %Print run-of-river hydropower potential information in ESRI ASCII ('asc')
     
%THIS IS REQUIRED FOR RUN-OF-RIVER HYDROPOWER POTENTIAL CALCULATIONS
output = {'flow', {'all'}}; %surface flowrate through cell

                   




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%IMPLEMENTATION CODE (DO NOT EDIT):    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%
%%%%%HPAT SPECIFIC
%%%%%
%Select hydro modeling formulation:
moduleSet = hydro_modules();

%Options from the Conceptual Cryosphere Hydrology Model:
writeType = 'ascii';
writeAll = 0;

useDebris = 0;
blClip = 0;
optType = 'hybrid';  
fitType = 'kge_Parajka';
nRuns = 8000;
methodSiteCombine = 'square';
nStage = 1; 
category = {'input','cryo', 'land','routing'};
useFdr = 1;
dateEval = endDate;
%%%%%
%%%%%HPAT SPECIFIC
%%%%%

%%%%
%%%% SAME AS CCHF:
%%%%
if iscell(region)
    nSites = numel(region(:));
elseif ischar(region)
    nSites = 1;
    region = {region};
else
    error('CCHF_main:nmRegions',['The model dies not recognize the '...
        'format of the region input variable.']);
end

% if nSites ~= numel(regions(:))
%    error('CCHF_main:diffNumRegions',['The model is set to evaluate for ']) 
% end


%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
cd(pathScript);
addpath(genpath(pathScript));

%Initialize structure arrays for passing variables to downscaling script:
clear global sHydro sAtm sCryo sLand %Possibly simply initialize instead of clearing
% global sHydro{ii}
pathGage = cell(nSites, 1);
sPath = cell(nSites, 1);
sHydro = cell(nSites, 1);
sObs = cell(nSites, 1);
sMeta = struct;
    sMeta.('runType')  = runType;
    sMeta.('dateStart')= startDate; 
    sMeta.('dateEnd')  = endDate;
    sMeta.('spinup')   = monthsSpin;
    sMeta.('dt')       = timeStep;
    sMeta.('dataTsRes')= dataRes;
    sMeta.('region')   = region;
    sMeta.('module')   = moduleSet;
    sMeta.('mode')     = blanks(0);
    sMeta.('coef')     = cell(1,2);
    sMeta.('wrtTyp')   = writeType;
    sMeta.('wrtAll')   = writeAll;
    sMeta.('output')   = cell(nSites, 1);
    sMeta.('iceGrid')  = iceGrid;
    sMeta.('rtDir')    = cell(nSites, 1);
    sMeta.varLd        = {'pr','tas','tasmin','tasmax'};

    
sOpt = struct;
    sOpt.('fitTest') = fitType;
    sOpt.('methodSiteCombine') = methodSiteCombine;

if regexpbl(sMeta.runType,'calibrate')
    localCluster = parcluster('local'); %Check number of possible "workers"
    sOpt.maxWorker = localCluster.NumWorkers;
    
    sOpt.('type')    = optType; 

    sOpt.nbest = 100; %Record best _n_ performing parameter sets
    
    %Determine number of generations based on model runs and a
    %fixed population of 30
    nPop = 30;
    nGen = round(nRuns / nPop);

    if nGen < 1
        nGen = 1;
    end
    if nGen < 20
       warning('CCHF_main:nGenCalibrate',['Currently the number of '...
           'generations is set to ' num2str(nGen) ...
           '.  The recommended minimum is 20.']); 
    end

    sOpt.('att') = {...
        'generations', nGen; ...
        'population', nPop};
        
    %Set parameters for number of calibration stages:
    sOpt.nStage = nStage;
    if nStage > 1
        sOpt.stageVar = category;
    end
    
    if regexpbl(sOpt.type, 'uniform') && nStage > 1
        disp(['The number of calibration stages is being set to 1 '...
            'because the otpimization type is ' char(39) sOpt.type ...
            char(39) '.'])
    end
    
    %Set number of generations until stagnation:
    if regexpbl(sOpt.type,{'genetic','GA'})
        sOpt.stagnate = 20;
    elseif regexpbl(sOpt.type,{'PSO','hybrid'})
        sOpt.stagnate = 10;
    else
        sOpt.stagnate = 15;
    end
    
    sOpt.dateEval = dateEval;
    %If genetic algorithm used, assign specific fields:
	if regexpbl(sOpt.type,{'genetic','GA'})
        nMemCross = ceil(rateCross*(nPop - nElite));
        if mod(nMemCross, 2) %Ensure even number of parameter sets crossed
            nMemCross = nMemCross + 1;
        end
        nMemMut = nPop - nElite - nMemCross;
        sOpt.att(end+1:end+6,:) = {...
            'fitness requirement', perUse; ...
            'crossover pop', nMemCross; ...
            'crossover sites', nCross; ...
            'mutation pop', nMemMut; ...
            'mutation sites', siteMutate; ...
            'elite pop', nElite};
	end
end


%Create time vector to loop over based on start date and time-step
%If this field populated here, it wont be populated inside each model call,
%which saves time
sMeta = dates_run(sMeta, 'spin');
sMeta.progress = 'year';

%Load global constants from funtion:
sMeta.global = global_params; %contains global parameter values (albedo of ice, etc.)


%LOCATE INPUTS FOR EACH SITE:
startPath = pwd;
%Loop over sites:
for ii = 1 : nSites
    %Initialize structure for current site
    sPath{ii} = struct;
        sPath{ii}.('dem')    = blanks(0);
        sPath{ii}.('pr')     = blanks(0);
        sPath{ii}.('tas')    = blanks(0);
        sPath{ii}.('tasmin') = blanks(0);
        sPath{ii}.('tasmax') = blanks(0);
        sPath{ii}.('regclip')= blanks(0);
        sPath{ii}.('output') = blanks(0);
        sPath{ii}.('coef')   = blanks(0);
    sHydro{ii} = struct;
    
    %Ask to load information from directory of calibration run that was
    %interrupted
    if ii == 1 && regexpbl(sMeta.runType,{'calib','resume'},'and')
        uiwait(msgbox(sprintf(['Select the folder containing containing '...
            'the unfinished calibration for ' sMeta.region{ii} '.\n']), ...
            '(Click OK to Proceed)','modal'));
        sPath{ii}.resume = uigetdir(startPath,['Select the folder containing '...
            'the unfinished calibration for ' sMeta.region{ii}]);
    end
    
    %Digital Elevation Model (DEM) selection and display:
    uiwait(msgbox(sprintf(['Select the Digital Elevation Model (DEM) for ' ...
        sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
    [fileDem, foldDem] = uigetfile({'*.asc';'*.txt'},['Select the Digital Elevation Model '...
        '(DEM) for ' sMeta.region{ii}], startPath);
    sPath{ii}.dem = fullfile(foldDem, fileDem);
    disp([char(39) sPath{ii}.dem char(39) ' has been chosen as the DEM.']);

    %Update search path:
    [startPath, ~, ~] = fileparts(sPath{ii}.dem);

    if isempty(fileDem) || isequal(fileDem, 0)
       error('CCHF_main:noDEM','A DEM has not been selected. Therefore, the program is aborting.'); 
    end

    %Load manual flow direction grid (ESRI ASCII Format):
    if useFdr == 1
        %Flow direction selection and display:
        uiwait(msgbox(sprintf(['Select the flow direction grid (typically created using ArcGIS) for ' ...
            sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
        [fileFdr, foldFdr] = uigetfile({'*.asc';'*.txt'},['Select the flow direction grid for ' ...
            sMeta.region{ii}], startPath);
        sPath{ii}.fdr = fullfile(foldFdr, fileFdr);
        disp([char(39) sPath{ii}.fdr char(39) ' has been chosen as the flow direction grid.']);

        if isempty(fileFdr) || isequal(fileFdr, 0)
           error('CCHF_main:noFDR',['No flow direction grid '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end
    end


    %FIND GLACIER DATA PATH:
    if regexpbl(sMeta.iceGrid, 'fine') && isempty(find_att(sMeta.module, 'glacier','no_warning'))
        warning('CCHF_main:fineDemNoGlacier',['The glacier grid resolution '...
            'is being set to be the same as the main grid because no '...
            'glacier dynamics module has been selected.']);
        sMeta.iceGrid = 'same';
    end

    if regexpbl(sMeta.iceGrid, 'same')
        uiwait(msgbox(sprintf(['Select the glacier presence '...
            'geo-referenced grid or shapefile for ' ...
            sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
        [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the binary glacier presence grid for ' ...
            sMeta.region{ii}], startPath);
        sPath{ii}.ice = fullfile(foldGlac, fileGlac);
        disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the binary glacier presence grid.']);

        if isempty(fileGlac) || isequal(fileGlac, 0)
           error('CCHF_main:noGlacier',['No glacier presence grid '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end
    elseif regexpbl(sMeta.iceGrid, 'fine')
        iceAuxDisp = ['on a finer grid that will determine the ' ...
            'spatial scale of glacier process evaluation'];

        uiwait(msgbox(sprintf(['Select the glacier presence geo-referenced grid or shapefile for ' ...
            sMeta.region{ii} iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
        [fileGlac, foldGlac] = uigetfile({'*.*'},['Select the glacier presence geo-referenced grid or shapefile for ' ...
            sMeta.region{ii} iceAuxDisp '.\n'], startPath);
        sPath{ii}.ice = fullfile(foldGlac, fileGlac);
        disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the glacier presence grid.']);

        if isempty(fileGlac) || isequal(fileGlac, 0)
           error('CCHF_main:noGlacier',['No glacier presence grid '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end

        uiwait(msgbox(sprintf(['Select the glacier surface DEM for ' ...
            sMeta.region{ii} iceAuxDisp '.\n']), '(Click OK to Proceed)','modal'));
        [fileGlacDem, foldGlacDem] = uigetfile({'*.*'},['Select the glacier surface DEM for ' ...
            sMeta.region{ii} iceAuxDisp '.\n'], startPath);
        sPath{ii}.iceDem = fullfile(foldGlacDem, fileGlacDem);
        disp([char(39) sPath{ii}.ice char(39) ' has been chosen as the glacier DEM grid.']);

        if isempty(fileGlacDem) || isequal(fileGlacDem, 0)
           error('CCHF_main:noGlacierDem',['No glacier DEM grid '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end
    elseif ~regexpbl(sMeta.iceGrid, 'none')
        error('CCHF_main:unknownIceGridType',['The ice grid method is ' sMeta.iceGrid ', which is not known.']);
    end


    %Load debris cover grid (if single value, assumed uniform thickness):
    if useDebris && ~regexpbl(sMeta.iceGrid, 'none')
        %Flow direction selection and display:
        uiwait(msgbox(sprintf(['Select the debris cover grid for ' ...
            sMeta.region{ii} '.\n']), '(Click OK to Proceed)','modal'));
        [fileDeb, foldDeb] = uigetfile({'*.asc';'*.txt'},['Select the debris cover grid for ' ...
            sMeta.region{ii}], startPath);
        sPath{ii}.debris = fullfile(foldDeb, fileDeb);
        disp([char(39) sPath{ii}.debris char(39) ' has been chosen as the debris cover grid.']);

        if isempty(fileDeb) || isequal(fileDeb, 0)
           error('CCHF_main:noDebris',['No glacier surface debris grid '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end
    end


    %If a binary grid is being used to clip the region:
    if blClip == 1
        uiwait(msgbox(sprintf(['Select the binary region grid that will be used to clip the ' ...
            sMeta.region{ii} ' region.\n']), '(Click OK to Proceed)','modal'));
        [fileRegClip, foldRegClip] = uigetfile({'*.asc';'*.txt'},['Select the binary region clip grid for ' ...
            sMeta.region{ii}], startPath);
        sPath{ii}.regionClip = fullfile(foldRegClip, fileRegClip);
        disp([char(39) sPath{ii}.regionClip char(39) ' has been chosen as the binary region clip grid.']);

        if isempty(fileRegClip) || isequal(fileRegClip, 0)
           error('CCHF_main:noClip',['No region clipping file '...
               'has been selected, even though the option was chosen.' ...
               ' Therefore, the program is aborting.']); 
        end
    end


    %Identify climate variables to load:
    [sMeta.varLd, sMeta.varLdDisp] = clm_var_load(sMeta.module);
    %Load climate variables:
    sPath{ii} = clm_path_ui(sPath{ii}, sMeta, sMeta.region{ii});


    %Load DEM and calculate watershed geometry fields:
    sDem = read_geodata(sPath{ii}.dem, sMeta, 'none');
    sHydro{ii}.dem = sDem.data;
    sHydro{ii}.lat = sDem.lat;
    sHydro{ii}.lon = sDem.lon;
    if isfield(sPath{ii}, 'fdr')
        sFdr = read_geodata(sPath{ii}.fdr , sMeta, 'none');
        sHydro{ii}.fdrESRI = sFdr.data;
        clear('sFdr');
    else
        warning('cchfMain:noFDR',['CCHF currently does not have the '...
            'capacity to calculate a flow direction grid. Instead, '...
            'this must be generated in a standalone program such as ArcMap.']);
    end
    clear('sDem');
    
    %Load observation data to use for calibration or validation:
    if regexpbl(runType,{'calibrat','validat','default'})
        %Calibration/validation data required:
        uiwait(msgbox(sprintf(['Select the file(s) with observation data for ' ...
            sMeta.region{ii} '.  If data does not include information such as location, you will '...
            'be required to enter that in the main panel.']), ...
            '(Click OK to Proceed)','modal'));

        pathGage{ii} = cell(nGage,1);
        for jj = 1 : nGage
            [fileGage, foldGage] = uigetfile({'*.*'}, ...
                ['Select observation data ' num2str(jj) ' of ' ...
                num2str(nGage) ' for ' sMeta.region{ii}], startPath);
            pathGage{ii}{jj} = fullfile(foldGage, fileGage);
            
            disp([pathGage{ii}{jj} ' has been selected as observation data file ' num2str(jj) ' of ' num2str(nGage) ' for ' sMeta.region{ii} '.']);
        end
    end

    %Create unique 'main' output directory based upon inputs:
    if ii == 1 
        if isfield(sPath{ii},'resume')
%             sPath{ii}.output = sPath{ii}.resume;
            foldOutputMain = sPath{ii}.resume;
            [~, sMeta.strModule] = CCHF_out_dir(sPath{ii}, sMeta);
        else
            [foldOutputMain, sMeta.strModule] = CCHF_out_dir(sPath{ii}, sMeta);
        end
    end
    
    sPath{ii}.outputMain = foldOutputMain;
    sPath{ii}.output = fullfile(sPath{ii}.outputMain, sMeta.region{ii});
end %End of loop to select sites 


%display message for location of output folder:
if numel(nSites) > 1
   uiwait(msgbox(sprintf(['Model outputs will be written to a subfolder of ' ...
       sMeta.regions{1} ' because this region is listed first.\n']), ...
       '(Click OK to Proceed)','modal'));
end


%GET AND LOAD PARAMETER FILE (IF VALIDATION OR SIMULATION RUN)
if regexpbl(runType,{'valid','sim'})
    %Load parameters from file:
    uiwait(msgbox(sprintf(['Select the set of parameter coefficients ' ...
        'written during a calibration run of the CCHF model for ' ...
        sMeta.region{1} '.\n']), '(Click OK to Proceed)','modal'));
    [fileCoef, foldCoef] = uigetfile({'*.txt'; '*.csv'}, ['Select the set of parameter coefficients for ' ...
        sMeta.region{1}], startPath);
    sPath{ii}.coef = fullfile(foldCoef, fileCoef);
    disp([char(39) sPath{ii}.coef char(39) ' has been chosen as the set of parameter coefficients.']);
end


%LOAD ALL OBSERVATION DATA
for ii = 1 : nSites
    if regexpbl(runType,{'calibrat','validat','default'})
        %Remove any duplicate or empty elements:
        pathGage{ii} = unique(pathGage{ii}(~cellfun('isempty',pathGage{ii})));  

        %Load observation data:
        sObs{ii} = read_gagedata(pathGage{ii}, sHydro{ii}.lon, sHydro{ii}.lat, ...
            'time',[sMeta.dateStart;sMeta.dateEnd], ...
            'mask',sHydro{ii}.dem);
    else
        sObs{ii} = struct;
    end

    if isempty(sMeta.output{ii})
       sMeta.output{ii} =  output;
    end
    [sMeta.output{ii}, sObs{ii}] = output_all_gage(sMeta.output{ii}, sObs{ii}, sHydro{ii}.lon, sHydro{ii}.lat);
end


for ii = 1 : nSites
    if blClip == 1 && isfield(sPath{ii}, 'regionClip') 
        sRegClip = read_geodata(sPath{ii}.regionClip, sMeta, 'none');
        sHydro{ii}.regionClip = sRegClip.data;
        clear('sRegClip');
    end

    %Calculate fdr, fac, flow calculation order, slope, aspect, distance between neighboring centroids
    sHydro{ii} = watershed(sHydro{ii});

    %Initialize ice grid:
    [sHydro{ii}.sIceInit, sHydro{ii}.icbl] = ice_grid_init_loc_v2(sHydro{ii}, sPath{ii}, sMeta);
end


%%Find paremeters needed for current version of model:
%Do this by running CCHF model in 'parameter' mode
sMeta.mode = 'parameter';
sMeta.coefAtt = CCHF_engine_v4(sPath, sHydro, sMeta);


%IF VALIDATION, LOAD COEFFICIENT SET AND COMPARE TO PARAMETERS NEEDED FOR CURRENT MODEL:
if regexpbl(runType,{'valid','sim'})
    if regexpbl(sPath{ii}.coef,'.txt')
        [cfTemp, hdrParam] = read_CCHF_coef(sPath{ii}.coef);

        %Check that loaded coefficients are same as those needed for current
        %model:
        if ~regexpbl(cfTemp(:,1), sMeta.coefAtt(:,1),'and')
            %Find missing parameter:
            cfBl = zeros(numel(sMeta.coefAtt(:,1)),1);
            for ii = 1 : numel(sMeta.coefAtt(:,1))
                if regexpbl(cfTemp(:,1),sMeta.coefAtt{ii,1})
                    cfBl(ii) = 1;
                end
            end
            indCf = find(cfBl == 0);
            if ~isempty(indCf)
                nmQuit = -9999;
                uiwait(msgbox(sprintf(['The loaded parameter set does not '...
                    'include ' num2str(numel(indCf)) ' parameters needed '...
                    'for the ' char(39) sMeta.module char(39) ' model '...
                    'version.' char(10) 'You will be prompted to enter '...
                    'values for these parameters. If you wish to cancel '...
                    'this model run, enter ' num2str(nmQuit) ' for one of the values.']),...
                    '(Click OK to Proceed)','modal'));

                for ii = 1 : numel(indCf)
                   cfCurr = input(['Enter a value for ' char(39) ...
                       sMeta.coefAtt{indCf(ii),1} char(39) '. The bounds are ' ...
                       num2str(sMeta.coefAtt{indCf(ii),2}) ' thru ' ...
                       num2str(sMeta.coefAtt{indCf(ii),3}) ' and the parameter is '...
                       'used in ' char(39) sMeta.coefAtt{indCf(ii),5} char(39) ':' char(10)]); 
                   cfTemp(end+1,:) = {sMeta.coefAtt{indCf(ii),1}, cfCurr};
                   if cfCurr == nmQuit
                       error('CCHF_main:missingParam',['The exit command '...
                           'was entered during the process of entinering a'...
                           ' missing parameter value.']);
                   end
                end
            end
        end

        sMeta.coef = cfTemp;
    elseif regexpbl(sPath{ii}.coef,'.csv')
        prmArray = csvread(sPath{ii}.coef,1,0);
        fidPrm = fopen(sPath{ii}.coef,'r');
        varPrm = textscan(fidPrm,'%s',numel(prmArray(1,:)),'Delimiter',',');  
            varPrm = varPrm{1};
        fclose(fidPrm);
        
        if numel(varPrm(:)) - 1 == numel(sMeta.coefAtt(:,1))
            varPrm = varPrm{1}(1:end-1)'; %Last column is evaluation metric
            prmArray = prmArray(:,1:end-1);
        elseif numel(varPrm(:)) ~= numel(sMeta.coefAtt(:,1))
            error('CCHF_main:nParam',['The number of parameters loaded '...
                'is not equal to the number required by this model formulation.']);
        end
        
        %Check that all parameters present:
        if ~regexpbl(varPrm, sMeta.coefAtt(:,1))
            error('CCHF_main:missingParam',['There is a missing '...
                'parameter in the current list of parameter sets.']);
        end
        
        %Check that order of parameters is correct:
        for ii = 1 : numel(varPrm)
            if ~regexpbl(varPrm{ii}, sMeta.coefAtt{ii,1})
                error('CCHF_main:outoforderParam',[varPrm{ii} ' in the '...
                    'loaded parameter set is out of order relative to '...
                    'the parameter set identified for the current model '...
                    'formulation.']);
            end
        end
        
        %If loaded parameters contain multiple sets, put each in seperate
        %cell array and calculate correlation between parameters:
        if all(size(prmArray) ~= 1)
            prmCellTemp = cell(numel(prmArray(:,1)),1);
            for mm = 1: numel(prmCellTemp(:))
                prmCellTemp{mm} = prmArray(mm, :);
            end
            
            sMeta.coef = prmCellTemp;
        else
            sMeta.coef = prmArray;
        end
        
    else
        [~, ~, extCf] = fileparts(sPath{ii}.coef);
        error('CCHF_main:unknownCoefFormat',['The coeffienct file is in a ' ...
            char(39) extCf char(39) ' format, which has not been programmed for.']);
    end
elseif ~regexpbl(sMeta.runType,'calibrate') && numel(sMeta.coefAtt(1,:)) >= 6
    %Reduce dimensions of sMeta.coef (loaded during 'parameter' run)
   sMeta.coef = sMeta.coefAtt(:,[1,4,5]); 
end



%%START PROCESSING LOG AND DISPLAY RELEVENT CONTENT
%Record all messages displayed on screen to file in output directory
[fileDiary] = downscale_diary(sPath{1}.outputMain);
disp(['All displayed text is being written to ' char(39) fileDiary ...
    char(39) '.']);
%Display modules chosen:
disp('The chosen set of module representations is: ');
for ii = 1 : numel(sMeta.module(:,1))
    disp([sMeta.module{ii,1} ' = ' sMeta.module{ii,2}]);
end
disp(blanks(1));

%Display message with citation information:
[~] = MHS_cite('','CCHF');

%Display modeling package being used and user choices:
[~,dFold1,dFold2] = fileparts(pwd);
disp([char(39) dFold1 dFold2 char(39) ' is being used.' char(10)]);
disp_CCHF_meta_v2(sPath, sMeta)  



%Find individual files to load each timestep (saves time later)
for ii = 1 : nSites
    sPath{ii} = path_find_files(sPath{ii}, sMeta);    
    %Record root directory in sMeta (for use during simulations)
    indOutRt = regexpi(sPath{ii}.output,filesep);
    sMeta.rtDir{ii} = sPath{ii}.output(1:indOutRt(end)-1);
end


%%RUN THE HYDROLOGIC MODEL:
tStrt = now;
% sOutput = cell(1,1);
if regexpbl(sMeta.runType,'calib')
    sMeta.mode = 'calibrate';
    sMeta.coef = sMeta.coefAtt;
        sMeta = rmfield_x(sMeta,'progress');
    [sOutput, sMeta.coef] = CCHF_calibrate_v3(sObs, sOpt, sPath, sHydro, sMeta);
elseif regexpbl(sMeta.runType,{'valid','sim'})
    sMeta.mode = 'validate';
    
    if iscell(sMeta.coef) && ~ischar(sMeta.coef{1}) && ~all2d( size(sMeta.coef{1}) == 1 )
        sMeta = rmfield_x(sMeta,'progress');
        
        coefAll = nan(numel(sMeta.coef(:,1)), numel(sMeta.coef{1}));
        for ii = 1 : numel(coefAll(:,1))
            coefAll(ii,:) = sMeta.coef{ii,:};
        end
        
        %Open "Matlab pool" for parallel validation:
        if isempty(gcp('nocreate'))
            localCluster = parcluster('local'); %Check number of possible "workers"
            if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
                workers = sOpt.maxWorker;
            else
                workers = localCluster.NumWorkers;
            end
            parpool(workers); %Dedicate all available cores to parallel calibration
        end
        
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('off', 'all') %Turn off warnings during parfor
        end
        
        sOutput = cell(numel(coefAll(:,1)),1);
        nRuns = numel(coefAll(:,1));
        parfor ii = 1 : numel(coefAll(:,1))
            if mod(ii,10) == 0
               display(['Currently on validation run ' num2str(ii) ' of ' num2str(nRuns)]); 
            end
            %sMeta.coef = coefAll(ii,:);
            sOutput{ii} = CCHF_engine_v4(sPath, sHydro, sMeta, coefAll(ii,:)');
        end
        
        if ~isempty(gcp('nocreate'))
           pctRunOnAll warning('on', 'all') %Turn off warnings during parfor
        end
        
%         %Close dedicated workers used for parallel processing:
%         poolobj = gcp('nocreate');
%         delete(poolobj);
    else
        sOutput = CCHF_engine_v4(sPath, sHydro, sMeta);
    end
elseif regexpbl(sMeta.runType,'default')
    sMeta.mode = 'default';
    sMeta.coef = sMeta.coefAtt(:,[1,4]);
    sOutput = CCHF_engine_v4(sPath, sHydro, sMeta);
else
    error('CCHF_main:runType', [char(39) sMeta.runType char(39) ...
        ' was selected as the type of model run, but this option is not recognized.']);
end
tEnd = now;
disp(['It took ' num2str(round2((tEnd-tStrt)*24*60,1)) ' minutes for the current model run.']);


if ~regexpbl(sMeta.runType,'sim')
    warning('off','all'); %Turn off warning that some time-series not being used.
    %Define additional stats to calculate on data
    statsExtra = {'KGEr', 'KGEs', 'KGEb', 'NSE', 'MAPE', 'MAE'};
    if strcmpi(fitType,'kge')
        cellStats = [fitType, statsExtra];
    else
        cellStats = [fitType, 'kge', statsExtra];
    end
    
    
    if numel(sOutput) == nSites
        for mm = 1 : nSites
            %Create output directory for assessment plots:
            dirModObs = fullfile(sPath{mm}.output, 'mod_v_obs_plots');
            mkdir(dirModObs);
        
            %Display all stats and write to file
            report_stats_v2(sObs{mm}, sOutput{mm}, cellStats, sPath{mm}.output, sMeta);

            %Create plots: modeled versus observed
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, fitType, 'plot', dirModObs, 'scatter','grid','combineType', 'lon', sHydro{mm}.lon, 'lat', sHydro{mm}.lat);
            mod_v_obs_v2(sObs{mm}, sOutput{mm}, fitType, 'plot', dirModObs,'combineType');
        end
    else %This is used when validation is run to assess equifinality
        
        %Open "Matlab pool" for parallel validation:
        if isempty(gcp('nocreate'))
            localCluster = parcluster('local'); %Check number of possible "workers"
            if isfield(sOpt,'maxWorker') && localCluster.NumWorkers > sOpt.maxWorker
                workers = sOpt.maxWorker;
            else
                workers = localCluster.NumWorkers;
            end
            parpool(workers); %Dedicate all available cores to parallel calibration
        end
        
        for mm = 1 : nSites
            %Loop in reverse order in order to write the best performing
            %calibration parameter set last
            scoreTemp = cell(numel(sOutput{mm}),1);
            [scoreTemp{1}, typeTemp] = report_stats(sObs{mm}, sOutput{mm}{1}, cellStats, sPath.output, sMeta);

            scoreOut = cell(numel(scoreTemp{1}{1}),1);
            [scoreOut{:}] = deal(nan(numel(sOutput{mm}(:,1)),numel(cellStats)));

            if ~isempty(gcp('nocreate'))
               pctRunOnAll warning('off', 'all') %Turn off warnings during parfor
            end

            parfor nn = 2 : numel(sOutput{mm})
                [scoreTemp{nn}, ~] = report_stats(sObs{mm}, sOutput{mm}{nn}, cellStats, sPath{mm}.output, sMeta, 'no_disp','no_write');
            end
        
        
            %Reformat statistics so that each obs type is in seperate cell
            %array
            for oo = 1 : numel(sOutput) %loop over model run
                %rows are model runs and columns are metrics
                for nn = 1 : numel(typeTemp{1}) %loop over observation type
                    for ll = 1 : numel(typeTemp) %loop over metric
                        scoreOut{nn}(oo,ll) = scoreTemp{oo}{ll}(nn);
                    end
                end
            end

            %Write all validation results to csv files (one for each obs type):
            %Create output directory for assessment plots:
            dirMult = fullfile(sPath{mm}.output, ['mod_stats_' num2str(numel(sOutput)) '_best_cal_runs']);
            mkdir(dirMult);

            nCol = numel(cellStats);
            headerFmt = repmat('%s,',1,nCol-1);
            numFmt = repmat('%f,',1,nCol-1);

            hdrStats = cellStats(:)';

            for nn = 1 : numel(scoreOut)
                pathStats = fullfile(dirMult, [typeTemp{1}{nn} '_obs.csv']);
                fStats = fopen(pathStats,'w+');

                fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
                for oo = 1 : numel(scoreOut{nn}(:,1))
                    fprintf(fStats, [numFmt,'%f\n'], scoreOut{nn}(oo,:));
                end
                fclose(fStats); 
            end
        
            %Calculate validation statistics and write to file:
            aggStats = {'best_cal_run','mean_val', 'high_val', 'low_val', 'median_val', 'mode_val', 'SD_val'};
            nCol = numel(cellStats) + 1;
            headerFmt = repmat('%s,',1,nCol-1);
            numFmt = repmat('%f,',1,nCol-2);

            hdrStats = [blanks(1), cellStats(:)'];

            for nn = 1 : numel(scoreOut)
                pathStats = fullfile(dirMult, [typeTemp{1}{nn} '_stats_4multruns.csv']);
                fStats = fopen(pathStats,'w+');

                fprintf(fStats, [headerFmt,'%s\n'], hdrStats{:});
                for oo = 1 : numel(aggStats)
                    switch aggStats{oo}
                        case 'best_cal_run'
                            statCurr = scoreOut{nn}(1,:);
                        case 'mean_val'
                            statCurr = mean(scoreOut{nn});
                        case 'high_val'
                            statCurr = max(scoreOut{nn});
                        case 'low_val'
                            statCurr = min(scoreOut{nn});
                        case 'median_val'
                            statCurr = median(scoreOut{nn});    
                        case 'mode_val'
                            statCurr = mode(scoreOut{nn});
                        case 'SD_val'
                            statCurr = std(scoreOut{nn});
                        otherwise
                            display(['The current loop is being skipped '...
                                'because ' aggStats{oo} ' is not a recognized case.']);
                            continue
                    end

                    fprintf(fStats, [numFmt,'%f\n'], statCurr);
                end
                fclose(fStats); 
            end
        end
    end
    
    if ~isempty(gcp('nocreate'))
       pctRunOnAll warning('on', 'all') %Turn off warnings during parfor
    end

    %Close dedicated workers used for parallel processing:
    poolobj = gcp('nocreate');
    delete(poolobj);
    warning('on','all'); %Turn off warning that some time-series not being used.
end



%Display message that processing complete and turn off diary:
if regexpbl(sMeta.runType,'calib')
    runTypeDisp = 'calibration';
elseif regexpbl(sMeta.runType,'valid')
    runTypeDisp = 'validation';
elseif regexpbl(sMeta.runType,'default')
    runTypeDisp = 'default parameter set';
elseif regexpbl(sMeta.runType,'sim')
    runTypeDisp = 'simulation';
else
    runTypeDisp = 'unknown';
end
regDisp = '';
for ii = 1 : nSites
   regDisp = [regDisp, sMeta.region{ii}];
   if ii ~= nSites
      regDisp = [regDisp '/']; 
   end
end
disp(['The ' runTypeDisp ' run(s) for ' regDisp ' have finished.']);
for ii = 1 : nSites
    disp(['Results have been written to ' char(39) sPath{ii}.output ...
        char(39)]);
end
diary off   %Stop recording command history.





%%%%%%%%%%%%%%%%%%%%%
%% HPAT SPECIFIC
%%%%%%%%%%%%%%%%%%%%%
%MAKE PLOTS OF MODEL OUTPUT
sOutputPlot = struct;
if iscell(sOutput) && numel(sOutput(:)) == 1
    flds = fieldnames(sOutput{1});
    if numel(flds) > 1 || ~regexpbl(flds{1},'all')
        sOutputPlot = rmfield(sOutput{1},'all');
    end
elseif isstruct(sOutput)
    flds = fieldnames(sOutput);
    if numel(flds) > 1 || ~regexpbl(flds{1},'all')
        sOutputPlot = rmfield(sOutput,'all');
    end
else
    warning('HPAT_main:unexpectedOutputType','No output plots will be produced because the output are of an unexpected type.')
end

%Make plot
if numel(fieldnames(sOutputPlot)) > 0 
    dirModPlots = fullfile(sPath.output, 'model_output_plots');
        mkdir(dirModPlots);
    pathOutRt = fullfile(dirModPlots,'model');
    plot_CCHF(sOutputPlot, {'flow'}, sMeta, [pathOutRt '_flow']);
end






%Do hydropower calculations only during validation or simulation modes:
if regexpbl(sMeta.runType,'calib')
    disp('Hydropower calculations will not be conducted because this is a calibration run.');
else
    for kk = 1 : nSites
        %Load/find flow data written to file during model runs:
        flow = nan([numel(sMeta.dateRun(:,1)), size(sHydro{kk}.dem)], 'single');
        if iscell(sOutput) && numel(sOutput(:)) == 1
            if isfield(sOutput{kk}.all, 'pathflow')
                dirFlow = sOutput{1}.all.pathflow;
            else
                dirFlow = '';
                if isfield(sOutput{kk}.all, 'flow')
                    flow = sOutput{kk}.all.flow;
                end
            end
        elseif isstruct(sOutput)
            if isfield(sOutput{kk}.all, 'pathflow')
                dirFlow = sOutput.all.pathflow;
            else
                dirFlow = '';
                if isfield(sOutput.all, 'flow')
                    flow = sOutput.all.flow;
                end
            end
        else
            warning('HPAT_main:unexpectedOutputType','No output plots will be produced because the output are of an unexpected type.');
        end

        if ~isempty(dirFlow) && all(all2d(flow))
            filesFlow = find_files(dirFlow,'nc');

            for ii = 1 : numel(sOutput{kk}.all.date(:,1))
                pathCurr = fullfile(dirFlow, filesFlow{ii});
                flow(ii,:,:) = ncread(pathCurr, 'flow');
            end
        end
        
        flowAvg = nanmean(flow,1);

        %CALCULATE HYDROPOWER POTENTIAL
        if iscell(sOutput) && numel(sOutput(:)) == 1
            [powerAvg, powerRho, powerSd, powerSm] = power_potential(sHydro{kk}.slopeFdr, sHydro{kk}.dlFdr, flow, sOutput{1}.all.date);
        elseif isstruct(sOutput)
            [powerAvg, powerRho, powerSd, powerSm] = power_potential(sHydro{kk}.slopeFdr, sHydro{kk}.dlFdr, flow, sOutput.all.date);
        else
            warning('HPAT_main:unexpectedOutputType','No output plots will be produced because the output are of an unexpected type.');
        end


        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%% DEFINE RANKING METRIC AND CALCULATE HERE
        %%%%%%%%%%%%%%%%%%%%%%%%
        powerQual = power_quality(powerRho, powerSm);
        %%%%%%%%%%%%%%%%%%%%%%%%
        %%%% NED OF RANKING METRIC
        %%%%%%%%%%%%%%%%%%%%%%%%


        %Write hydropower statistics to file:
        dirHydroStats = fullfile(sPath{kk}.output,'hydro_stats');
            mkdir(dirHydroStats);
        pathFlowAvg = fullfile(dirHydroStats,'flow_mean');
        pathPwrAvg = fullfile(dirHydroStats,'power_mean');
        pathPwrRho = fullfile(dirHydroStats,'power_density');
        pathPwrSd  = fullfile(dirHydroStats,'power_standard_deviation');
        pathPwrSm  = fullfile(dirHydroStats,'power_stability_metric');
        pathPwrQual = fullfile(dirHydroStats,'power_quality');

        if regexpbl(printHydro,{'nc', 'netcdf'})
            print_grid_NC([pathFlowAvg, '.nc'], flowAvg,      'flow_avg', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 1);
            print_grid_NC([pathPwrAvg, '.nc'], powerAvg,       'pwr_avg', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 1);
            print_grid_NC([pathPwrRho, '.nc'], powerRho,   'pwr_density', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 2);
            print_grid_NC([ pathPwrSd, '.nc'],  powerSd,        'pwr_sd', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 2);
            print_grid_NC([ pathPwrSm, '.nc'],  powerSm, 'pwr_stability', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 2);
            print_grid_NC([pathPwrQual, '.nc'], powerQual, 'pwr_quality', sHydro{kk}.lon, sHydro{kk}.lat, nan(1,3), 2);
        elseif regexpbl(printHydro, 'asc')
            hdrWrt = ESRI_hdr(sHydro{kk}.lon, sHydro{kk}.lat, 'corner');
            write_ESRI_v4(  flowAvg, hdrWrt, [pathFlowAvg, '.asc'], 1);
            write_ESRI_v4( powerAvg, hdrWrt, [ pathPwrAvg, '.asc'], 1);
            write_ESRI_v4( powerRho, hdrWrt, [ pathPwrRho, '.asc'], 2);
            write_ESRI_v4(  powerSd, hdrWrt, [  pathPwrSd, '.asc'], 2);
            write_ESRI_v4(  powerSm, hdrWrt, [  pathPwrSm, '.asc'], 2);
            write_ESRI_v4(powerQual, hdrWrt, [pathPwrQual, '.asc'], 2);
        end
    end
end