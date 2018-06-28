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



%This script downscales a low-resolution (e.g. 0.5 degree) gridded monthly 
%time-series dataset with a high-resolution (e.g. 30 arc-second) set of
%monthly climatologies.  The code has been setup and tested with three
%low-res input datasets (Willmott & Matsuura, Climate Research Unit, and 
%Global Precipitation Climatology Centre) and one high-res climatology
%dataset (WorldClim), although it may also work with others.  If other
%high-res climatologies are used, the 'yrsClim' argument must be changed to
%reflect the years used in the climatologies.

%The first section, 'USER INPUTS' contains all of the options for the user
%to edit.  The second section 'IMPLEMENTATION CODE' should not be altered.
%Note: All strings are case insensitive.

clear all   %Clears variables in MATLAB, ensuring no conflicts.
clc         %Clears all existing text in the command window.

%%USER INPUTS:
%Select the time period to downscale:
tPeriod = 'his';   %Either 'his' (for historical downscaling) or 
                    %'project' (for future projections)

%Select the time-series elements to produce:
yrs = [1970, 2014]; %Ordered vector stating the first and last year of data 
                    %to be downscaled.
mnths = (1:12);  %List of all months to be used (e.g. if Jan, Feb, and Jun
                 %are to be downscaled, it should read 'mnths = [1,2,6];' 
                 %or 'mnths = (1:12);').

%Meteorological variable: 
metVar = {'pre','tmp','tmn','tmx'};   %Cell array of either/cobmination (Syntax for multiple is '{'pre', 'tmp'}'): 
                    %'pre' = monthly total precipitation
                    %'tmp' = monthly mean temperature
                    %'tmn' = monthly mean of minimum daily temperatures
                    %'tmx' = monthly mean of maximum daily temperatures

%Downscaling Method:
dMethod = 'delta';
                    %'QM_delta' = for downscaling GCM projection data 
                        %(implements empirical quantile mapping bias-correction and 
                        %Delta downscaling method).
                    %'JBC_delta' = for downscaling GCM projection data 
                        %(implements empirical joint bias-correction and 
                        %Delta downscaling method).
                    %'delta' = implements delta method.  This is
                        %recommended for hindcast downscaling.
                    %'direct' = directly interpolates low-resolution 
                        %time-series grid to reference climatology's 
       
                        
%REGION OPTIONS: Preset for Oregon, USA
%Select the name for the region (used for naming output files):
% region{1} = 'HPATtest';
% %Select the bounding box for the region of interest (in units of decimal 
% %degrees).  Coordinates of West and South correspond to negative values
% %here.  For example, '116 degrees West' should be input as '-116'.
% lonW{1} = -123;    %Western longitude defining region of interest
% lonE{1} = -121;    %Eastern longitude defining region of interest
% latS{1} = 43;      %Southern latitude defining region of interest
% latN{1} = 45;      %North

region{1} = 'USGS14158790';
%Select the bounding box for the region of interest (in units of decimal 
%degrees).  Coordinates of West and South correspond to negative values
%here.  For example, '116 degrees West' should be input as '-116'.
lonW{1} = -122.14172363281;    %Western longitude defining region of interest
lonE{1} = -122.02500;    %Eastern longitude defining region of interest
latS{1} =   44.324829101563;      %Southern latitude defining region of interest
latN{1} =   44.3892;      %North


%PRINTING OPTIONS: 
%Determine output file type:
writeType = 'ascii';    %'ascii' = individual ascii files (ESRI format)
                        %'netCDF' = A single netCDF file (Analogous to CMIP5 format)
%Determine which data to print to file ('writeOpts' is a cell array that 
%can contain any permutation of the below options): 
writeOpts = {'downscale_ts','input_ts','ref_clim','anomaly_lr','anomaly_hr'};   
%Options for data to write:
%'downscale_ts' - print high-resolution monthly downscaled data
%'input_clim' - print low-res climatology from time-series
%'input_ts' - print current low-res time-series element
%'downscale_clim' - print average of high-res downscaled data
%'anomaly_lr' - print low-res anomaly
%'anomaly_hr' - print high-res anomaly
%'ref_clim' - print high-res, clipped, Worldclim data
%'cdf' - print CDFs of all input data (only applicable when downscaling GCM data
%'downscale_seasonal' - print time-series downscaled output seasonsally averaged (in the
                %case of temperature) or integrated (in the case of 
                %precipitation) for all months being downscaled (i.e. if 
                %months = (1:12), these output are annual time-series)
%'snowfall' - option only available if both mean temperature and
                %precipitation downscaled in same package run.  %Calculates 
                %snow as any precipitation where surface air temperature 
                %for cell is at or below 0 degrees Celsius


%GEO-STAT CALCULATION (ONLY USED WHEN DOWNSCALING HISTORIC GRIDDED OBSERVATION DATA):   
ghcnStats = {'joint','adj_cdf','non_cdf'};
%ghcnStats = {}; 
                                %'adj' = calculate GHCN statistics using 
                                    %adjusted GHCN data
                                %'non' =  calculate GHCN statistics using 
                                    %unadjusted GHCN data
                                %'_cdf' = Add this extension to calculate 
                                    %stats in 1D probability space instead 
                                    %of temporal space.
                                %'_joint' = Add this extension to 
                                    %calculate stats in joint probability 
                                    %space (pre and tmp must both be 
                                    %simultaneously downscaled).
                                    
                            
                                    
%OVERWRITE DEFAULT SETTINGS:
yrsClimMan = [1960,1990]; %Format = [1st year, end year] of climatology. 
                    %By default script assumes WorldClim climatologies are 
                    %used and automatically assigns years [1950,2000].  In
                    %regions where PRISM available, this may be preferred
                    %and this variable must be set to years of PRISM
                    %climatology.
                    
%Anomaly interpolation method:
interpMeth = 'pchip';   %'pchip' (for Piecewise Cubic Hermite Interpolation 
                        %Polynomial) or 'spline' or 'linear'
stnExclude = []; %GHCN station numbers to exclude from
                            %statistical analysis.  Leave empty is no 
                            %stations need to be excluded.  
%Clip to non-rectangular region of interest using reference file:
boolClip = 'n'; %('y' = yes or 'n' = no).  If 'y', user prompted to 
                %locate reference file in ESRI ACII format.  All downscaled 
                %output will only be produced for cells in reference file 
                %where value is not NaN.

%Inverse Distance weighting (adjust downscaled historic data to match
%available GHCN stations):
%If selected, the downscaled grids will be compared to available GHCN
%station data for the region and time-series element, and a Shepard's
%Weighting bias correction method will be implemented.  Both the
%uncorrected and the corrected downscaled grids will be written to files.
shepard = 0;    %0 = Don't implement Shepard's bias correction
                %1 = Do implement (this adds significant time to the 
                %downscaling process) 

 
                
%%IMPLEMENTATION CODE (DO NOT EDIT):     
%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
cd(pathScript);
addpath(genpath(pathScript));

%Load message to user
uiwait(msgbox(sprintf(['You will be requested to load several sources of climate data. \n'...
    'Ensure that each variable each input type is in a separate sub-folder. \n\n' ...
    'In cases where the input consists of multiple files, select one file from the sub-folder. ' ...
    'The program will detect the others.\n']), ...
        '(Click OK to Proceed)','modal'));

% %Add time period to downscaling method string:
% if ~regexpbl(methodDs,{'proj','his'})
%     if regexpbl(simPeriod,'proj')
%         methodDs = [methodDs '_proj'];
%     elseif regexpbl(simPeriod,'his')
%         methodDs = [methodDs '_his'];
%     end
% end




%Put input coordinates into cell (this allows multiple regions to be
%downscaled
if isnumeric(lonDs)
   lonDs = {lonDs}; 
end
if isnumeric(latDs)
   latDs = {latDs}; 
end
if ischar(region)
   region = {region}; 
end

nR = numel(latDs(:));


%Ensure that months are in ascending order and only unique elements:
mnthsDs = sort(mnthsDs,'ascend');
mnthsDs = unique(mnthsDs);

%If only one year entry, assume this is the only year being downscaled.
if numel(yrsDs) == 1
    yrsDs = [yrsDs, yrsDs];
end

%Ensure downscaling variable is cell:
if ischar(metVarDs)
    metVarDs = {metVarDs};
end
nVarDs = numel(metVarDs);



%Initialize structure arrays for passing variables to downscaling script:
sPath = cell(nVarDs, 1);
    [ sPath{:} ] = deal(struct);
sMeta = struct;
    sMeta.('varOutDisp') = var_full_name(metVarDs);
    sMeta.('yrsOut')     = yrsDs;
    sMeta.('mnthsOut')   = mnthsDs;
    sMeta.('lonOut')     = lonDs;
    sMeta.('latOut')     = latDs;
    sMeta.('yrsBase')    = yrsClim;
    sMeta.('region')     = region;
    sMeta.('varOut')     = metVarDs;
    sMeta.('method')     = methodDs;
    sMeta.('methodDisp') = cell(nVarDs, 1);
    sMeta.('intrp')      = methodInterp;
    sMeta.('resample')   = methodResample;
    sMeta.('wrtData')    = writeOpts;
    sMeta.('wrtTyp')     = writeType; 
    sMeta.('shepard')    = shepard;
    sMeta.('period')     = simPeriod;
    sMeta.('nSim')       = nSim;
    sMeta.('timestep')   = dataTimestep;
    if exist('unitsUse', 'var')
        sMeta.('unitsUse')   = unitsUse;
    end

    
%Querry downscaling methods catalog to determine which inputs to load
dataLd = cell(nVarDs);
for ii = 1 : nVarDs
    [dataLd{ii}, sMeta.('methodDisp')] = ds_catalog(methodDs, metVarDs{ii}, simPeriod);
end
clear ii


%GUI TO LOCATE INPUT FOLDERS:
strUnknown = '';
uiInput = cell(nVarDs,1);
[ uiInput{:} ] = deal(cell(0,2));
for ii = 1 : nVarDs
    disp(['Selecting inputs for ' sMeta.('varOut'){ii}]);
    
    %Check if there are any user inputs (strings or numbers) to request:
    for jj = numel(dataLd{ii}(:,1)) : -1 : 1
        if regexpbl(dataLd{ii}{jj,2}, 'uiget')
            uiInput{ii}(end+1,:) = {[dataLd{ii}{jj,2} '_' dataLd{ii}{jj,1}],''};
            
            if regexpbl(dataLd{ii}{jj,1}, 'string')
                uiInput{ii}{end,2} = input(dataLd{ii}{jj,3},'s');
            elseif regexpbl(dataLd{ii}{jj,1}, 'number')
                uiInput{ii}{end,2} = input(dataLd{ii}{jj,3});
            else
                error('gcdMain:unknownUIRequest',['The UI input request for type ' ...
                    dataLd{ii}{jj,1} ' is not recognized and has not been programmed for.']);
            end
            
            dataLd{ii}(jj,:) = []; 
        end
    end
    clear jj
    
    %Loop over the required inputs for current downscaling variable:
    for jj = 1 : numel(dataLd{ii}(:,1))
        
        if regexpbl(dataLd{ii}{jj,1}, 'unknown') && isempty(strUnknown)
           strUnknown = input(['What is the variable to be used for ' ...
               'defining the bias (' char(39) 'tmp' char(39) ', ' ...
               char(39) 'tmx' char(39) ', ' char(39) 'tmn' char(39) ...
               ', or ' char(39) 'pre' char(39) ')?'], 's');
        end
        
        if regexpbl(dataLd{ii}{jj,1}, 'unknown')
            dataLd{ii}{jj,1} = strUnknown;
            fieldCurr = [dataLd{ii}{jj,2} '_' strUnknown];
        else
            fieldCurr = [dataLd{ii}{jj,2} '_' dataLd{ii}{jj,1}];
        end
        
        if regexpbl(fieldCurr, 'sim') %For projections, can have multiple independent simulations
            %Determine if any previous simulation fields loaded
            prevSim = nan;
            if jj ~= 1
                for mm = 1 : jj-1
                    if regexpbl(dataLd{ii}{mm,2}, 'sim')
                        prevSim = mm;
                    end
                end
                clear mm
            end
            
            sPath{ii}.(fieldCurr) = cell(nSim,1);
            
            %Find start path:
            for mm = 1 : nSim
                if ~isnan(prevSim)
                    fieldPrev = [dataLd{ii}{prevSim,2} '_' dataLd{ii}{prevSim,1}];

                    if isfield(sPath{ii}, fieldPrev) && ~isempty(sPath{ii}.(fieldPrev))
                        strtPath = sPath{ii}.(fieldPrev){mm};
                        
                        indSepClim = regexpi(strtPath,filesep);
                        strtPath = strtPath(1:indSepClim(end-1)-1);
                    else
                        strtPath = pwd;
                    end
                else
                    if mm ~= 1
                        strtPath = sPath{ii}.(fieldCurr){mm-1};
                        indSepClim = regexpi(strtPath,filesep);
                        strtPath = strtPath(1:indSepClim(end-2)-1);
                    elseif jj ~= 1
                        fieldPrev = [dataLd{ii}{jj-1,2} '_' dataLd{ii}{jj-1,1}];
                    
                        if isfield(sPath{ii}, fieldPrev) && ~isempty(sPath{ii}.(fieldPrev))
                            if iscell(sPath{ii}.(fieldPrev))
                                strtPath = sPath{ii}.(fieldPrev){mm};
                            else
                                strtPath = sPath{ii}.(fieldPrev);
                            end

                            indSepClim = regexpi(strtPath,filesep);
                            strtPath = strtPath(1:indSepClim(end-2)-1);
                        else
                            strtPath = pwd;
                        end
                    else
                        strtPath = pwd;
                    end
                end
                
                strLdCurr = ['Select a ' dataLd{ii}{jj,1} ' ' dataLd{ii}{jj,3} ' file (' num2str(mm) ' of ' num2str(nSim) ')'];

                uiwait(msgbox(sprintf([strLdCurr '\n']), ...
                    '(Click OK to Proceed)','modal'));
    
                [fileCurr, foldCurr, ~] = uigetfile('*.*', ...
                    strLdCurr, strtPath);
                if isequal(fileCurr, 0)
                    error('downscaleMain:noFileSelect', [dataLd{ii}{jj,3} ' not selected']);
                end
                sPath{ii}.(fieldCurr){mm} = fullfile(foldCurr, fileCurr);

                %Try to identify name of historic time-series for bias-correction:
                [nmCurr, nmCurrDisp] = clim_data_name(sPath{ii}.(fieldCurr){mm});
                sPath{ii}.(['nm_' fieldCurr]){mm} = {nmCurr, nmCurrDisp};

                disp([sPath{ii}.(fieldCurr){mm} ' has been chosen as the ' dataLd{ii}{jj,3} ...
                    ' (' num2str(mm) ' of ' num2str(nSim) '), which '... 
                    'has been identified as being from the ' nmCurrDisp ' dataset']);

            end %End of simulation selection loop
            clear mm
            
        else
            prevVar = nan;
            if jj ~= 1
                for mm = 1 : jj-1
                    if regexpbl(dataLd{ii}{mm,2},dataLd{ii}{jj,2}) && ~regexpbl(dataLd{ii}{mm,1},dataLd{ii}{jj,1}) 
                        prevVar = mm;
                    end
                end
                clear mm
            end
            
            if ~isnan(prevVar)
                fieldPrev = [dataLd{ii}{prevVar,2} '_' dataLd{ii}{prevVar,1}];

                if isfield(sPath{ii}, fieldPrev) && ~isempty(sPath{ii}.(fieldPrev))
                    strtPath = sPath{ii}.(fieldPrev);
                    
                    if iscell(strtPath)
                        strtPath = strtPath{1};
                    end

                    indSepClim = regexpi(strtPath,filesep);
                    strtPath = strtPath(1:indSepClim(end-1)-1);
                else
                    strtPath = pwd;
                end
            else
                if jj ~= 1
                    fieldPrev = [dataLd{ii}{jj-1,2} '_' dataLd{ii}{jj-1,1}];
                    
                    if isfield(sPath{ii}, fieldPrev) && ~isempty(sPath{ii}.(fieldPrev))
                        strtPath = sPath{ii}.(fieldPrev);
                        
                        if iscell(strtPath)
                            strtPath = strtPath{1};
                        end

                        indSepClim = regexpi(strtPath,filesep);
                        strtPath = strtPath(1:indSepClim(end-2)-1);
                    else
                        strtPath = pwd;
                    end
                else
                    strtPath = pwd; 
                end
            end
            
            %Create/update search path:
            if ii == 1 && jj == 1
                strtPath = pwd;
            elseif ii == 1 && jj ~= 1
                fieldPrev = [dataLd{ii}{jj-1,2} '_' dataLd{ii}{jj-1,1}];

                if isfield(sPath{ii}, fieldPrev) && ~isempty(sPath{ii}.(fieldPrev))
                    strtPath = sPath{ii}.(fieldPrev);
                    if iscell(strtPath)
                        strtPath = strtPath{1};
                    end
                    indSepClim = regexpi(strtPath,filesep);
                    strtPath = strtPath(1:indSepClim(end-2)-1);
                else
                    strtPath = pwd;
                end
            else
                fieldPrev = [dataLd{ii-1}{jj,2} '_' dataLd{ii-1}{jj,1}];

                if isfield(sPath{ii-1}, fieldPrev) && ~isempty(sPath{ii-1}.(fieldPrev))
                    strtPath = sPath{ii-1}.(fieldPrev);
                    if iscell(strtPath)
                        strtPath = strtPath{1};
                    end
                    indSepClim = regexpi(strtPath,filesep);
                    strtPath = strtPath(1:indSepClim(end-1)-1);
                else
                    strtPath = pwd;
                end
            end

            strLdCurr = ['Select a ' dataLd{ii}{jj,1} ' ' dataLd{ii}{jj,3} ' file'];

            uiwait(msgbox(sprintf([strLdCurr '\n']), ...
                    '(Click OK to Proceed)','modal'));
                
            [fileCurr, foldCurr, ~] = uigetfile('*.*', ...
                strLdCurr, strtPath);
            if isequal(fileCurr, 0)
                error('downscaleMain:noFileSelect', [dataLd{ii}{jj,3} ' not selected']);
            end
            sPath{ii}.(fieldCurr) = fullfile(foldCurr, fileCurr);

            %Try to identify name of historic time-series for bias-correction:
            [nmCurr, nmCurrDisp] = clim_data_name(sPath{ii}.(fieldCurr));
            sPath{ii}.(['nm_' fieldCurr]) = {nmCurr, nmCurrDisp};
            
            disp([sPath{ii}.(fieldCurr) ' has been chosen as the ' dataLd{ii}{jj,3} ...
                ', which has been identified as being from the ' nmCurrDisp ' dataset']);
        end
    end
    clear jj
    
%     %Check if data selected seem to match meteological variable:
%     var_consistency_chk(sMeta,sPath{ii})
end
clear ii

%Set any UI inputs:
sMeta.('uiInput') = uiInput;


%Deselect Shepard's bias correction if projected time-series being used
if regexpbl(sMeta.period,'proj') && sMeta.shepard == 1
   sMeta.shepard = 0;
   warning('downscale_main:ShepardWeightFut',['Shepards weighting '...
       'has been unselected because it is only applicable to '...
       'downscaling historic data.']);
end


%%RUN DOWNSCALING AND STATISTICAL FUNCTIONS FOR EACH METEOROLOGICAL VARIABLE   
%Loop over input models
for rr = 1 : nR
    sMeta.indRegion = rr;
    
    for mm = 1 : nSim
        sMeta.indSim = mm;
            
        for jj = 1 : nVarDs
            sMeta.indDs = jj;
            
            %Create output directory
            foldSubRt = char([sMeta.method '_' sMeta.intrp '_' sMeta.varOut{sMeta.indDs} '_' sMeta.region{sMeta.indRegion}]);
            sPath{jj}.output = ds_dir_out(sPath{jj}, sMeta.indSim, foldSubRt);
            
            %Record all messages displayed on screen to file in output directory
            [fileDiary] = command_log(sPath{jj}.output);
            disp(['All displayed text is being written to ' char(39) fileDiary ...
                char(39) '.']);

            %Display message with citation information:
            [~] = MHS_cite_v2('downscale');
            
            %Display downscaling package being used and user choices:
            [~,dFold1,dFold2] = fileparts(pwd);
            disp([char(39) dFold1 dFold2 char(39) ' is being used.' char(10)]);
            ds_meta_disp(sPath{jj}, dataLd{jj}, sMeta.indSim, sMeta, sMeta.varOut{sMeta.indDs})  

            %Execute downscaling script:
            disp([upper(sMeta.varOutDisp{jj}(1)) sMeta.varOutDisp{jj}(2:end) ...
                ' grid processing has begun.']);
            ds_engine_v10(sPath{jj}, sMeta);
            
            
            %Only compare output data to GHCN station records if historic data
            %being processed.
            if ~isempty(ghcnStats) && regexpbl(sMeta.period,'his')
                pathGridData = fullfile(sPath{jj}.output, 'output');

                for kk = 1 : numel(ghcnStats)
                    if regexpbl(ghcnStats{kk},{'adj','non'}) && ~regexpbl(ghcnStats{kk},{'joint'})
                        disp(['The downscaled grids are now being compared to available '...
                            'GHCN station data.']);
                        GHCN_2_grid_cmpr_v3(pathGridData, sMeta, ghcnStats{kk}, 1, stnExclude);
                    end
                end
            end

            %Display message that processing complete and turn off diary:
            disp(['Finished processing ' sMeta.region{sMeta.indRegion} ' ' sMeta.varOutDisp{jj} ' data.' char(10) ...
                'Output is in the directory ' char(39) sPath{jj}.output char(39)'.' char(10)]);
            diary off   %Stop recording command history.
        end
        clear jj
        
        
        %Check if joint variable validation requested:
        if ~isempty(ghcnStats) && regexpbl(sMeta.period,'his') && regexpbl(ghcnStats,'joint')
            if regexpbl(sMeta.varOut,'pre') &&  regexpbl(sMeta.varOut,'tmp')
        %         pathGHCNOut = regexpi(sPathC.output{jj},filesep);

                pathGridData = cell(2,1);
                pathGridData{1} = fullfile(sPath{strcmpi(sMeta.varOut,'tmp')}.output, 'output');
                pathGridData{2} = fullfile(sPath{strcmpi(sMeta.varOut,'pre')}.output, 'output');
                sMetaCGHCN = sMeta;
                sMetaCGHCN.metVar = {'tmp','pre'};
                disp(['The downscaled grids for precipitation and mean '...
                    'temperature are now being jointly compared to available '...
                    'GHCN station data.']);

                for kk = 1 : numel(ghcnStats)
                    if regexpbl(ghcnStats{kk},'joint')
%                     if regexpbl(ghcnStats{kk},{'adj','joint'},'and') || regexpbl(ghcnStats{kk},{'non','joint'},'and') 
                        GHCN_2_grid_cmpr_v3(pathGridData, sMetaCGHCN, ghcnStats{kk}, 1, stnExclude);
                    end
                end
                clear kk
            else
                warning('downscale_main:noJointGHCN',['Joint GHCN comparison '...
                    'is not being carried out because either precipitation or '...
                    'mean temperature were not downscaled.']);
            end
        end %End of joint variable check
        
    end
    clear mm
end
clear rr 