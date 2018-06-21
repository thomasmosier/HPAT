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

%If snowfall calculated, use these parameters:
snowMod = 'ramp';
snowPrm = [-7.9, 5.2];
                
                
%%IMPLEMENTATION CODE (DO NOT EDIT):     
%LOAD MODULES:
%Add path of modules that may be called.
pathScript = mfilename('fullpath');
[pathScript, ~, ~] = fileparts(pathScript);
cd(pathScript);
addpath(genpath(pathScript));

%Load message to user
uiwait(msgbox(sprintf(['You will be requested to load several sources of climate data. \n'...
    'Ensure that each variable each input type is in a separate sub-folder. \n' ...
    'In cases where the input consists of multiple files, select one file from the sub-folder. \n' ...
    'The program will detect the others.\n']), ...
        '(Click OK to Proceed)','modal'));

%Add time period to downscaling method string:
if ~regexpbl(dMethod,{'proj','his'})
    if regexpbl(tPeriod,'proj')
        dMethod = [dMethod '_proj'];
    elseif regexpbl(tPeriod,'his')
        dMethod = [dMethod '_his'];
    end
end


%Select climatology years based on time period and downscaling method:
if ~isempty(yrsClimMan)
	disp(['Years of climatology used are being set to ' num2str(yrsClimMan(1)) '-' num2str(yrsClimMan(2))])
    yrsClim = yrsClimMan;
else
    yrsClim = [1950,2000]; %These are the WorldClim years.
end


nR = 1;

%Loop over regions to check order of bounding coordinates:
for rr = 1 : nR
%Check that lat-lon bounding box entered correctly
if latN{rr} <= latS{rr}
    display(['The northern bounding latitude entered is ' num2str(latN{rr}) ...
        ' and the southern latitude is ' num2str(latN{rr}) '.' ]);
    if latN{rr} < latS{rr}
        lonX = input(['The bounding box currently has a northern ' ...
            'extent of ' num2str(latN{rr}) ' degrees and southern extent of ' ...
            num2str(latN{rr}) ' degrees.' 'Are these supposed to be opposite? Enter ' char(39) 'yes' char(39) ' or ' ...
            char(39) 'no' char(39) ':' char(10)],'s');
        if regexpbl(lonX,'y')
            lonTemp = latN{rr};
            latN{rr} = latS{rr};
            latS{rr} = lonTemp;
            clear('latTemp');
        elseif regexpbl(lonX,'n')
            error('down_main:latXNo','Incorrect latitudes entered.');
        else
            error('down_main:latXUnknown',['Unkown option' lonX 'entered.']);
        end
    elseif latS{rr} == latN{rr}
        error('down_main:latS{rr}ame',['The southern and northern bounding '...
            'latitudes cannot be the same.']);
    end
end
if lonE{rr} <= lonW{rr}
    display(['The western bounding longitude entered is ' num2str(lonW{rr}) ...
        ' and the eastern longitude is ' num2str(lonE{rr}) '.' ]);
    if lonE{rr} < lonW{rr}
        lonX = input(['The bounding box currently has a western ' ...
            'extent of ' num2str(lonW{rr}) ' degrees and eastern extent of ' ...
            num2str(lonE{rr}) ' degrees.' 'Are these supposed to be opposite? Enter' char(39) 'yes' char(39) ' or ' ...
            char(39) 'no' char(39) ':' char(10)],'s');
        if regexpbl(lonX,'y')
            lonTemp = lonW{rr};
            lonW{rr} = lonE{rr};
            lonE{rr} = lonTemp;
            clear('lonTemp');
        elseif regexpbl(lonX,'n')
            error('down_main:lonXNo','Incorrect longitudes entered.');
        else
            error('down_main:lonXUnknown',['Unkown option' lonX 'entered.']);
        end
    elseif lonE{rr} == lonW{rr}
        error('down_main:lonSame',['The western and eastern bounding '...
            'longitudes cannot be the same.']);
    end
end
end

%Ensure that months are in ascending order and only unique elements:
mnths = sort(mnths,'ascend');
mnths = unique(mnths);

%If only one year entry, assume this is the only year being downscaled.
if numel(yrs) == 1
    yrs = [yrs, yrs];
end

%If copula bias-correction being used, need directories for both 
%temperature and precipitation data:
if regexpbl(dMethod, {'copula','JBC'})
    if ~any(cellfun(@(x) strcmpi(x,'pre') | strcmpi(x,'tmp'), metVar))
        error('downscale_main:copulaVar',['Currently copula bias-'...
            'correction is only available for precipiation and temperature']);
    end
    if ~all(cellfun(@(x) strcmpi(x,'pre') | strcmpi(x,'tmp'), metVar))
        warning('downscale_main:copulaVar',['Copula bias-'...
            'correction may or may not work with the current selection '...
            'of variables.']);
    end
    if any(cellfun(@(x) strcmpi(x,'tmp'), metVar)) && ~any(cellfun(@(x) strcmpi(x,'pre'), metVar))
        metVar = [metVar, 'pre'];
    end
    if any(cellfun(@(x) strcmpi(x,'pre'), metVar)) && ~any(cellfun(@(x) strcmpi(x,'tmp'), metVar))
        metVar = [metVar, 'tmp'];
    end
end

nMod = 1;
%Initialize structure arrays for passing variables to downscaling script:
sPath = struct;
    sPath.('hisObs')   = cell(length(metVar),nMod);
    sPath.('hisMod')  = cell(length(metVar),nMod);
    sPath.('hrClim')  = cell(length(metVar),1);
    sPath.('projMod')  = cell(length(metVar),nMod);
    sPath.('regClip') = blanks(0);
    sPath.('output')  = cell(length(metVar),1);
sMeta = struct;
    sMeta.('metVar')     = metVar;
    sMeta.('metVarDisp') = var_full_name(metVar);
    sMeta.('yrsOut')     = yrs;
    sMeta.('mnths')      = mnths;
%     sMeta.('crd')        = [lonW; lonE; latS; latN];
    sMeta.('yrsClim')     = yrsClim;
    sMeta.('region')      = region;
    sMeta.('hisObs')      = cell(numel(metVar),1);
    sMeta.('hisObsDisp')  = cell(numel(metVar),1);
    sMeta.('hisMod')      = cell(numel(metVar),nMod);
    sMeta.('hisModDisp')  = cell(numel(metVar),nMod);
    if regexpbl(tPeriod,'pro')
        sMeta.('projMod')      = cell(numel(metVar),nMod);
        sMeta.('projModDisp')  = cell(numel(metVar),nMod);
    end
    sMeta.('currVar')     = blanks(0);
sDownscale = struct;
    sDownscale.('method')     = dMethod;
    sDownscale.('methodFull') = blanks(0);
    sDownscale.('intrp')      = interpMeth;
    sDownscale.('wrtData')    = writeOpts;
    sDownscale.('wrtTyp')     = writeType; 
    sDownscale.('shepard')    = shepard;



if regexpbl(boolClip,'y')      
    uiwait(msgbox(sprintf(['Select the ESRI formatted ASCII reference ' ...
        'file, which will be used to clip downscaled data. ' ...
        '- Points that are NaN in the reference file will be set to ' ...
        'NaN in the Delta Data. \n']), '(Click OK to Proceed)','modal'));
    [fileRefOR, pathRefOR, ~] = uigetfile('*.asc');
    sPath.regClip = fullfile(pathRefOR,fileRefOR);
    disp(['The file ' char(39) sPath.regClip char(39) ' will be used as a '...
        'reference for clipping all downscaled grids to.']);
else
    sPath.regClip = '';
end
disp(blanks(1));   %Creates line of space between entries on different subjects


%GUI TO LOCATE INPUT FOLDERS:
for ii = 1 : length(sMeta.metVar(:))
    if ii == 1
        startPathClim = pwd;
    else
        if ~isempty(sPath.hrClim{ii-1})
            [startPathClim, ~, ~] = fileparts(sPath.hrClim{ii-1});
            indSepClim = regexpi(startPathClim,filesep);
            startPathClim = startPathClim(1:indSepClim(end)-1);
        else
            startPathClim = pwd;
        end
    end
    
    if regexpbl(sDownscale.method, 'direct')
        %Climatology selection and display (can be either file or folder depending on method used):
        strHrMsg = 'Select the gridded file that the time-series will be resampled to.';
        strHrLd = strHrMsg;
    else
        %Climatology selection and display (can be either file or folder depending on method used):
        strHrMsg = ['Select the high-resolution reference ' ...
            'climatologies for ' sMeta.metVarDisp{ii} ' (e.g. file in the folder ' ...
            'containing only WorldClim binary data and headers).'];
        strHrLd = ['Select file in the directory containing high-res reference climatologies for ' ...
            sMeta.metVarDisp{ii} '.'];
    end
    
    uiwait(msgbox(sprintf([strHrMsg '\n']), ...
        '(Click OK to Proceed)','modal'));
    
    [fileHrClim, pathHrClim, ~] = uigetfile('*.*', ...
        strHrLd, startPathClim);

    if regexpbl(fileHrClim,{'nc','grb'})
       sPath.hrClim{ii} = fullfile(pathHrClim, fileHrClim);
    else
        sPath.hrClim{ii} = pathHrClim;
    end

    disp([sPath.hrClim{ii} ' has been chosen as the high-resolution ' ...
        'reference climatology source of ' sMeta.metVarDisp{ii} ' data.']);

    %Try to identify name of high-resolution reference climatology:
    if regexpbl(sPath.hrClim{ii},{'world','clim'},'and')

    elseif regexpbl(sPath.hrClim{ii},{'will','matsuura'})
        sMeta.hrClim{ii} = 'WM';
        sMeta.hrClim{ii} = 'Willmott and Matsuura';
    elseif regexpbl(sPath.hrClim{ii},{'CRU'}) || regexpbl(sPath.hrClim{ii},{'climate','research','unit'},'and')
        sMeta.hrClim{ii} = 'CRU';
        sMeta.hrClim{ii} = 'Climate Research Unit';
    elseif regexpbl(sPath.hrClim{ii},'GPCC') || regexpbl(sPath.hrClim{ii},{'global','precipitation','climatology'},'and')
        sMeta.hrClim{ii} = 'GPCC';
        sMeta.hrClimDisp{ii} = 'Global Precipitation Climatology Centre';
    else
        [~, sPath.hrClim{ii}, ~] = fileparts(sPath.hrClim{ii});
        indHG = regexpi(sMeta.hrClim{ii}, '_');
        sPath.hrClim{ii} = sPath.hrClim{ii}(indHG(2)+1:indHG(end)-1);
        sMeta.hrClimDisp{ii} = ['CMIP5 GCM ' sMeta.hrClim{ii}];

        if ~regexpbl(sMeta.hrClim{ii},{'Amon','historic'})
            sMeta.hrClim{ii} = '?';
            sMeta.hrClimDisp{ii} = 'unknown';
            warning('main:input_TS_type', ['The program was not able to '...
                'detect the name of the historic observation time-series loaded.']);
        end
    end
    
    %For time-series data to downscale, loop over models to be used.
    for mm = 1 : nMod
        %Update search path:
        if ~isempty(sPath.hrClim{ii})
            if ii == 1
                [strtPathHis, ~, ~] = fileparts(sPath.hrClim{ii});
                indPathTs = regexpi(strtPathHis, filesep);
                strtPathHis = strtPathHis(1:indPathTs(end-1)-1);
            else
                [strtPathHis, ~, ~] = fileparts(sPath.hisMod{ii-1,mm});
                indSepHis = regexpi(strtPathHis,filesep);
                strtPathHis = strtPathHis(1:indSepHis(end)-1);
            end
        else
            strtPathHis = pwd;
        end
        
        %If the time-series being downscaled is projection data, need historic
        %time-series for same dataset.
        if regexpbl(sDownscale.method,'proj')
            projSrc = 'CMIP 5 GCM projection run';

            uiwait(msgbox(sprintf(['Select the ' sMeta.metVarDisp{ii} ' ' projSrc ' (' num2str(mm) ' of ' num2str(nMod) ') '...
                 '..' char(10)]), '(Click OK to Proceed)','modal'));
            [fileProjMod, pathProjMod, ~] = uigetfile({'*.nc'; '*.asc'; '*.grb'}, ...
                ['Select ' projSrc ' ' num2str(mm) ' of ' num2str(nMod) ' for ' sMeta.metVarDisp{ii}], strtPathHis);

            if regexpbl(fileProjMod,{'nc','grb'})
               sPath.projMod{ii,mm} = fullfile(pathProjMod, fileProjMod);
            else
                sPath.projMod{ii,mm} = pathProjMod;
            end

            disp(['You have chosen to use the file ' char(39) ...
                sPath.projMod{ii,mm} char(39) ' as the historical run output.']);  


            %Try to identify name of historic time-series being downscaled:
            if regexpbl(sPath.projMod{ii,mm},{'will','matsuura'})
                sMeta.projMod{ii,mm} = 'WM';
                sMeta.projModDisp{ii.mm} = 'Willmott and Matsuura';
            elseif regexpbl(sPath.projMod{ii,mm},{'CRU'}) || regexpbl(sPath.projMod{ii,mm},{'climate','research','unit'},'and')
                sMeta.projMod{ii,mm} = 'CRU';
                sMeta.projModDisp{ii.mm} = 'Climate Research Unit';
            elseif regexpbl(sPath.projMod{ii,mm},'GPCC') || regexpbl(sPath.projMod{ii,mm},{'global','precipitation','climatology'},'and')
                sMeta.projMod{ii,mm} = 'GPCC';
                sMeta.projModDisp{ii.mm} = 'Global Precipitation Climatology Centre';
            else
                [~, sMeta.projMod{ii,mm}, ~] = fileparts(sPath.projMod{ii,mm});
                indHG = regexpi(sMeta.projMod{ii,mm}, '_');
                sMeta.projMod{ii,mm} = sMeta.projMod{ii,mm}(indHG(2)+1:indHG(end)-1);
                sMeta.projModDisp{ii,mm} = ['CMIP5 GCM ' sMeta.projMod{ii,mm}];

                if ~regexpbl(sMeta.projMod{ii,mm},{'rcp'})
                    sMeta.projMod{ii,mm} = '?';
                    sMeta.projModDisp{ii.mm} = 'unknown';
                    warning('main:input_TS_type', ['The program was not able to '...
                        'detect the name of the historic observation time-series loaded.']);
                end
            end
        end


        %Time-series selection and display:
        if regexpbl(tPeriod,'proj')
            hisModStr = 'historic run for the projection run being downscaled';
            
            %Update start path based on projecation dataset selected:
            strtPathHis = sPath.projMod{ii,mm};
        elseif regexpbl(tPeriod,'his')
            hisModStr = 'data to downscale';
        end

        uiwait(msgbox(sprintf(['Select the ' sMeta.metVarDisp{ii} ' '...
            hisModStr ' (' num2str(mm) ' of ' num2str(nMod) ')' '. \n']), '(Click OK to Proceed)','modal'));

        [fileHisMod, pathHisMod,~] = uigetfile({'*.nc'; '*.asc'; '*.grb'}, ...
            ['Select ' sMeta.metVarDisp{ii} ' ' hisModStr ' model ' num2str(mm) ' of ' num2str(nMod) '.'], strtPathHis);  

        if regexpbl(fileHisMod,{'nc','grb'})
           sPath.hisMod{ii,mm} = fullfile(pathHisMod, fileHisMod);
        else
            sPath.hisMod{ii,mm} = pathHisMod;
        end

        disp([sPath.hisMod{ii,mm} ' has been chosen as the ' ...
            sMeta.metVarDisp{ii} ' time-series to downscale.']);  


        %Try to identify name of historic time-series being downscaled:
        if regexpbl(sPath.hisMod{ii,mm},{'will','matsuura'})
            sMeta.hisMod{ii,mm} = 'WM';
            sMeta.hisModDisp{ii,mm} = 'Willmott and Matsuura';
        elseif regexpbl(sPath.hisMod{ii,mm}, {'CRU'}) || regexpbl(sPath.hisMod{ii,mm}, {'climate','research','unit'},'and')
            sMeta.hisMod{ii,mm} = 'CRU';
            sMeta.hisModDisp{ii,mm} = 'Climate Research Unit';
        elseif regexpbl(sPath.hisMod{ii,mm}, 'GPCC') || regexpbl(sPath.hisMod{ii,mm}, {'global','precipitation','climatology'},'and')
            sMeta.hisMod{ii,mm} = 'GPCC';
            sMeta.hisModDisp{ii,mm} = 'Global Precipitation Climatology Centre';
        else
            [~, sMeta.hisMod{ii,mm}, ~] = fileparts(sPath.hisMod{ii,mm});
            indHG = regexpi(sMeta.hisMod{ii,mm}, '_');
            sMeta.hisMod{ii,mm} = sMeta.hisMod{ii,mm}(indHG(2)+1:indHG(end)-1);
            sMeta.hisModDisp{ii,mm} = ['CMIP5 GCM ' sMeta.hisMod{ii,mm}];

            if ~regexpbl(sMeta.hisMod{ii,mm}, {'Amon','historic'})
                sMeta.hisMod{ii,mm} = '?';
                sMeta.hisModDisp{ii,mm} = 'unknown';
                warning('main:input_TS_type', ['The program was not able to '...
                    'detect the name of the historic observation time-series loaded.']);
            end
        end
    end
    
    
    %If bias correcting, need 'observation' dataset to bias-correct against:
    if regexpbl(sDownscale.method,{'QM','Piani','JBC'})
        %Update search path:
        if ii == 1
            indPathNC = regexpi(strtPathHis, filesep);
            strtPathHis = strtPathHis(1:indPathTs(end-1)-1);
        else
            [strtPathHis, ~, ~] = fileparts(sPath.hisObs{ii-1});
            indPathNC = regexpi(strtPathHis, filesep);
            strtPathHis = strtPathHis(1:indPathNC(end)-1);
        end
        
        uiwait(msgbox(sprintf(['Select the ' sMeta.metVarDisp{ii} ...
            ' dataset to be used for bias-correcting the data being '...
            'downscaled.' char(10)]), '(Click OK to Proceed)','modal'));
        [fileHisObs, pathHisObs, ~] = uigetfile({'*.nc'; '*.grb'}, ...
            ['Select the ' sMeta.metVarDisp{ii} ...
            ' dataset to be used for bias-correction.'], strtPathHis);
        
        if regexpbl(fileHisObs,{'nc','grb'})
           sPath.hisObs{ii} = fullfile(pathHisObs, fileHisObs);
        else
            sPath.hisObs{ii} = pathHisObs;
        end
        
        disp(['You have chosen to use ' char(39) ...
            sPath.hisObs{ii} char(39) ' as the dataset to '...
            'bias-correct the data being downscaled to.']);  
    
        %Try to identify name of historic time-series for bias-correction:
        if regexpbl(sPath.hisObs{ii},{'will','matsuura'})
            sMeta.hisObs{ii} = 'WM';
            sMeta.hisObs{ii} = 'Willmott and Matsuura';
        elseif regexpbl(sPath.hisObs{ii},{'CRU'}) || regexpbl(sPath.hisObs{ii},{'climate','research','unit'},'and')
            sMeta.hisObs{ii} = 'CRU';
            sMeta.hisObsDisp{ii} = 'Climate Research Unit';
        elseif regexpbl(sPath.hisObs{ii},'GPCC') || regexpbl(sPath.hisObs{ii},{'global','precipitation','climatology'},'and')
            sMeta.hisObs{ii} = 'GPCC';
            sMeta.hisObsDisp{ii} = 'Global Precipitation Climatology Centre';
        else
            [~, sMeta.hisObs{ii}, ~] = fileparts(sPath.hisObs{ii});
            indHG = regexpi(sMeta.hisObs{ii}, '_');
            sMeta.hisObs{ii} = sMeta.hisObs{ii}(indHG(2)+1:indHG(end)-1);
            sMeta.hisObsDisp{ii} = ['CMIP5 GCM ' char(sMeta.hisObs{ii})];

            if ~regexpbl(sMeta.hisObs{ii},{'Amon','historic'})
                sMeta.hisObs{ii} = '?';
                sMeta.hisObsDisp{ii} = 'unknown';
                warning('main:input_TS_type', ['The program was not able to '...
                    'detect the name of the historic observation time-series loaded.']);
            end
        end 
    end
end

%Check if data selected seem to match meteological variable:
var_consistency_chk(sMeta,sPath)

%Create field for full name of the downscaling method being used
if regexpbl(sDownscale.method,'delta')
    sDownscale.methodDisp = 'Delta downscaling';
elseif regexpbl(sDownscale.method,'direct')
    sDownscale.methodDisp = 'Direct interpolation';
else
%     dMethDisp = 'An unknown ';
    error('main:downscale_type', ['The program is about to crash ' ...
        'because an unknown downscaling method has been selected.']);
end

if regexpbl(sDownscale.method,'QM')
    sDownscale.methodDisp = [sDownscale.methodDisp ' with univariate quantile mapping bias correction'];
elseif regexpbl(sDownscale.method,'Piani')
    sDownscale.methodDisp = [sDownscale.methodDisp ' with bivariate copula bias correction'];
elseif regexpbl(sDownscale.method,'JBC')
    sDownscale.methodDisp = [sDownscale.methodDisp ' with bivariate joint bias correction'];
end

%Deselect Shepard's bias correction if projected time-series being used
if regexpbl(sDownscale.method,'proj') && sDownscale.shepard == 1
   sDownscale.shepard = 0;
   warning('downscale_main:ShepardWeightFut',['Shepards weighting '...
       'has been unselected because it is only applicable to '...
       'downscaling historic data.']);
end


%%RUN DOWNSCALING AND STATISTICAL FUNCTIONS FOR EACH METEOROLOGICAL VARIABLE    
for mm = 1 : nMod
    sPathC = sPath;
    if isfield(sPathC,'hisMod') && ~isempty(sPath.hisMod{1,mm})
       sPathC.hisMod = sPath.hisMod(:,mm);  
    end
    if isfield(sPathC,'projMod') && ~isempty(sPath.projMod{1,mm})
       sPathC.projMod = sPath.projMod(:,mm);  
    end
    if isfield(sPathC,'hisObs') && ~isempty(sPath.projMod{1,1})
       sPathC.hisObs = sPath.hisObs(:,1);  
    end
    
    sMetaC = sMeta;
    if isfield(sMetaC,'hisMod') && ~isempty(sMeta.hisMod{1,mm})
       sMetaC.hisMod     =     sMeta.hisMod(:,mm);  
       sMetaC.hisModDisp = sMeta.hisModDisp(:,mm);
    end
    if isfield(sMetaC,'projMod') && ~isempty(sMeta.projMod{1,mm})
       sMetaC.projMod     =     sMeta.projMod(:,mm);  
       sMetaC.projModDisp = sMeta.projModDisp(:,mm);
    end
    
    for rr = 1 : nR
        sMetaC.crd = [lonW{rr}, lonE{rr}, latS{rr}, latN{rr}];
        sMetaC.region = region{rr};

        for jj = 1 : length(sMetaC.metVar(:))
            %Determine current variable to downscale:
            sMetaC.currVar = sMetaC.metVar{jj};

            %Create unique output directory based upon inputs:
            sPathC.output{jj} = downscale_out_dir(sPathC,sMetaC,sDownscale);

            %Record all messages displayed on screen to file in output directory
            [fileDiary] = downscale_diary(sPathC.output{jj});
            disp(['All displayed text is being written to ' char(39) fileDiary ...
                char(39) '.']);

            %Display message with citation information:
            [~] = MHS_cite(sDownscale.method,'downscale');

            %Display downscaling package being used and user choices:
            [~,dFold1,dFold2] = fileparts(pwd);
            disp([char(39) dFold1 dFold2 char(39) ' is being used.' char(10)]);
            disp_dowscale_meta(sPathC, sMetaC, sDownscale)  

            %Execute downscaling script:
            disp([upper(sMetaC.metVarDisp{jj}(1)) sMetaC.metVarDisp{jj}(2:end) ...
                ' grid processing has begun.']);
            downscale_v9(sPathC, sMetaC, sDownscale);

            %Only compare output data to GHCN station records if hindcast data
            %being produced.
            if ~isempty(ghcnStats) && regexpbl(sDownscale.method,'his')
                pathGridData = fullfile(sPathC.output{jj}, 'Downscaled_TS');

                for kk = 1 : numel(ghcnStats)
                    if regexpbl(ghcnStats{kk},{'adj','non'}) && ~regexpbl(ghcnStats{kk},{'joint'})
                        disp(['The downscaled grids are now being compared to available '...
                            'GHCN station data.']);
                        GHCN_2_grid_cmpr_v3(pathGridData, sMetaC, ghcnStats{kk}, 1, stnExclude);
                    end
                end
            end

            %Display message that processing complete and turn off diary:
            disp(['Finished processing ' sMetaC.region ' ' sMetaC.metVarDisp{jj} ' data.  ' ...
                'Output is in the directory ' char(39) sPathC.output{jj} char(39)'.' char(10)]);
            diary off   %Stop recording command history.
        end

        %Check if joint variable validation requested:
        if ~isempty(ghcnStats) && regexpbl(sDownscale.method,'his') && regexpbl(ghcnStats,'joint')
            if regexpbl(sMetaC.metVar,'pre') &&  regexpbl(sMetaC.metVar,'tmp')
        %         pathGHCNOut = regexpi(sPathC.output{jj},filesep);

                pathGridData = cell(2,1);
                pathGridData{1} = fullfile(sPathC.output{strcmpi(sMetaC.metVar,'tmp')}, 'Downscaled_TS');
                pathGridData{2} = fullfile(sPathC.output{strcmpi(sMetaC.metVar,'pre')}, 'Downscaled_TS');
                sMetaCGHCN = sMetaC;
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
            else
                warning('downscale_main:noJointGHCN',['Joint GHCN comparison '...
                    'is not being carried out because either precipitation or '...
                    'mean temprature were not downscaled.']);
            end
        end





        %%OPTIONAL: CALCULATE MONTHLY SNOWFALL
        if regexpbl(sDownscale.wrtData,'snowfall') && regexpbl(sMetaC.metVar,{'pre','tmp'},'and')
            for ii = 1 : numel(sMetaC.metVar)
               if regexpbl(sMetaC.metVar{ii},'pre')
                   indPre = ii;
               elseif regexpbl(sMetaC.metVar{ii},'tmp')
                   indTmp = ii;
               end
            end

            %Determine iterations based on expected downscaling outputs:
            if regexpbl(sDownscale.method,'proj')
                nTyp = 2;
            else
                nTyp = 1;
            end

            for ii = 1 : nTyp %Loop over original downscaled data and bias-corrected:
                if ii == 1 %Original data:
                    foldPre = fullfile(sPathC.output{indPre},'Downscaled_TS');
                    foldTmp = fullfile(sPathC.output{indTmp},'Downscaled_TS');
                elseif ii == 2 %bias-corrected data:
                    foldPre = fullfile(sPathC.output{indPre},'Downscaled_TS_BC');
                    foldTmp = fullfile(sPathC.output{indTmp},'Downscaled_TS_BC');
                end

                if regexpbl(sDownscale.wrtTyp,'asc')
                    dirOutSnow = snowfall(foldPre,foldTmp, snowMod, snowPrm, sDownscale.wrtTyp,1);
                elseif regexpi(sDownscale.wrtTyp,'cdf')
                    dirOutSnow = snowfall(foldPre,foldTmp, snowMod, snowPrm, sDownscale.wrtTyp,1,sMetaC.crd, sMetaC.yrsOut);
                else
                    error('downscale_main:snowDataTyp',[sDownscale.wrtTyp ...
                    ' is an unknown data type.  Add clause for this.']);
                end
            end
        elseif regexpbl(sDownscale.wrtData,'snowfall') && ~regexpbl(sMetaC.metVar,{'pre','tmp'},'and')
            warning('downscale_main:metVarNeed2',['Snowfall processing requires'...
            ' both precipitation and mean temperature be proccessed in the '...
            'same run.']);
        end
        
        %Close all files being used
        fclose all;
    end
end