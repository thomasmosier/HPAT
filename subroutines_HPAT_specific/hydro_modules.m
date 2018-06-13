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

function moduleSet = hydro_modules()

moduleSet = { ...
    'heat',  'SDI'; ... %surface heat flux module representation
    'mass', 'step'; ... %snow energy and mass module representation
    'icmlt', 'ratio'; ... %ice melt module
    'runoff', 'bucket'; ... %runoff module representation
    'toa', 'DeWalle'; ... %top-of-atmosphere radiation module representation
    'trans', 'DeWalle'; ... %atmospheric transmissivity module representation
    'albedo', 'Pelli'; ... %albedo module representation
    'pet', 'Hammon'; ... %potential evapotranspiration module representation
    'timelag', 'Johnstone'; ... %timelag (of flow through each grid cell) module representation
    'flow', 'lumped'; ... %flow routing module
    'density', 'Liston'; ... %densification of snow (does not participate in firn to ice processes)
    'snlq', 'percent'; ... %snow holding capacity
    'sca', 'max'; ...
    'glacier0', 'shea'; ...
	'glacier', 'static'; ... %glacier dynamic processes
    };