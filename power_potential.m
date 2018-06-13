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

function [powerAvg, powerRho, powerSd, sm] = power_potential(slope, dl, flow, date)
 
%Basic power equation is P = e * V_flow * dH * rho_w * g, where e (the efficiency) is
%assumed to be one here.
%Power density is calculated as P_rho = (e * V_flow * dH) / dl, where dl is the flow path
%between grid cell centroids

power = nan(size(flow),'single');

days = days_since([1900,1,1], date, 'gregorian');
nTstep = numel(days);

szGrid = size(slope);


% dDays = diff(days);

% if nanmean(dDays) > 0.9/24 && nanmean(dDays) < 1.1/24 %Assume hourly time step
%     dt = ones(nTstep,1,'single')/24;
% elseif nanmean(dDays) > 0.9 && nanmean(dDays) < 1.1 %Assume daily time step
%     dt = ones(nTstep,1,'single');
% elseif all(dDays > 26) && all(dDays > 32)  %Assume monthly time step
%     dt = nan(nTstep,1,'single');
%     for  ii = 1 : nTstep
%         dt(ii) = eomday(date(ii,1:2), 'gregorian');
%     end
% else %Not sure what the time step is...
%     error('power_potential:unknownTimeStep',['The average time step is ' ...
%         num2str(nanmean(dDays)) ' days, which has not been programmed for.']);
% end

%Calculate power distribution at each time step:
for ii = 1 : nTstep
    power(ii,:,:) = 1000*9.81*(abs(slope) .* squeeze(flow(ii,:,:))); 
%    power(ii,:,:) = 1000*9.81*dt(ii)*(dem .* squeeze(flow(ii,:,:))); 
end

%Calculate average daily power and standard deviation:
powerAvg = squeeze(nanmean(power, 1));
powerRho = powerAvg ./ dl;
powerSd  = squeeze(nanstd(power));

%Calculate stability metric (based on monthly statistics):
tsMonth = unique(date(:,1:2),'rows');
powerMnth = nan([numel(tsMonth(:,1)), szGrid],'single');
for ii = 1 : numel(tsMonth(:,1))
    indCurr = ismember(date(:,1:2), tsMonth(ii,:),'rows');
    powerMnth(ii,:,:) = nanmean(power(indCurr,:,:), 1); 
end

powerMnthAvg = nan([12, szGrid],'single');
powerMnthDev = nan([12, szGrid],'single');

for ii = 1 : 12
    %indCurr = find(tsMonth(:,2) == ii); 
    powerMnthAvg(ii,:,:) = squeeze(nanmean(power(tsMonth(:,2) == ii,:,:), 1)); 
    powerMnthDev(ii,:,:) = abs(squeeze(powerMnthAvg(ii,:,:)) - powerAvg);
end
powerMonthSd = squeeze(nanstd(powerMnthAvg));

d = zeros(szGrid,'single');

for ii = 1 : 12
    d = d + (squeeze(powerMnthAvg(ii,:,:)) > powerMonthSd);
end

sm = 1 - squeeze(nansum(powerMnthDev, 1)) ./ (d.*squeeze(abs(nanmax(powerMnthAvg, 1))) + 11*squeeze(nanmean(powerMnthAvg, 1)));
