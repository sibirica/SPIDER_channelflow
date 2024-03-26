%
% Turbmat - a Matlab library for the JHU Turbulence Database Cluster
%   
% Get*.m, part of Turbmat
%

%
% This file is part of Turbmat.
% 
% Turbmat is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free
% Software Foundation, either version 3 of the License, or (at your option)
% any later version.
% 
% Turbmat is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
% 
% You should have received a copy of the GNU General Public License along
% with Turbmat.  If not, see <http://www.gnu.org/licenses/>.
%


function GetThresholdResult = GetThreshold(obj,authToken,dataset,field, ...
                                           time,threshold,spatialInterpolation, ...
                                           x_start,y_start,z_start,x_end,y_end,z_end)
%GetThreshold(obj,authToken,dataset,field,time,threshold,spatialInterpolation,X,Y,Z,Xwidth,Ywidth,Zwidth)
%
%   Retrieve the grid locations with norms for the specified field above the given threshold for the given time.
%   
%     Input:
%       authToken = (string)
%       dataset = (string)
%       field = (string)
%       time = (float)
%       threshold = (float)
%       spatialInterpolation = (SpatialInterpolation)
%       x_start = (int)
%       y_start = (int)
%       z_start = (int)
%       x_end = (int)
%       x_end = (int)
%       x_end = (int)
%   
%     Output:
%       GetThresholdResult = (ArrayOfThresholdInfo)

% Build up the argument lists.
data.val.authToken.name = 'authToken';
data.val.authToken.val = authToken;
data.val.authToken.type = '{http://www.w3.org/2001/XMLSchema}string';
data.val.dataset.name = 'dataset';
data.val.dataset.val = dataset;
data.val.dataset.type = '{http://www.w3.org/2001/XMLSchema}string';
data.val.field.name = 'field';
data.val.field.val = field;
data.val.field.type = '{http://www.w3.org/2001/XMLSchema}string';
data.val.time.name = 'time';
data.val.time.val = time;
data.val.time.type = '{http://www.w3.org/2001/XMLSchema}float';
data.val.threshold.name = 'threshold';
data.val.threshold.val = threshold;
data.val.threshold.type = '{http://www.w3.org/2001/XMLSchema}float';
data.val.spatialInterpolation.name = 'spatialInterpolation';
data.val.spatialInterpolation.val = spatialInterpolation;
data.val.spatialInterpolation.type = '{http://turbulence.pha.jhu.edu/}SpatialInterpolation';
data.val.x_start.name = 'x_start';
data.val.x_start.val = x_start;
data.val.x_start.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.y_start.name = 'y_start';
data.val.y_start.val = y_start;
data.val.y_start.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.z_start.name = 'z_start';
data.val.z_start.val = z_start;
data.val.z_start.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.x_end.name = 'x_end';
data.val.x_end.val = x_end;
data.val.x_end.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.y_end.name = 'y_end';
data.val.y_end.val = y_end;
data.val.y_end.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.z_end.name = 'z_end';
data.val.z_end.val = z_end;
data.val.z_end.type = '{http://www.w3.org/2001/XMLSchema}int';

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://turbulence.pha.jhu.edu/', ...
    'GetThreshold', ...
    data,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://turbulence.pha.jhu.edu/GetThreshold', ...
    soapMessage);
GetThresholdResult = parseSoapResponse(response);

% Fault message handling
if isfield(GetThresholdResult, 'faultstring')
    error('faultcode: %s\nfaultstring: %s\n', ...
        GetThresholdResult.faultcode, ...
        GetThresholdResult.faultstring);
end
