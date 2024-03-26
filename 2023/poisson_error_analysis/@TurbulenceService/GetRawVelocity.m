%
% Turbmat - a Matlab library for the JHU Turbulence Database Cluster
%   
% Get*.m, part of Turbmat
%

%
% Written by:
%  
% Zhao Wu
% The Johns Hopkins University
% Department of Mechanical Engineering
% zhao.wu@jhu.edu
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

function GetRawVelocityResult = GetRawVelocity(obj,authToken,dataset,time,X,Y,Z,Xwidth,Ywidth,Zwidth)
%GetPressure(obj,authToken,dataset,time,spatialInterpolation,temporalInterpolation,points)
%
%   Spatially interpolate the pressure field at a number of points for a given time.
%   
%     Input:
%       authToken = (string)
%       dataset = (string)
%       time = (float)
%       spatialInterpolation = (SpatialInterpolation)
%       temporalInterpolation = (TemporalInterpolation)
%       points = (ArrayOfPoint3)
%   
%     Output:
%       GetPressureResult = (ArrayOfPressure)

% Build up the argument lists.
data.val.authToken.name = 'authToken';
data.val.authToken.val = authToken;
data.val.authToken.type = '{http://www.w3.org/2001/XMLSchema}string';
data.val.dataset.name = 'dataset';
data.val.dataset.val = dataset;
data.val.dataset.type = '{http://www.w3.org/2001/XMLSchema}string';
data.val.time.name = 'time';
data.val.time.val = time;
data.val.time.type = '{http://www.w3.org/2001/XMLSchema}float';
data.val.X.name = 'X';
data.val.X.val = X;
data.val.X.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.Y.name = 'Y';
data.val.Y.val = Y;
data.val.Y.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.Z.name = 'Z';
data.val.Z.val = Z;
data.val.Z.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.Xwidth.name = 'Xwidth';
data.val.Xwidth.val = Xwidth;
data.val.Xwidth.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.Ywidth.name = 'Ywidth';
data.val.Ywidth.val = Ywidth;
data.val.Ywidth.type = '{http://www.w3.org/2001/XMLSchema}int';
data.val.Zwidth.name = 'Zwidth';
data.val.Zwidth.val = Zwidth;
data.val.Zwidth.type = '{http://www.w3.org/2001/XMLSchema}int';

% Create the message, make the call, and convert the response into a variable.
soapMessage = createSoapMessage( ...
    'http://turbulence.pha.jhu.edu/', ...
    'GetRawVelocity', ...
    data,'document');
response = callSoapService( ...
    obj.endpoint, ...
    'http://turbulence.pha.jhu.edu/GetRawVelocity', ...
    soapMessage);
GetRawVelocityResult = parseSoapResponse(response);

% Fault message handling
if isfield(GetRawVelocityResult, 'faultstring')
    error('faultcode: %s\nfaultstring: %s\n', ...
        GetRawVelocityResult.faultcode, ...
        GetRawVelocityResult.faultstring);
end
