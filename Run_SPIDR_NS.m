addpath turbmat-master/

% File save location
if exist("data", "dir")~=7
   mkdir("data") 
end

reg_opts = {}; % Regression options
reg_opts.domain_loc = "edge"; % default location of domains (only change for NS in center)

% (Navier-Stokes in center of domain)
%file = "data/NS_center";
%file_data = "data/NS_center_data";
%reg_opts.domain_loc = "center";
%file = "data/NS_center_n0.015"; % with noise

% (Navier-Stokes near -1 edge)
file = "data/NS_edge";
file_data = "data/NS_edge_data";
%file = "data/NS_edge_n0.2"; % with noise

% (Divergence-free condition)
%file = "data/div";
%file_data = "data/NS_edge_data";
%file = "data/div_n0.47"; % with noise

% (Boundary conditions)
%file = "data/BC";
%file_data = "data/BC_data";
%file = "data/BC_n50"; % with noise

file_out = file+"_output.mat";
reg_opts.file = file_data;

% Choose mode: BC, NS, or div
mode = "NS";
reg_opts.mode = mode;

% Determine which library terms to search for
glossary = glossaryNS(mode);     % second col contains descriptions
term_names = glossary(:, 1);

% Other inputs:
% Options for library matrix computation
reg_opts.track = 1;              % command window tracking
reg_opts.noise = 0;           % noise level (fraction of total variation)
if mode == "BC"
    N_d = 30;
    N_h = [1, 0, 1, 1];
    reg_opts.env_exp = [6, 0, 6, 6]; % envelope exponents
else
    N_d = 30;
    N_h = [1, 1, 1, 1];
    reg_opts.env_exp = [6, 6, 6, 6]; % envelope exponents
end

% Options for sparse regression
opts = {};

% "pareto": pareto-front method; 
% "threshold": threshold by relative residual
% "multiplicative": threshold by error as fraction of largest term norm
opts.threshold = "pareto";

opts.verbose = 1; % 0 for no cmd window output
% thresholds for "threshold"/"multiplicative" thresholding methods
opts.gamma = 3; opts.epsilon = 0.01;
%% RUN CODE
% Can comment this out if data is loaded

[Q, P, H, char_sizes, valid_single, lib_names] = ...
         SPIDR_NS(term_names, N_d, N_h, reg_opts);
     
%% ----------------------------------------------------------------------
%                               REGRESSION
% -----------------------------------------------------------------------

if mode == "BC"
    % (have both tangential and normal lib_names)
    lib_names_t = lib_names+"_tangential"; 
    lib_names_n = lib_names+"_normal";

    % Threshold regression
    [c_t, lambda_t, best_term_t, lambda1_t] = SparseReg(Q.t, char_sizes, valid_single.t, opts);
    [c_n, lambda_n, best_term_n, lambda1_n] = SparseReg(Q.n, char_sizes, valid_single.n, opts);
    STR_t = lambda1_t/lambda_t; % ratio of single term residual to many-term residual
    best_name_t = lib_names_t(best_term_t);
    best_name_n = lib_names_n(best_term_n);
    STR_n = lambda1_n/lambda_n;
    c_t = c_t/max(abs(c_t));
    c_n = c_n/max(abs(c_n));

    nrmlz_t = max(vecnorm(Q.t));
    nrmlz_n = max(vecnorm(Q.n));
    res = {};
    res.t = norm(c_t'.*Q.t)/nrmlz_t;
    res.n = norm(c_n'.*Q.n)/nrmlz_n;

    % Store results in structs
    coeffs = {};
    nrms = {};

    ind_t = find(c_t);
    c_cell_t = num2cell(c_t'); 
    coeffs.t = cell2struct(c_cell_t(ind_t), lib_names_t(ind_t), 2); % coefficient values

    ind_n = find(c_n);
    c_cell_n = num2cell(c_n'); 
    coeffs.n = cell2struct(c_cell_n(ind_n), lib_names_n(ind_n), 2); % coefficient values

    nrm_vec_t = vecnorm(Q.t.*c_t');
    nrm_cell_t = num2cell(nrm_vec_t/max(abs(nrm_vec_t)));
    nrms.t = cell2struct(nrm_cell_t(ind_t), lib_names_t(ind_t), 2); % norm values

    nrm_vec_n = vecnorm(Q.n.*c_n');
    nrm_cell_n = num2cell(nrm_vec_n/max(abs(nrm_vec_n)));
    nrms.n = cell2struct(nrm_cell_n(ind_n), lib_names_n(ind_n), 2); % norm values
else
    % Threshold regression
    [c, lambda, best_term, lambda1] = SparseReg(Q, char_sizes, valid_single, opts);
    STR = lambda1/lambda;
    best_name = lib_names(best_term);
    c = c/max(abs(c));
    nrmlz = max(vecnorm(Q));
    res = norm(c'.*Q)/nrmlz;

    % Store results in structs
    ind = find(c);
    c_cell = num2cell(c'); 
    coeffs = cell2struct(c_cell(ind), lib_names(ind), 2); % coefficient values

    nrm_vec = vecnorm(Q.*c');
    nrm_cell = num2cell(nrm_vec/max(abs(nrm_vec)));
    nrms = cell2struct(nrm_cell(ind), lib_names(ind), 2); % norm values
end
     
%% SHOW RESULTS

if mode == "BC"
    coeffs.t, coeffs.n        % identified models
    nrms.t, nrms.n            % normalized term strengths
    res.t, res.n              % residual, relative to norm of largest term
    best_name_t, best_name_n  % single best term
    STR_t, STR_n              % ratio of residuals between single and multi-term
else
    coeffs                    % identified models
    nrms                      % normalized term strengths
    res                       % residual, relative to norm of largest term
    best_name                 % single best term
    STR                       % ratio of residuals between single and multi-term
end
% Save output
save(file_out, 'coeffs', 'nrms', 'res', 'Q', 'P', 'H', 'char_sizes', 'valid_single', 'lib_names', '-v7.3');

