% v2: Limiting to just one city (Delhi), to help simulation time

clear all; 

gps.age      = {'ch','ad','el'};
gps.quar     = {'q0','q1'};
gps.severity = {'mld','svr'};

states1 = {'S','E','A','P','R'};                                           % Stratified by age, quar status
states2 = {'I','Dx1','Dx2'};                                               % Stratified by age, quar status, severity
states3 = {'H','D'};                                                       % Stratified by age
states4 = {'sU','sDx1','sDx2','sQ'};                                       % Not stratified

% states1      = {'S','E','A','P','QS','QE','QA','QP','R','H','D'};
% states2      = {'I','Dx1','Dx2','QI'};

[i, s, d, lim] = get_addresses({states1, gps.quar, gps.age}, [], [], [], 0);
[i, s, d, lim] = get_addresses({states2, gps.quar, gps.age, gps.severity}, i, s, d, lim);
[i, s, d, lim] = get_addresses({states3, gps.age}, i, s, d, lim);
i.nIstates = i.nstates;
[i, s, d, lim] = get_addresses({states4}, i, s, d, lim);
d = char(d);

% This step is to create numerical indices for each of the state variables,
% which makes it easy to set up the governing equations when the number of
% state variables becomes large. This approach depends heavily on Matlab 
% 'structures', which are the same as 'lists' in R. 
% In the above lines:
% - 's' contains the numerical indices of different groups of states. For
% example, s.E gives you all stratifications of E, and s.mld give you all
% states relating to mild illness (see pdf for list of states)
% - 'i' gives the numerical indices of specific states. For example, type
% i.S.el to find that the index of susceptible, elderly people is 3.
% - 'd' is a full list of the state variables, in numerical order

s.infectious = intersect([s.A, s.P, s.I, s.Dx1, s.Dx2], s.q0);             % Indices of all infectious compartments
s.prevalent  = [s.I, s.Dx1, s.Dx2];                                        % Indices of all symptomatic prevalence compartments 
s.S1         = intersect(s.S,s.q0);
s.E1         = intersect(s.E,s.q0);

% --- Include the auxiliaries ---------------------------------------------

numcoms = length(gps.age);
names   = {  'inc',     'hosp',  'pcr', 'rdt'};
lgths   = [numcoms,    numcoms,      2,     2];
for ii = 1:length(names)
    inds = lim + [1:lgths(ii)];
    i.aux.(names{ii}) = inds;
    lim = inds(end);
end
i.nx = lim;

% The above are additional variables that we include in the governing
% equations, to help count 'incidence-like' terms, e.g. the daily
% symptomatic incidence, and the daily admissions to hospital. They're not
% state variables, instead we refer to them as 'auxiliaries', with the
% label 'aux' incorporated in their numerical indices. When we include 
% counts of the numbers of tests being performed per day, they'll be 
% included here as more 'incidence-like' terms. 


% --- Counting symptomatic incidence, by age
mat = zeros(i.nstates);
mat(s.I,:)     = 1;
mat(s.q1,s.q0) = 0; mat(s.q0,s.q1) = 0;
sel.inc = sparse(mat - diag(diag(mat)));

mat = zeros(numcoms,i.nstates);
for ia = 1:length(gps.age)
    mat(ia, intersect(s.I,s.(gps.age{ia}))) = 1;
end
agg.inc = sparse(mat);


% --- Counting hospitalisations
mat = zeros(i.nstates);
mat(s.H,:) = 1;
sel.hosp = sparse(mat - diag(diag(mat)));

mat = zeros(numcoms,i.nstates); row = 1;
for ia = 1:length(gps.age)
    mat(ia, intersect(s.H,s.(gps.age{ia}))) = 1;
end
agg.hosp = sparse(mat);


% --- Counting PCR tests being used amongst those seeking care
mat = zeros(i.nstates);
mat(s.Dx1,:)  = 1;
mat(s.sDx1,:) = 1;
sel.pcr = sparse(mat - diag(diag(mat)));

mat = zeros(2,i.nstates);
mat(1,s.Dx1)  = 1;
mat(2,s.sDx1) = 1;
agg.pcr = sparse(mat);


% --- Counting RDT tests being used to screen general population
vec = zeros(i.nstates,2);
vec(intersect([s.A, s.P, s.I],s.q0),1) = 1;
vec(intersect([s.S, s.E, s.R],s.q0),2) = 1;
sel.rdt = sparse(vec);

% The above lines are setting up matrices to help track the 'incidence-like' 
% terms, for example the incidence of symptomatic illness. Don't worry too 
% much about them for now, will explain


% --- Set up the parameters -----------------------------------------------

R0  = 1;                                                                   % Basic reproduction number
incub_period  = 5;                                                         % Average incubation period (days)
presym_period = 1;                                                         % Average pre-symptomatic period (days)
infec_period  = 5;                                                         % Average infectious period (days)

prop_hosp    = [0.04 0.278 0.599];                                         % Age-specific proportions needing hospitalisation (severe disease)
cfr          = [0.01, 0.3, 6.4]/100;                                       % Age-specific case fatality rates
mort_on_hosp = cfr./prop_hosp;                                             % Amongst those being hospitalised, what proportion would die
dur_in_hosp  = 10;                                                         % Average days in hospital


% --- Set up parameter values
r.incub    = 1/incub_period;                                               % Rate of progression to infectiousness
p.sympto   = 2/3;                                                          % Proportion developing symptomatic infection
p.c        = 2/3;                                                          % Relative infectiousness of asymptomatic vs symptomatic infection
r.eta      = 1/presym_period;                                              % Rate from pre-symptomatic to symptomatic infection
r.gamma    = 1/infec_period;
p.hosp     = prop_hosp;
r.hosp     = 1/5;                                                          % Amongst those with severe disease, average rate of progression to needing hospitalisation
r.mu       = mort_on_hosp*(1/dur_in_hosp);                                 % Rate of mortality amongst those needing hospitalisation
r.gamma_h  = (1-mort_on_hosp)*(1/dur_in_hosp);                             % Rate of recovery amongst those needing hospitalisation

p.imm            = 0;
r.waning         = 0;
prm.PCR_capacity = Inf;
prm.pNCS         = 0.05;
% p.sympto = 0.05;                                                           % What proportion of the population have COVID symptoms, but no COVID

% Interventions
p.sens     = [0.9, 0.8];                                                   % Sensitivity of PCR and RDT testing
p.spec     = [0.98, 0.95];
r.hold     = 2;     %r.hold = 1e3;
r.Dx       = 1;     %r.Dx   = 1e3;
r.screen   = 0;                                                            % Rate of screening asymptomatics
r.careseek = 0;                                                            % Rate at which symptomatics seek care (potential contact with PCR)
r.quar     = 1/7;                                                          % Rate of release from quarantine (used only for non-COVID population)

% In the above lines, notice that we use structures once again, i.e. p and
% r, as a coding convenience. In general, all proportions are listed under
% p and all rates are listed under r. Type 'p' and you'll see the values of
% all proportions involved in the model.


% Get India parameters
prm.N       = [7554531	11954901	818671];
prm.connmat = 1;
prm.mixmat  = [5.59	2.57	0.08
              2.18	5.56	0.08
              0.06	0.12	0.01];       
% 'mixmat' is the contact matrix between different age groups. Note: this is
% only provisional, to be updated in light of recent data
           
r.beta = 1;
r.beta = 1/find_R0(p, r, i, s, gps, prm);
% This step is just to ensure that beta is chosen so that R0 = 1, at this 
% preparatory stage. When we come round to simulating the model, we can 
% then simply multiply beta by the desired value of R0.

% find_R0(p, r, i, s, gps, prm)
% r.careseek = 100; prm.r = r;
% find_R0(p, r, i, s, gps, prm)

prm.r = r; prm.p = p;

save Model_setup;