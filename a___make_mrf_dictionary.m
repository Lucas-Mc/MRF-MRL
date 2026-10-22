%% Make T1/T2/B1+ dictionary
% Created by: Jason Otenson 2024
% Modified by: Lucas McCullum 2024-2025
clear all;
close all;
clc;
addpath('./src/')


%% Parameters
% NIST CODATA
gammabar_si = 42.57747852e6;
% Number of dynamics
nExc = 2048;
% TR (s)
tr_s = 20 * 1e-3;
% First echo time (s)
te_s = 3.21 * 1e-3;
% Inversion time (s)
ti_s = 100 * 1e-3;
% Determines FA pattern (other variations in csv file currently ignored)
dict_num = 300;
faMax_deg = 60;
fnMrf = fullfile("data_in", strcat("MRF", num2str(dict_num), ".csv"));

% Fit discretizations
n_bins_T1 = 75;
n_bins_T2 = 75;
n_bins_B1 = 1;
min_T1 = 100;
max_T1 = 2000;
min_T2 = 5;
max_T2 = 500;
min_B1 = 1;
max_B1 = 1;

% No spoiling
phi0_deg = 0;
fnDict = fullfile(...
    "data_out", ...
    strcat(...
        "dictMRF", num2str(dict_num), ...
        "_", num2str(nExc), ...
        "_TR", num2str(tr_s*1e3), ...
        "_TE", replace(num2str(te_s*1e3), ".", "-"), ...
        "_T1-", num2str(min_T1), "-", num2str(max_T1), "-", num2str(n_bins_T1), ...
        "_T2-", num2str(min_T2), "-", num2str(max_T2), "-", num2str(n_bins_T2), ...
        "_B1-", num2str(100*min_B1), "-", num2str(100*max_B1), "-", num2str(n_bins_B1) ...
    ) ...
);
fnDict = strcat(fnDict, ".mat");

% T1 subspace
t1Dict_s = 10.^(linspace(log10(min_T1),log10(max_T1),n_bins_T1)) / 1000;
% T2 subspace
t2Dict_s = 10.^(linspace(log10(min_T2), log10(max_T2), n_bins_T2)) / 1000;
% B1 subspace, make sure this includes 1
b1s = linspace(min_B1, max_B1, n_bins_B1);


%% Derivative params
mrfData = importdata(fnMrf);
b1Scalings = mrfData.data(1:nExc,1);
[T1,T2] = ndgrid(t1Dict_s, t2Dict_s);


%% Dictionary
disp("Generating or loading dictionary...")
tic;

phi0_rad = phi0_deg / 180*pi;
rfSpoilInds = mod(0:(nExc-1), 360);
phi0Max_rad = 360 * (phi0_deg+1) / 180 * pi;
phiRfSpoil_rad = phi0_rad * cumsum(rfSpoilInds);
while any( phiRfSpoil_rad > phi0Max_rad )
    phiRfSpoil_rad( phiRfSpoil_rad > phi0Max_rad ) = phiRfSpoil_rad( phiRfSpoil_rad > phi0Max_rad ) - phi0Max_rad;
end

teSched_s = te_s * ones(nExc,1);
trSched_s = tr_s * ones(nExc,1);
faSched_deg = faMax_deg * b1Scalings;
szOmega = 100;

nT1 = numel(t1Dict_s);
nT2 = numel(t2Dict_s);
nB1 = numel(b1s);
dictTot = zeros(nExc, nT1, nT2, nB1);
for ii = 1:nB1
    fprintf('Generating sub-dictionary for B1 %d of %d...\n',ii,nB1)
    myB1 = b1s(ii);
    S = EPG_MRF_SSFP_par(T1, T2, teSched_s, trSched_s, myB1*faSched_deg, nExc, szOmega, phiRfSpoil_rad*180/pi, ti_s);
    S = S .* repmat(exp(-1i * phiRfSpoil_rad(:)), [1, size(S,2)]);
    dictTot(:,:,:,ii) = reshape(S, [nExc, nT1, nT2]);
end
dictTot = dict_norm(dictTot);

modelDict.t1 = t1Dict_s;
modelDict.t2 = t2Dict_s;
modelDict.T1grid = T1;
modelDict.T2grid = T2;
modelDict.b1 = b1s;
modelDict.te = teSched_s;
modelDict.tr = trSched_s;
modelDict.fa = faSched_deg;
modelDict.szOmega = szOmega;
modelDict.rfSpoil_rad = phiRfSpoil_rad;
modelDict.ti = ti_s;
modelDict.D = dictTot;
save(fnDict, 'modelDict', '-v7.3')

t = toc;
fprintf("Generating or loading dictionary complete in %.2f s.\n",t);
