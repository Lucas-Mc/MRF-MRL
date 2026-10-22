%% Reconstruct images from the raw data
% Created by: Jason Otenson 2024
% Modified by: Lucas McCullum 2024-2025
clear all;
close all;
clc;

diary mylogfile.log

addpath("./src/")
% See: https://finufft.readthedocs.io/en/latest/matlab.html
addpath("./finufft/matlab")
addpath("./contrib/utils")
% See: https://github.com/marcsous/pdff.git
addpath("./contrib/presco/pdff-master/")
% See: https://github.com/marcsous/unwrap.git
addpath("./contrib/presco/unwrap-master/")

% Colormaps
% See: https://www.fabiocrameri.ch/colourmaps/
load('./colormaps/lipari.mat');
load('./colormaps/navia.mat');


%% Settings
disp("Loading definitions...")
tic;

% Whether to save all the figures or not
SAVE_ALL = false;

% Limits, ticks, and units for figures
all_info = struct(...
    'T1', struct(...
        'range', [200,2000],...
        'ticks', [0:500:2000],...
        'units', 'ms'...
    ),...
    'T2', struct(...
        'range', [25,125],...
        'ticks', [0:25:125],...
        'units', 'ms'...
    ),...
    'R2star', struct(...
        'range', [10,100],...
        'ticks', [0:25:100],...
        'units', '1/s'...
    ),...
    'PDFF', struct(...
        'range', [0,100],...
        'ticks', [0:25:100],...
        'units', '%'...
    )...
);

% Scan parameters
object_types = {
    "Calimetrix";
    %"ISMRM-NIST";
    %"FF";
};

for sc = 1:size(object_types,1)
    fprintf("Starting: %d / %d\n",sc,size(object_types,1));
    object_type = object_types{sc,1};
    dict_name = "dictMRF300_2048_TR20_TE3-21_T1-100-2000-75_T2-5-500-75_B1-100-100-1.mat";
    max_T1 = regexp(dict_name, "T1-\d+-(\d+)-\d+", "tokens");
    max_T1 = str2double(max_T1{1});
    max_T2 = regexp(dict_name, "T2-\d+-(\d+)-\d+", "tokens");
    max_T2 = str2double(max_T2{1});
    % Number of dynamics
    nExc = 2048;
    % TR (s)
    tr_s = 20e-3;
    % Slice thickness (mm)
    slice_thickness = 10;
    % Number of echoes
    nTe = 9;
    % Minimum TE (s)
    teMin_s = 3.21e-3;
    % If using M2D
    multislice = false;
    if multislice == true
        % Set to number of slices acquired during scanning
        n_slices = 5;
    end

    % Constants
    % B0 (T)
    b0_T = 1.5;
    % In-vivo
    csH2O_ppm = 4.65;
    csRef_ppm = csH2O_ppm;
    % doi: 10.1002/mrm.24196
    % [terminal methylene, bulk methylene, methylene beta to carbonyl, allylic methylene, methylene alpha to carbonyl, diallylic methylene, glycerol head, glycerol head, water, glycerol head, olefinic methine]
    csTag_ppm =  csRef_ppm - [0.9, 1.3, 1.59, 2.03, 2.25, 2.77, 4.1, 4.3, csH2O_ppm, 5.21, 5.31];
    weightsWAT = [0.0897756,  0.623441, 0, 0.150249/2, 0.150249/2, 0.0205736, 0.0399002/2, 0.0399002/2, 0, 0, 0.0760599];
    % NIST CODATA
    gammabar_si = 42.57747852e6;
    % For PRESCO fat/water separation
    ndb = 2.8;
    % Double check delayed recon DICOM data when changing number of echoes
    deltaTe_s = 1.53e-3; % ms

    if multislice == true
        file_name = strcat(num2str(nExc), "_TR", num2str(tr_s*1000), "_", n_slices, "+", num2str(slice_thickness), "mm_", num2str(nTe), "ec");
    else
        file_name = strcat(num2str(nExc), "_TR", num2str(tr_s*1000), "_", num2str(slice_thickness), "mm_", num2str(nTe), "ec");
    end

    fnDict = fullfile("data_out", dict_name);
    disp("Generating or loading dictionary...")
    tic;
    load(fnDict);
    toc;
    disp("Generating or loading dictionary complete.")
    fnRaw = fullfile("data_in", "scans", object_type, file_name);
    if SAVE_ALL == true
        if ~exist(fullfile("data_out", "scans", object_type), "dir")
           mkdir(fullfile("data_out", "scans", object_type))
        end
    end

    % Sequence: SPGR IR w/variable flip angle
    % Radial angle (rad)
    radialAng_rad = 2/(1+sqrt(5))*pi;
    % Inversion time (s)
    ti_s = 100 * 1e-3;
    % Maximum flip angle (deg)
    faMax_deg = 60;
    % Recovery time b/w dyanamics
    tRest_s = Inf;

    % Image dimension adjustments
    nR = 172;
    nC = 172;

    % Number of iterations used in fitting
    nIter = 40;

    % Low rank factor
    rankFactor = 1 - 1e-4;

    toc;
    disp("Loading definitions complete.")


    %% Derivative parameters
    disp("Generating derivative params...")
    tic;

    te_s = teMin_s:deltaTe_s:((nTe-1) * deltaTe_s + teMin_s);
    nImgSqr = nR * nC;
    freqTag_Hz = gammabar_si * csTag_ppm * 1e-6 * b0_T;
    % doi: 10.1002/mrm.27557
    A0 = exp(1i * 2 * pi * te_s(:) * freqTag_Hz(:)');
    % For FSF / PDFF
    fatModel = A0 * weightsWAT(:);
    fatNorm = norm(fatModel);
    % For UI / MIDBI
    reducedFatInds = [2,6,11];
    % Magnitude, NDB, NMIDB
    % Ctag = [   0    2    0;   % 11 - olefinic
    %            1    0    0;
    %            0    0    0;
    %            2    0    0;
    %            2    0    0;
    %            0    0    2;   % 6  - diallylic
    %            6    0    0;
    %            0    4   -4;
    %            6    0    0;
    %         76.8 -6.5    2;   % 2  - bulk methylene
    %            9    0    0];
    fatModelExpanded = A0(:,reducedFatInds) * diag(weightsWAT(reducedFatInds));
    fatNormExpanded = vecnorm(fatModelExpanded);
    % Avoid divide by 0 error
    fatNormExpanded(fatNormExpanded == 0) = 1;
    % For water
    waterModel = exp(1i * 2 * pi * te_s(:) * (csRef_ppm - csH2O_ppm) * gammabar_si/1e6 * b0_T);
    waterNorm = sqrt(nTe);

    toc;
    disp("Generating derivative params complete.")


    %% Load data
    disp("Loading data...")
    tic;

    % Need to convert Philips LIST/DATA into MATLAB object using your own converter
    load(fullfile("data_in", "scans", object_type, "raw_data.mat"));
    [myRd,myExc,~,myLoc,~,~,myTe,~,~,nCh,~,~] = size(dataRaw, 1:12);
    dataRaw = permute(dataRaw,[[1,2,4,7,10],[3,5,6,8,9,11,12]]);
    nRd = myRd;

    if (myExc ~= nExc) || (myTe ~= nTe)
        error("k-space parameters don't match raw data.")
    end
    % Handle k-space information
    disp("Getting k-space trajectory...")
    tic;
    kSpaceComplex = (-nRd/2:nRd/2-1)'/nRd * exp(1i * (pi/2 + radialAng_rad * (0:nExc-1)));
    kSpace = cat(3, real(kSpaceComplex), imag(kSpaceComplex));
    toc;
    disp("Getting k-space trajectory complete.")


    %% Iterate across all slices
    for tempLoc=1:myLoc
        disp(strcat("Slice = ", num2str(tempLoc), " / ", num2str(myLoc)))
        % [myRd, myExc, myLoc, myTe, nCh]
        temp_dataRaw = dataRaw(:,:,tempLoc,:,:);
        temp_dataRaw = permute(temp_dataRaw,[[1,2,4,5],[3]]);
        temp_dataRaw = permute(temp_dataRaw,[1,3,2,4]);
        % Adjust for bipolar correction?
        temp_dataRaw(:,2:2:end,:,:) = flip(temp_dataRaw(:,2:2:end,:,:),1);
        toc;
        disp("Loading data complete.")

        % Get inter-excitation basis
        disp("Calculating approx. inter-excitation basis ...")
        tic;
        % t long (inter-excitation)
        [u,s,v] = svd(modelDict.D(:,:), 'econ');
        sNrg = cumsum(diag(s).^2) / sum(diag(s).^2);
        % Number of components that explain (100*rankFactor)% of the total energy
        lLong = sum(sNrg <= rankFactor);
        Ulong = u(:, 1:lLong);
        Slong = s(1:lLong, 1:lLong);
        Slong_weighted = Slong / vecnorm(Slong);
        Slong_weighted_reshaped = reshape(Slong_weighted,1,1,[]);
        Vlong = v(:, 1:lLong);
        dCompress = Ulong' * modelDict.D(:,:);
        toc;
        disp("Calculation of inter-excitation basis complete.")

        % Get SENSE maps and recon estimates by TE
        % Y = U * E * C
        % For: U * E
        % Adj: E' * U'

        disp("Calculating coil temp_dataRaw sensitivity maps...")
        tic;

        opt = struct();
        opt.nthread = 0;
        eps = 1e-10;
        plansFor = finufft_plan(2, [nR, nC], -1, lLong*nCh, eps, opt);
        plansFor.setpts(reshape(2*pi*kSpace(:,:,1), [nRd*nExc, 1]), reshape(2*pi*kSpace(:,:,2), [nRd*nExc, 1]));
        plansAdj = finufft_plan(1, [nR, nC], 1, lLong*nCh, eps, opt);
        plansAdj.setpts(reshape(2*pi*kSpace(:,:,1), [nRd*nExc, 1]), reshape(2*pi*kSpace(:,:,2), [nRd*nExc, 1]));
        dimsFor = [nRd, lLong, nCh];
        dimsAdj = [nR, nC, lLong, nCh];
        myCgOp = @(x,dir) afun(x, dir, Ulong, plansFor, plansAdj, dimsFor, dimsAdj);

        % For TV
        usegpu = 0; % set to 1 for gpu acceleration
        TV = TV_operator('2D', usegpu);

        % For LRTV optimization parameters
        param.maxit  = nIter;
        param.step   = 1 / 2^10;    % initial stepsize
        param.lambda = 5e-3;        % 0.5e-3; %2e-3; % TV regularisation weight
        param.backtrack = true;     % backtracking adaptive step size search
        param.dimsAdj = dimsAdj;
        % End TV

        imgStackCompress = zeros(nR, nC, nTe, lLong);
        principal_images_max = 0;

        for iiTe = 1:nTe
            fprintf("Reconstructing image stack for echo %d of %d...\n", iiTe, nTe);

            y = double(temp_dataRaw(:,iiTe,:,:));
            x = solve_LRTV_mod(y(:), myCgOp, TV, param);

            if iiTe == 1

                % End TV
                imgReconMch = reshape(x, dimsAdj);

                disp("Calculating coil sensitivity using matched filtering...")
                % Make SENSE maps using matched filter?
                senseMaps = zeros(nR, nC, nCh);
                for ii = 1:nR
                    for jj = 1:nC
                        % [nR, nC, lLong, nCh]
                        s = imgReconMch(ii,jj,1,:);         % use first singular image
                        x = s(:)*s(:)';                     % rank-1 approximation of coil covariance matrix per-pixel
                        x1 = x(:,1);                        % use first coil as phase reference (coherent phase across coils)
                        senseMaps(ii,jj,:) = x1 / norm(x1); % normalize everything
                    end
                end
                disp("Calculating coil sensitivity using matched filtering complete.")

                figure;
                imagesc(abs(reshape(senseMaps,[nR,nC*nCh])));
                axis image;
                colormap(gray);
                title(strcat("SENSE Maps, Slice ", num2str(tempLoc)))
                if SAVE_ALL == true
                    saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___SENSE-maps___slice-", num2str(tempLoc), ".png")))
                end
                drawnow()

                imgRecon = zeros(nR, nC, lLong);
                for ii = 1:lLong
                    myImg = reshape(imgReconMch(:,:,ii,:),[nR,nC,nCh]) .* conj(senseMaps);
                    imgRecon(:,:,ii) = sum(myImg,3);
                end

            else

                % If multiple coil channels
                if dimsAdj(4) > 1
                    imgRecon = zeros(nR, nC, lLong);
                    x = reshape(x, nR, nC, lLong, []);
                    for ii = 1:lLong
                        myImg = reshape(x(:,:,ii,:),[nR,nC,nCh]) .* conj(senseMaps);
                        imgRecon(:,:,ii) = sum(myImg,3);
                    end
                else
                    imgRecon = reshape(x, dimsAdj(1:end-1));
                end

            end

            imgStackCompress(:,:,iiTe,:) = reshape(imgRecon, [nR,nC,1,lLong]);

            figure;
            clf;
            myImTi = imtile(imgRecon);
            temp_data = abs(myImTi);
            if principal_images_max == 0
                principal_images_max = max(temp_data,[],'all');
            end
            imagesc(temp_data);
            clim([0,principal_images_max]);
            axis image;
            colormap(gray);
            colorbar();
            title(sprintf(strcat("Principal Images, Echo ", num2str(iiTe), ", Slice ", num2str(tempLoc))));
            if SAVE_ALL == true
                saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___principal-images___echo-", num2str(iiTe), "___slice-", num2str(tempLoc), ".png")))
            end
            drawnow()

        end

        toc;
        disp("Calculating coil sensitivity maps complete.")

        % Read in the B0 map
        load(fullfile("data_in", "scans", object_type, "B0_map.mat"));

        % Do fat water separation with PRESCO using the first singular image
        for ii = 1:1
            [fw,~] = presco(te_s, imgStackCompress(:,:,:,ii), b0_T, 'ndb', ndb, 'nonnegFF', 1, 'nonnegR2', 1, 'maxit', [8 10 3]);
        end
        t2StarFitInit = 1 ./ fw.R2;

        % Start: New method for better FSF
        phiFit_rad = fw.PH;
        % For just fat and water
        Areg = conj([waterModel/waterNorm, fatModel/fatNorm]);
        % For UI and MIDBI
        % To include water
        % Areg_expanded = conj([waterModel/waterNorm, fatModelExpanded./fatNormExpanded]);
        % For just fat
        Areg_expanded = conj([fatModelExpanded./fatNormExpanded]);
        waterFits = zeros(nR*nC, lLong);
        fatFits = zeros(nR*nC, lLong);
        bulkFits = zeros(nR*nC, lLong);
        UIFits = zeros(nR*nC, lLong);
        MIDBIFits = zeros(nR*nC, lLong);
        ndbFits = zeros(nR*nC, lLong);
        nmidbFits = zeros(nR*nC, lLong);

        for ii = 1:lLong
            % [nR, nC, nTe, lLong]
            myL = imgStackCompress(:,:,:,ii);
            myL = reshape(myL,[],nTe).';

            for jj = 1:size(myL,2)
                J = diag(exp(1i * (2 * pi * te_s(:) * b0Fit_Hz(jj) + phiFit_rad(jj)) - te_s(:) ./ t2StarFitInit(jj)));
                A = J * Areg;
                c = A \ myL(:,jj);
                % Less noisy fit
                waterFits(jj,ii) = c(1);
                fatFits(jj,ii) = c(2);
                bulkFits(jj,ii) = c(1) + c(2);
                % Noisier, but truer fit? Includes several fat peaks
                A_expanded = J * Areg_expanded;
                c_expanded = A_expanded \ myL(:,jj);
                % waterFits(jj,ii) = c_expanded(1); %c(1);
                % fatFits(jj,ii) = totalFatAmp; %c(2);
                % bulkFits(jj,ii) = c_expanded(1) + sum(c_expanded(2:end)); %c(1) + c(2);
                % Determine UI and MIDBI
                % To include water
                % totalFatAmp = sum(c_expanded(2:end));
                % diallylicAmp = c_expanded(1+2); %c_expanded(1+6);
                % olefinicAmp = c_expanded(1+3); %c_expanded(1+11);
                % For just fat
                totalFatAmp = sum(c_expanded); % c(2)
                diallylicAmp = c_expanded(2); % 6
                olefinicAmp = c_expanded(3); % 11
                temp_UI = abs(real(olefinicAmp) / real(totalFatAmp));
                temp_MIDBI = abs(real(diallylicAmp) / real(totalFatAmp));
                UIFits(jj,ii) = temp_UI;
                MIDBIFits(jj,ii) = temp_MIDBI;
                % Another attempt with closer nmidb but more assumptions
                temp_R1 = abs(real(olefinicAmp) / real(totalFatAmp));
                temp_solution = roots([2*temp_R1*0.093, -6.5*temp_R1-2, 76.8*temp_R1]);
                temp_solution = min(abs(temp_solution(temp_solution > 0)));
                ndbFits(jj,ii) = temp_solution;
                nmidbFits(jj,ii) = 0.093*(temp_solution^2);
            end
        end
        waterFits = reshape(real(waterFits), [nR,nC,lLong]);
        fatFits = reshape(real(fatFits), [nR,nC,lLong]);
        bulkFits = reshape(real(bulkFits), [nR,nC,lLong]);
        % 2.8 ndb in WAT (doi: 10.1016/j.plipres.2008.03.003)
        % 4.59 ndb in soybean oil (doi: 10.1002/mrm.24297)
        UIFits = reshape(real(UIFits), [nR,nC,lLong]);
        % Reconstruct all SVD weights
        UIFits = sum(UIFits.*Slong_weighted_reshaped, 3);
        % Based on diallylic proton contribution
        % 2.07 nmidb in soybean oil (doi: 10.1002/mrm.24297)
        MIDBIFits = reshape(real(MIDBIFits), [nR,nC,lLong]);
        % Reconstruct all SVD weights
        MIDBIFits = sum(MIDBIFits.*Slong_weighted_reshaped, 3);
        % nmidbMap = 3 .* MIDBIFits;
        % ndbMap = 3 .* UIFits + nmidbMap;
        ndbFits = reshape(real(ndbFits), [nR,nC,lLong]);
        % Reconstruct all SVD weights
        ndbFits = sum(ndbFits.*Slong_weighted_reshaped, 3);
        nmidbFits = reshape(real(nmidbFits), [nR,nC,lLong]);
        % Reconstruct all SVD weights
        nmidbFits = sum(nmidbFits.*Slong_weighted_reshaped, 3);
        % doi: 10.1002/mrm.27557
        FufaFits = 100 .* (ndbFits - nmidbFits) ./ 3;
        FufaFits(FufaFits < 0) = NaN;
        % doi: 10.1002/mrm.28471
        FsfaFits = 100 - FufaFits;
        FsfaFits(FsfaFits < 0) = NaN;
        % Assume Ftufa = 0% (otherwise 2%)
        % doi: 10.1194/jlr.D700024-JLR200
        Ftufa = 0;
        FmufaFits = (100 .* ((ndbFits - 2.*nmidbFits) ./ 3)) + Ftufa;
        FmufaFits(FmufaFits < 0) = NaN;
        FpufaFits = (100 .* (nmidbFits ./ 3)) - Ftufa;
        FpufaFits(FpufaFits < 0) = NaN;
        % doi: 10.1016/j.mri.2011.07.004
        % CL = 16.8 + 0.25 * ndb
        CLFits = 16.8 + 0.28 .* ndbFits;
        disp("Fitting for fat-water in T2* basis complete.")

        % Get water and fat mask
        disp("Calculating fat-water masks...")
        tic;

        maskWater = false(nR,nC);
        maskFat = false(nR,nC);
        % New FSF method
        m0WaterInit = reshape(vecnorm(reshape(waterFits, [nR*nC,lLong])'), [nR,nC]);
        m0FatInit = reshape(vecnorm(reshape(fatFits, [nR*nC,lLong])'), [nR,nC]);
        maskWater(m0WaterInit > multithresh(m0WaterInit,1)) = true;
        maskFat(m0FatInit > multithresh(m0FatInit,1)) = true;

        % FSF-based noise mask
        fsfInit = m0FatInit ./ (m0FatInit + m0WaterInit);
        Imtot = abs(m0WaterInit) + abs(m0FatInit);
        noiseMask = false(size(Imtot));
        Mythresh = multithresh(Imtot,20);
        noiseMask(Imtot > Mythresh(1)) = 1;
        noiseMask = imfill(noiseMask, 'holes');

        toc;
        disp("Calculating fat-water masks complete.")

        figure;
        imagesc(1000*1./fw.R2.*noiseMask);
        axis off;
        colormap(navia);
        col = colorbar();
        ylabel(col, 'T2* (ms)');
        clim([0,200]);
        title(strcat("T2* Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T2star__slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(fw.R2.*noiseMask);
        axis off;
        colormap(navia);
        col = colorbar();
        ylabel(col, 'R2* (1/s)');
        clim(all_info.R2star.range);
        title(strcat("R2* Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___R2star-PRESCO__slice-", num2str(tempLoc), ".png")))
        end

        figure;
        PDFF_map = 100*(fw.F./(fw.W + fw.F));
        imagesc(PDFF_map.*noiseMask);
        axis off;
        colormap(turbo);
        col = colorbar();
        ylabel(col, 'PDFF (%)');
        clim(all_info.PDFF.range);
        title(strcat("PDFF Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___PDFF___slice-", num2str(tempLoc), ".png")))
        end

        large_gcf = figure;
        large_gcf.WindowState = 'maximized';
        subplot(131);
        imagesc(m0WaterInit.*noiseMask);
        axis off;
        caxis([0 1]);
        colormap(turbo);
        colorbar();
        title(strcat("Water Fit, Slice ", num2str(tempLoc)));
        subplot(132);
        imagesc(m0FatInit.*noiseMask);
        axis off;
        caxis([0 1]);
        colormap(turbo);
        colorbar();
        title(strcat("Fat Fit, Slice ", num2str(tempLoc)));
        subplot(133);
        imagesc(fsfInit.*noiseMask);
        axis off;
        caxis([0 1]);
        colormap(turbo);
        colorbar();
        title(strcat("FSF, Slice ", num2str(tempLoc)));
        drawnow()
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___FSF___slice-", num2str(tempLoc), ".png")))
        end

        % Fit T1 & B1+
        disp("Estimating T1 and T2 using B1+ map...")
        tic;

        % Constrained fit
        % For constraining by B1+
        pList{1} = modelDict.t1;
        pList{2} = modelDict.t2;
        pList{3} = modelDict.b1;
        cList{1} = [2,3];
        cList{2} = [4];
        cMaps{1} = ones(nR, nC);

        szD = size(modelDict.D);
        Dc = reshape(dCompress, [lLong, szD(2:end)]);
        % New FSF method
        [pMaps_water,m0Map_water,subsMax_water,cMapsDisc_water,cMapsDiscInds_water] = ip_fit(reshape(waterFits,[],lLong),Dc,pList,cList,cMaps);
        [pMaps_fat,m0Map_fat,subsMax_fat,cMapsDisc_fat,cMapsDiscInds_fat] = ip_fit(reshape(fatFits,[],lLong),Dc,pList,cList,cMaps);
        [pMaps_bulk,m0Map_bulk,subsMax_bulk,cMapsDisc_bulk,cMapsDiscInds_bulk] = ip_fit(reshape(bulkFits,[],lLong),Dc,pList,cList,cMaps);

        toc;
        disp("Estimating T1 and T2 using dictionary complete.")

        t1MapWater = 1000 .* reshape(pMaps_water{1}, [nR,nC]);
        t2MapWater = 1000 .* reshape(pMaps_water{2}, [nR,nC]);
        t1MapFat = 1000 .* reshape(pMaps_fat{1}, [nR,nC]);
        t2MapFat = 1000 .* reshape(pMaps_fat{2}, [nR,nC]);
        t1MapBulk = 1000 .* reshape(pMaps_bulk{1}, [nR,nC]);
        t2MapBulk = 1000 .* reshape(pMaps_bulk{2}, [nR,nC]);

        figure;
        imagesc(t1MapWater.*noiseMask);
        axis off;
        colormap(lipari);
        col = colorbar();
        clim(all_info.T1.range);
        ylabel(col, 'T1 (ms)');
        title(strcat("T1 Water Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T1-water___slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(t1MapFat.*noiseMask);
        axis off;
        colormap(lipari);
        col = colorbar();
        clim(all_info.T1.range);
        ylabel(col, 'T1 (ms)');
        title(strcat("T1 Fat Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T1-fat___slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(t1MapBulk.*noiseMask);
        axis off;
        colormap(lipari);
        col = colorbar();
        clim(all_info.T1.range);
        ylabel(col, 'T1 (ms)');
        title(strcat("T1 Bulk Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T1-bulk___slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(t2MapWater.*noiseMask);
        axis off;
        colormap(navia);
        col = colorbar();
        clim(all_info.T2.range);
        ylabel(col, 'T2 (ms)');
        title(strcat("T2 Water Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T2-water___slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(t2MapFat.*noiseMask);
        axis off;
        colormap(navia);
        col = colorbar();
        clim(all_info.T2.range);
        ylabel(col, 'T2 (ms)');
        title(strcat("T2 Fat Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T2-fat___slice-", num2str(tempLoc), ".png")))
        end

        figure;
        imagesc(t2MapBulk.*noiseMask);
        axis off;
        colormap(navia);
        col = colorbar();
        clim(all_info.T2.range);
        ylabel(col, 'T2 (ms)');
        title(strcat("T2 Bulk Fit, Slice ", num2str(tempLoc)));
        if SAVE_ALL == true
            saveas(gcf, fullfile("data_out", "scans", object_type, strcat(file_name, "___T2-bulk___slice-", num2str(tempLoc), ".png")))
        end

        % Save results
        recon_img = abs(imgStackCompress(:,:,1,1));

        results = struct();
        results.b0Map = b0Fit_Hz;
        results.T2starMap = t2StarFitInit;
        results.UIFits = UIFits;
        results.MIDBIFits = MIDBIFits;
        results.ndbFits = ndbFits;
        results.nmidbFits = nmidbFits;
        results.FufaFits = FufaFits;
        results.FsfaFits = FsfaFits;
        results.Ftufa = Ftufa;
        results.FmufaFits = FmufaFits;
        results.FpufaFits = FpufaFits;
        results.CLFits = CLFits;
        results.fsfMap = fsfInit;
        results.imgCoreByTe = imgStackCompress;
        results.Ulong = Ulong;
        results.Slong = Slong;
        results.Slong_weighted = Slong_weighted;
        results.Slong_weighted_reshaped = Slong_weighted_reshaped;
        results.t1MapsWater = t1MapWater;
        results.t2MapsWater = t2MapWater;
        results.t1MapsFat = t1MapFat;
        results.t2MapsFat = t2MapFat;
        results.t1MapsBulk = t1MapBulk;
        results.t2MapsBulk = t2MapBulk;
        results.PDFF_map = PDFF_map;
        results.recon_img = recon_img;
        results.csRef = csRef_ppm;
        results.csTAG = csTag_ppm;
        if exist("Ctag","var")
            results.Ctag = Ctag;
        end
        results.senseMaps = senseMaps;
        results.te_s = te_s;
        if exist("satCoeffs","var")
            results.satCoeffs = satCoeffs;
        end
        results.fnDict = fnDict;
        results.rankFactor = rankFactor;
        if exist("b0s_Hz","var")
            results.b0s_Hz = b0s_Hz;
        end
        if exist("nMfiBasis","var")
            results.nFreqBasis = nMfiBasis;
        end
        if exist("t2Stars_s","var")
            results.t2Stars_s = t2Stars_s;
        end
        % B0, R2, FF, PH, W, F
        if exist("fw","var")
            results.fw = fw;
        end
        results.noiseMask = noiseMask;
        results.maskWater = maskWater;
        results.maskFat = maskFat;

        % Saved file name
        fnSave = fullfile("data_out", "scans", object_type, strcat(file_name, "___recon_test"));

        % All results
        if SAVE_ALL == true
            save(strcat(fnSave, "___slice-", num2str(tempLoc), ".mat"), 'results');
        end

    end

end
