function Ft = EPG_MRF_SSFP_par( T1s,T2s,TEs,TRs,FAs_deg,nReps,szOmega,phis_deg,TI )
    % returns the complex net transverse magnetization using the extended phase graph approach
    % parallelized to compute a large subset of atoms simultaneously
    % Jason Ostenson 2024

    %   INPUT: T1s/T2s = lists of T1/T2 long./transverse time
    %           constants, the lists should repeat so that [T1s(:)
    %           T2s(:)] list every combination of T1 and T2
    %          TEs = vector of echo times
    %          TRs = vector of repetition times
    %          FAs_deg = vector of flip angles in degrees
    %          nReps = number of excitations
    %          szOmega = pos integer length of the state vectors used to
    %               calculate the magnetization by the EPG formalism (too short
    %               a vector will change the answer since some
    %               magnetization with be lost by a truncated state vector)
    %          phis_deg = vecotr of RF phase angles in degress
    %          TI = inversion time
    %          nT1T2 = number of elements in the subset (number of elements in
    %          T1_list)

    %   OUTPUT: Ft = array containting the complex net trans
    %                   magnetization

    nT1T2 = numel( T1s );
    o1 = 3; o2 = szOmega; o3 = nT1T2;
    o123 = o1*o2*o3; o12 = o1*o2; o23 = o2*o3;

    % lookups
    idx0 = reshape( 1:o123, [o1 o2 o3] );
    idx = idx0;
    idx(1,:,:) = circshift( idx0(1,:,:), 1, 2 );
    idx(2,:,:) = circshift( idx0(2,:,:), -1, 2 );
    idx(3,:,:) = idx0(3,:,:);

    idx_conj_to = 1:o12:o123;
    idx_conj_from = 2:o12:o123;
    idx_zero =  ( o12 - 1 ):o12:o123;

    % convert phi_v to radians
    phis_rad = phis_deg * pi/180;

    % initialize magnetization w/inversion
    omega = single( complex( zeros( o1,o2,o3 ) ) );
    omega(3,1,:) = -1;

    % inversion delay
    omega = relax( omega, TI, T1s, T2s );

    Ft = complex(zeros(nReps, nT1T2, 'single'));
    for ii=1:nReps

        alpha = FAs_deg(ii);
        TR = TRs(ii);
        phi = phis_rad(ii);
        TE = TEs(ii);

        % flip angle transition
        T = rfTransEPG(alpha,phi);
        omega = reshape( T * reshape( omega, [o1 o23] ), [o1 o2 o3] );

        % decay for readout at TE
        omega = relax( omega, TE, T1s, T2s );

        % DEBUG
        % figure(99); clf;
        % plot(log10(abs(omega(1,:))))
        % hold on
        % plot(log10(abs(omega(2,:))))
        % hold off
        % drawnow

        % store readout transverse magnetization
        Ft(ii,:) = omega(1,1,:);

        % apply gradient dephasing specified by delk
        tmp = omega( idx );
        omega = reshape( tmp, [o1 o2 o3] );
        omega( idx_conj_to ) = conj( omega( idx_conj_from ) );
        omega( idx_zero ) = 0;

        % decay for flip angle at TR
        tau = TR-TE;
        omega = relax( omega, tau, T1s, T2s );

    end

end


function w = relax( w,tau,T1s,T2s )
    e2s = exp( -tau./T2s(:).' );
    e1s = exp( -tau./T1s(:).' );
    % w(1,:,:) = squeeze( w(1,:,:) ).* e2s;
    % w(2,:,:) = squeeze( w(2,:,:) ).* e2s;
    % w(3,:,:) = squeeze( w(3,:,:) ).* e1s;
    w(1,:,:) = w(1,:,:) .* shiftdim(e2s,-1);
    w(2,:,:) = w(2,:,:) .* shiftdim(e2s,-1);
    w(3,:,:) = w(3,:,:) .* shiftdim(e1s,-1);
    w(3,1,:) = squeeze( w(3,1,:) ) + (1-e1s(:));
end
