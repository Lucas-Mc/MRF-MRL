function r = afun(x,opt,U,plansFor,plansAdj,dimsFor,dimsAdj)
    % x: nR x nC x L x nCh (for)
    % x: nRd x nT x nCh (adj)
    % U: nT x L
    % Jason Otenson 2024
    nT = size(U,1);

    if strcmp(opt,'notransp')
        x = reshape(x,dimsAdj);
        r = forward( x,U,plansFor,dimsAdj,dimsFor );
        r = r(:);
    else
        x = reshape(x,[dimsFor(1) nT dimsFor(3)]);
        r = adjoint( x,U,plansAdj,dimsFor,dimsAdj );
        r = r(:);
    end

end


function y = forward(x,U,plans,dimsIn,dimsOut)
    m1 = dimsIn(1); % nR
    m2 = dimsIn(2); % nC
    m3 = dimsIn(3); % L
    m4 = dimsIn(4); % nCh

    n1 = dimsOut(1); %nRd
    n2 = dimsOut(2); %L
    n3 = dimsOut(3); %nCh
    nT = size(U,1);

    x = reshape(x,[m1 m2 m3*m4] );

    Y = plans.execute( x );
    Y = reshape( Y, [n1 nT n2*n3] );

    y = zeros( n1, nT, n3 );
    for ii = 1:nT
        % p = plans{ii};
        % myY = p.execute(x); % nRd x L x nCh
        myY = Y(:,ii,:);
        myY = reshape( myY,[n1 n2 n3] );
        u = U(ii,:);
        u = repmat( u, [n1 1 n3] );
        myY = sum( myY .* u, 2 );
        y(:,ii,:) = myY;
    end
end

function x = adjoint(y,U,plans,dimsIn,dimsOut)
    m1 = dimsIn(1); %nRd
    m2 = dimsIn(2); %L
    m3 = dimsIn(3); %nCh

    n1 = dimsOut(1); %nR
    n2 = dimsOut(2); %nC
    n3 = dimsOut(3); %L
    n4 = dimsOut(4); %nCh
    nT = size(U,1);

    % x = zeros( n1, n2, n3, n4 );
    Y = zeros( m1, nT, m2*m3 );
    for ii = 1:nT
        u = U(ii,:);
        myY = y(:,ii,:);
        myY = myY .* conj(u);
        Y(:,ii,:,:) = reshape( myY, [m1 1 m2*m3] );
        % myY = reshape( myY, [m1 m2*m3] );
        % p = plans{ii};
        % myX = p.execute(myY);
        % myX = reshape( myX, [n1 n2 n3 n4] );
        % x = x + myX;
    end

    x = plans.execute( Y );
    x = reshape( x, [n1 n2 n3 n4] );

end
