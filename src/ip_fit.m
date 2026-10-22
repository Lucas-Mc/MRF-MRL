function [pMaps,m0Map,subsMax,cMapsDisc,cMapsDiscInds] = ip_fit(imgData,D,pList,cInds,cMaps)
    % Inner product matching and parameter map generation between imgData and
    % dictionary D using paramter list in pList

    % IN:   imgData = # vox x signal dim
    %       D = signal dim x param1 x param2 ..., columns should be normalized
    %       to  unity
    %       pList = (1 x nP cell) discretized parameter lists, nP = # of
    %       parameters to fit plus number of constraint parameters
    %       cInds = (1 x 2 cell) list of fit inds and constraint inds spec'd as dictionary dimensions
    %       cMaps = (1 x nC cell) constraint maps
    %
    % OUT:  pMaps = (1 x nPc cell) vectorized parameter maps
    %       m0Map = vectorized map of inner product of best-fit dictionary element and signal
    %       subsMax = (1 x nPc cell) index maps of parameters (indices
    %       correspond to ordering in pList cells)
    %       cMapsDic = (1 x nC cell) vectorized and disretized contraint parameter maps
    % Jason Otenson 2024

    chnkSz = 50;
    nV = size( imgData,1 );
    dimsD = size( D );
    nP = numel( pList );
    %
    % if nP ~= (numel(dimsD)-1)
    %     error('Number of parameters in list do not match dictionary dimension.')
    % end

    % dictionary setup
    dInds = cell(1,nP+1);
    for ii = 1:nP+1
        dInds{ii} = 1:size(D,ii);
    end
    nPc = numel(cInds{1});

    if ~isempty(cInds{2})

        % discretize constraint maps
        cMapsDisc = cell( size(cMaps) );
        cMapsDiscInds = cell( size(cMaps) );

        for ii = 1:numel(cInds{2})

            myList = pList{cInds{2}(ii)-1};
            myDist = abs( cMaps{ii}(:) - myList(:).' );
            [~,myInds] = min( myDist,[],2 );
            cMapsDiscInds{ii} = myInds;
            cMapsDisc{ii} = myList(myInds).';

        end

        % for each combination of constraints find matches in constraint maps (i.e.
        % get image indices)
        % form sub-dict
        % check chunk size
        % get ip
        % find max params
        % get param indices of maxima
        % store in param maps at appropriate indices

        indsMap = zeros( nV,1 );
        m0Map = zeros( nV,1 );

        cArray = [cMapsDiscInds{:}];
        cComb = my_comb_fun( dInds(cInds{2}) );
        nComb = size(cComb,1);

        for ii = 1:nComb

            [myMatch,myRows] = ismember( cArray,cComb(ii,:),'rows' );

            if any(myMatch)
                dInds(cInds{2}) = num2cell(cComb(ii,:).');
                subD = D(dInds{:});
                indRows = find(myRows);
                nRows = sum(myRows);
                s = imgData(indRows,:);
                nChnks = ceil( nRows/chnkSz );

                for jj = 1:nChnks

                    myInds = (jj-1)*chnkSz + 1:jj*chnkSz;
                    if jj == nChnks
                        myInds( myInds > nRows ) = [];
                    end
                    ip = s(myInds,:) * conj(subD(:,:));
                    [myM0,myIndsMax] = max(ip,[],2);
                    indsMap(indRows(myInds)) = myIndsMax;
                    m0Map(indRows(myInds)) = myM0;

                end
            end

        end

    else

        % no constraint
        nChnks = ceil(nV/chnkSz);
        indsMap = zeros( nV,1 );
        m0Map = zeros( nV,1 );

        for ii = 1:nChnks
            myInds = (ii-1)*chnkSz + 1:ii*chnkSz;
            if ii == nChnks
                myInds( myInds > nV ) = [];
            end
            ip = imgData(myInds,:) * conj(D(:,:));
            [m0,indsMax] = max(ip,[],2);
            indsMap(myInds,:) = indsMax;
            m0Map(myInds,:) = m0;
        end

        cMapsDisc = [];
        cMapsDiscInds = [];

    end

    subsMax = cell(1,nPc);
    dimsSubD = [dimsD(1) dimsD(cInds{1})];
    [subsMax{:}] = ind2sub(dimsSubD(2:end),indsMap);
    pMaps = cell(1,nPc);

    for ii = 1:nPc
        pMaps{ii} = pList{cInds{1}(ii)-1}(subsMax{ii});
    end

end


function cComb = my_comb_fun( cInds )

    nC = numel(cInds);
    indsTot = zeros(nC,1);

    for ii = 1:nC
        indsTot(ii) = numel(cInds{ii});
    end

    nRow = prod( indsTot(:).' );
    cComb = zeros( nRow,nC );
    myRow = cell( 1,nC );

    for ii = 1:nRow
        [myRow{:}] = ind2sub( indsTot,ii );
        cComb(ii,:) = [myRow{:}];
    end

end
