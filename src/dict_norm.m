function [Dnorm,v] = dict_norm(D)
    szD = size(D);
    v = vecnorm(D(:,:));
    Dnorm = D(:,:) ./ repmat( v, [szD(1) 1] );
    Dnorm = reshape( Dnorm, szD );
end
