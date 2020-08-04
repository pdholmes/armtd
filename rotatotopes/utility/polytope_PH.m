function [PA, Pb, C] = polytope_PH(Z, options)
% based on the polytope function available in the CORA 2018 toolbox, which
% was written by Matthias Althoff. turns a zonotope representation into a
% polytopic representation.

%obtain number of generators, dimensions
c = Z(:, 1);
G = Z(:, 2:end);

% reduce small generators for numerical reasons:
% compute metric of generators
h = vnorm(G, 1, 2);
% sort generators according to metric
[h_sort, indices] = sort(h, 'descend');
threshold = 0.001;
first_reduce_idx = find(h_sort < threshold, 1, 'first');
if ~isempty(first_reduce_idx)
    Gunred = G(:, indices(1:first_reduce_idx-1));
    Gred = G(:, indices(first_reduce_idx:end));
    % box remaining generators
    d=sum(abs(Gred),2);
    %build box Gbox from interval hull vector d
    Gbox=diag(d);
    G = [Gunred, Gbox];
end
    
[dim,nrOfGenerators]=size(G);

if dim > 2
    if ~exist('options', 'var') || nrOfGenerators > options.maxcombs
        %get number of possible facets
        comb=combinator(nrOfGenerators,dim-1,'c');
    else
        comb = options.combs{nrOfGenerators};
    end

    C=[];
    
    % pairs of generator vectors:
    Q = [G(:, comb(:, 1)); G(:, comb(:, 2))];
    if dim == 3
        % implements a cross product:
        C = [Q(2, :).*Q(6, :) - Q(3, :).*Q(5, :); -(Q(1, :).*Q(6, :) - Q(3, :).*Q(4, :)); Q(1, :).*Q(5, :) - Q(2, :).*Q(4, :)];
    else
        error('Dimension not supported.');
    end
    C = (C./sqrt(sum(C.^2)))';

    %remove NaN rows due to rank deficiency
    index = find(sum(isnan(C),2));
    if ~isempty(index)
        C(index,:) = [];
    end
else
    C = G;
    %... get perpendicular vector
    C = [-C(2, :); C(1, :)];
    C = (C./sqrt(sum(C.^2)))';
end

deltaD = sum(abs((C*G)'))';

d = C*c;

PA = [C; -C];
Pb = [d+deltaD; -d+deltaD];

end