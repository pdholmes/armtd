classdef rotatotope_v2
    %ROTATOTOPE construction and manipulation
    %   Rotatotopes (Definition 4) form the backbone of ARMTD. In short,
    %   they are used to represent a zonotope rotated by sets of rotation
    %   matrices. ARMTD uses rotatotopes to overapproximate the forward
    %   occupany map of a robot link for a set of possible configurations.
    
    properties
        rotation_axes = []; % column array of vectors representing 3D axes of rotation in body-fixed frame
        Jit = {}; % cell of zonotopes representing JRS of i-th joint at time t
        Li; % zonotope representing i-th link (or joint) volume
        
        n_generators = 15; % maximum number of rotatotope generators to keep after reduction
        
        trig_dim = [1, 2]; % cos(q_i) and sin(q_i) dimensions in each Jit
        k_dim = 4; % traj. param dimensions in each Jit, shannon changed from 3

        Vit; % output rotatotope center and generators
        
        dimension; % 3D or 2D space
        
        link_predecessor_joints = {}; %want to keep this info during stacking
        
        % the following are used for bookkeeping: need to remember which 
        % rotatotopes generators were created by multiplication with 
        % k-sliceable generators (k_slc), matrices created from Jit centers, 
        % and the center of Li, as these determine which
        % rotatotope generators are "fully-k^a-sliceable" (eqn. 19) which
        % we keep track of with fully_slc
        
        % keep track of which generators are sliceable by which parameters,
        % and which are fully sliceable:
        k_slc = []; 
        fully_slc = [];
        
        % keep track of intervals of trajectory parameters K
        c_k;
        delta_k;
    end
    
    methods
        function obj = rotatotope_v2(varargin)
            %ROTATOTOPE Construct an instance of this class
            %   obj = rotatotopes(rotation_axes, Jit, Li) where rotation
            %   axes is an array of column vectors specifying 3D axes of
            %   rotation (in body-fixed frame), Jit is a cell of zonotopes
            %   specifying the JRS of predecessor joints for the current
            %   time step, and Li is a zonotope representing link (or
            %   joint) volume to be rotated.
            if nargin == 4
                % parse the input arguments:
                obj.rotation_axes = varargin{1};
                obj.Jit = varargin{2};
                obj.Li = varargin{3}; % see eqn. 1, lemma 7
                obj.link_predecessor_joints = varargin{4};
                if length(obj.Jit) ~= size(obj.rotation_axes, 2)
                    error('Specify as many JRS zonotopes as rotation axes');
                end
                if ~isa(obj.Li, 'zonotope')
                   error('Specify the link volume Li as a zonotope'); 
                end
            else
                error('rotatotope requires 4 arguments');
            end
            
            % infer dimension from link volume
            obj.dimension = size(obj.Li.Z, 1);
            if ~(obj.dimension == 2 || obj.dimension == 3)
                error('Specify the link volume as 2D or 3D');
            end
            
            % initialize outputs
            Vit_tmp = obj.Li.Z;
            if obj.dimension == 2
                Vit_tmp = [Vit_tmp; zeros(1, size(Vit_tmp, 2))];
            end
            n_vec = size(Vit_tmp, 2);
            fully_slc_tmp = zeros(1, n_vec);
            fully_slc_tmp(1) = 1;
            k_slc_tmp = [];
            
            % apply the rotations specified in the JRS:
            for i = length(obj.Jit):-1:1
                
                obj.c_k(i) = obj.Jit{i}.Z(obj.k_dim, 1);

                fully_slc_new = [];
                k_slc_new = [];
                Vit_new = [];

                % multiply link volume by rotation matrix created from
                % center of JRS; no indeterminate is added (see eqn. (13))
                % also see line 9 of Alg. 2
                M = obj.make_matrix(obj.rotation_axes(:, i),...
                    obj.Jit{i}.Z(obj.trig_dim, 1), true);
                Vit_new = M*Vit_tmp;
                fully_slc_new = [fully_slc_new, fully_slc_tmp];
                k_slc_new = [k_slc_new, [zeros(1, n_vec); k_slc_tmp]];
                
                % multiply link volume by rotation matrix created from
                % generators of JRS; indeterminate is added (see eqn. (13)).
                % also see line 9 of Alg. 2
                
                G = obj.Jit{i}.Z(:, 2:end);
                G(:, ~any(G)) = []; % delete zero columns of G
                for j = 1:size(G, 2)
                    M = obj.make_matrix(obj.rotation_axes(:, i),...
                        G(obj.trig_dim, j), false);
                    Vit_new(:, (end+1):(end+size(Vit_tmp, 2))) = M*Vit_tmp;
                    if any(G(obj.k_dim, j)) ~= 0
                        % if generator is k-sliceable, then we can still
                        % evaluate indeterminate later on.
                        obj.delta_k(i) = G(obj.k_dim, j); % save K interval information
                        fully_slc_new = [fully_slc_new, fully_slc_tmp];
                        k_slc_new = [k_slc_new, [ones(1, n_vec); k_slc_tmp]];
                    else
                        fully_slc_new = [fully_slc_new, zeros(1, n_vec)];
                        k_slc_new = [k_slc_new, [zeros(1, n_vec); k_slc_tmp]];
                    end
                end
                
                % reduce number of generators
                % see Appendix D.D 
                [Vit_tmp, fully_slc_tmp, k_slc_tmp] = obj.reduce(Vit_new,...
                    fully_slc_new, k_slc_new);
                n_vec = size(Vit_tmp, 2);
            end 
            
            % check if 2D is desired; project if so
            % SHANNON needed to comment this out because we're keeping the
            % row of zeros in Vit
%             if obj.dimension == 2
%                Vit_tmp(3, :) = []; 
%             end
            
            % store rotatotope
            obj.Vit = Vit_tmp;
                        
            % disregard first column (rotatotope center has no indeterminates)
            obj.fully_slc = fully_slc_tmp(1, 2:end);
            obj.k_slc = k_slc_tmp(:, 2:end);
        end
        
        function M = make_matrix(obj, rotation_axis, x, is_center)
            % given a rotation axis and JRS generator, construct the 3D
            % rotation matrix describing a rotation by cos(q_i) and
            % sin(q_i).
            % this is essentially line 5 of Alg. 2, see also eqn. (12) and
            % Appendix D.C
            
            cq = x(1); % cosine dimension
            sq = x(2); % sine dimension
            if (length(rotation_axis) ~= 3)
                error('Specify a 3D rotation axis. If in 2D, use [0;0;1]');
            end
            
            % use axis-angle formula to get rotation matrix:
            % https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
            e = rotation_axis./norm(rotation_axis);
            K = [0 -e(3) e(2);...
                 e(3) 0 -e(1);...
                 -e(2) e(1) 0];
            if is_center
                M = eye(3) + sq*K + (1 - cq)*K^2;
            else
                M = sq*K + -cq*K^2;
            end
        end
        
        function [Vit_tmp, fully_slc_tmp, k_slc_tmp] = reduce(obj,...
                Vit_new, fully_slc_new, k_slc_new)
            % look at Appendix D.D for more details
            % based off of the "reduceGirard.m" function included in CORA
            
            % only reduce if number of gens is greated than desired
            if size(Vit_new, 2)-1 > obj.n_generators
                c = Vit_new(:, 1);
                G = Vit_new(:, 2:end);
                
                fully_slc_new_G = fully_slc_new(:, 2:end);
                k_slc_new_G = k_slc_new(:, 2:end);

                %compute metric of generators
%               h = vnorm(G,1,1)-vnorm(G,1,inf);
                h = vnorm(G, 1, 2);
                
                % want to add a fix here that prioritizes keeping
                % k-sliceable generators, and re-reducing previously reduced
                % generators (since these are already an axis-aligned box)
                
                % sort generators according to metric
                [~, indices] = sort(h);
                
                n_reduced_generators = length(indices) - obj.n_generators + length(c);

                % pick generators that are reduced
                G_reduced = G(:, indices(1:n_reduced_generators));
                
                % unreduced generators
                G_unreduced = G(:, indices((n_reduced_generators + 1):end));

                % box generators to reduce
                d = sum(abs(G_reduced),2);
                
                % build box Gbox from interval hull vector d
                G_box = diag(d);
                
                % reconstruct new rotatotope
                Vit_tmp = [c, G_unreduced, G_box];
                fully_slc_tmp = [fully_slc_new(1, 1), fully_slc_new_G(1, indices((n_reduced_generators+1):end)), zeros(1, size(G_box, 2))];
                k_slc_tmp = [k_slc_new(:, 1), k_slc_new_G(:, indices((n_reduced_generators+1):end)), zeros(size(k_slc_new, 1), size(G_box, 2))];
            else
                Vit_tmp = Vit_new;
                fully_slc_tmp = fully_slc_new;
                k_slc_tmp = k_slc_new;
            end
        end
        
        function obj = redefine_k_slc( obj, rot_axes_new, Jit_new )
            
            % initialize outputs
            Vit_tmp = obj.Li.Z;
            if obj.dimension == 2
                Vit_tmp = [Vit_tmp; zeros(1, size(Vit_tmp, 2))];
            end
            n_vec = size(Vit_tmp, 2);
            fully_slc_tmp = zeros(1, n_vec);
            fully_slc_tmp(1) = 1;
            k_slc_tmp = [];
            
            % apply the rotations specified in the JRS:
            for i = length(Jit_new):-1:1
                
                % save Kai interval information:
                obj.c_k(i) = Jit_new{i}.Z(obj.k_dim, 1);
                
                fully_slc_new = [];
                k_slc_new = [];
                Vit_new = [];
                
                % multiply link volume by rotation matrix created from
                % center of JRS; no indeterminate is added (see eqn. (13))
                % also see line 9 of Alg. 2
                M = obj.make_matrix(rot_axes_new(:, i),...
                    Jit_new{i}.Z(obj.trig_dim, 1), true);
                Vit_new = M*Vit_tmp;
                fully_slc_new = [fully_slc_new, fully_slc_tmp];
                k_slc_new = [k_slc_new, [zeros(1, n_vec); k_slc_tmp]];
                
                % multiply link volume by rotation matrix created from
                % generators of JRS; indeterminate is added (see eqn. (13)).
                % also see line 9 of Alg. 2
                
                G = Jit_new{i}.Z(:, 2:end);
                G(:, ~any(G)) = []; % delete zero columns of G
                for j = 1:size(G, 2)
                    M = obj.make_matrix(rot_axes_new(:, i),...
                        G(obj.trig_dim, j), false);
                    Vit_new(:, (end+1):(end+size(Vit_tmp, 2))) = M*Vit_tmp;
                        if any(G(obj.k_dim, j)) ~= 0
                            % if generator is k-sliceable, then we can still
                            % evaluate indeterminate later on.
                            obj.delta_k(i) = G(obj.k_dim, j); % save K interval information
                            fully_slc_new = [fully_slc_new, fully_slc_tmp];
                            k_slc_new = [k_slc_new, [ones(1, n_vec);...
                                k_slc_tmp]];
                        else
                            fully_slc_new = [fully_slc_new, zeros(1, n_vec)];
                            k_slc_new = [k_slc_new, [zeros(1, n_vec);...
                                k_slc_tmp]];
                        end
                end
                
                % reduce number of generators
                % see Appendix D.D 
                [Vit_tmp, fully_slc_tmp, k_slc_tmp] = obj.reduce(Vit_new,...
                    fully_slc_new, k_slc_new);
                n_vec = size(Vit_tmp, 2);
            end 
                        
            % disregard first column (rotatotope center has no indeterminates)
            obj.fully_slc = fully_slc_tmp(1, 2:end);
            obj.k_slc = k_slc_tmp(:, 2:end);
            
        end
        
        function obj = stack(obj, prev)
            % "stack" (a.k.a, take Minkowski sum) of this rotatotope and
            % the rotatotope specified by prev (usually, the rotatotope
            % describing possible positions of the predecessor joint)
            % see eqn. (16) for more details, and line 13 of Alg. 2
            
            c = obj.Vit(:, 1) + prev.Vit(:, 1);
            G = [obj.Vit(:, 2:end), prev.Vit(:, 2:end)];
            obj.Vit = [c, G];
            
            obj.fully_slc = [obj.fully_slc, prev.fully_slc];
            
            % if k_slc not same size, add nans
            k_rows_1 = size(obj.k_slc, 1);
            [k_rows_2, k_cols_2] = size(prev.k_slc);

            obj.k_slc = [obj.k_slc, [prev.k_slc; nan(k_rows_1 -...
            k_rows_2, k_cols_2)]];
            
        end
        
        function obj = stack_com(obj, prev)
            % "stack" (a.k.a, take Minkowski sum) of this rotatotope and
            % the rotatotope specified by prev (usually, the rotatotope
            % describing possible positions of the predecessor joint)
            % see eqn. (16) for more details, and line 13 of Alg. 2
            
            c = obj.Vit(:, 1) + prev.Vit(:, 1);
            G = [obj.Vit(:, 2:end), prev.Vit(:, 2:end)];
            obj.Vit = [c, G];
            
            obj.fully_slc = [obj.fully_slc, prev.fully_slc];
            
            % if k_slc not same size, add nans
            k_rows_1 = size(obj.k_slc, 1);
            [k_rows_2, k_cols_2] = size(prev.k_slc);

            obj.k_slc = [obj.k_slc, [prev.k_slc; nan(k_rows_1 -...
                k_rows_2, k_cols_2)]];
            
        end
        
        function Z = slice(obj, k, fully_slc_only)
            % slice all generators, including the not-fully-k-sliceable
            % generators, unless the fully_slc_only flag is true.
            % this algorithm is described in Alg. 1
            if ~exist('fully_slc_only', 'var')
                fully_slc_only = false;
            end
            if length(k) ~= length(obj.c_k)
                error('Slice point not correct dimension');
            end
            c = obj.Vit(:, 1);
            G_sliced = obj.Vit(:, 2:end);
            slice_coefficient = zeros(length(k), 1);
            for i = 1:length(k)
                if abs(k(i) - obj.c_k(i)) > obj.delta_k(i)
                    error('Slice point is out of bounds');
                end
                slice_coefficient(i) = (k(i) - obj.c_k(i))/obj.delta_k(i);
                % see Alg. 1, line 6:
                if ~fully_slc_only
                    G_sliced(:, obj.k_slc(i, :) == 1) =...
                        G_sliced(:, obj.k_slc(i, :) == 1)*slice_coefficient(i); % slice gens
                else
                    G_sliced(:, (obj.k_slc(i, :) == 1 &...
                        obj.fully_slc == 1)) = G_sliced(:,...
                        (obj.k_slc(i, :) == 1 & obj.fully_slc == 1))*slice_coefficient(i); % slice gens
                end
            end
            
            % take the fully sliced generators, add to center
            % see Alg. 1, line 11
            c_out = c + sum(G_sliced(:, obj.fully_slc == 1), 2);
            G_out = G_sliced(:, ~(obj.fully_slc == 1));
            
            Z = [c_out, G_out];
        end
        
        function [p] = plot(obj, color, buffer, p)
            if ~exist('color', 'var')
                color = 'b';
            end
            if ~exist('buffer', 'var')
                buffer = 0;
            end
            if ~exist('p', 'var')
                p = [];
            end
            %Shannon added: back down to 2!!
            if obj.dimension == 2
                Z = zonotope([obj.Vit(1:2, :), buffer/2*eye(obj.dimension)]);
            else
                Z = zonotope([obj.Vit, buffer/2*eye(obj.dimension)]);
            end
            switch obj.dimension
                case 2
                    p = plotFilled(Z, [1, 2], color);
                    p.FaceAlpha = 0.1;
                case 3
                    if isempty(p)
                        Z = reduce(Z, 'girard', 5);
                        V = vertices(project(Z, [1, 2, 3]));
                        shp = alphaShape(V(1, :)', V(2, :)', V(3, :)', inf);
                        p = plot(shp);
                        p.FaceAlpha = 0.05;
                        p.EdgeAlpha = 0.15;
                        p.EdgeColor = 'k';
                        p.FaceColor = color;
                    else
                        Z = reduce(Z, 'girard', 5);
                        V = vertices(project(Z, [1, 2, 3]));
                        shp = alphaShape(V(1, :)', V(2, :)', V(3, :)', inf);
                        [tri, P] = shp.alphaTriangulation();
                        p.Faces = tri;
                        p.Vertices = P;
                    end
            end
        end
        
        function [p] = plot_slice(obj, k, color, buffer, p)
            if ~exist('color', 'var')
                color = 'g';
            end
            if ~exist('buffer', 'var')
                buffer = 0;
            end
            if ~exist('p', 'var')
                p = [];
            end
            
            % use this flag to specify slicing all generators, or only the
            % fully sliceable generators:
            fully_slc_only = true;
            
            Z = obj.slice(k, fully_slc_only);
            %Shannon added: back down to 2!!
            if obj.dimension == 2
                Z = Z(1:2, :);
            end
            Z = zonotope([Z, buffer/2*eye(obj.dimension)]);
            switch obj.dimension
                case 2
                    p = plotFilled(Z, [1, 2], color);
                    p.FaceAlpha = 0.1;
                case 3
                    if isempty(p)
                        Z = reduce(Z, 'girard', 5);
                        V = vertices(project(Z, [1, 2, 3]));
                        shp = alphaShape(V(1, :)', V(2, :)', V(3, :)', inf);
                        p = plot(shp);
                        p.FaceAlpha = 0.05;
                        p.EdgeAlpha = 0.15;
                        p.EdgeColor = 'k';
                        p.FaceColor = color;
                    else
                        Z = reduce(Z, 'girard', 5);
                        V = vertices(project(Z, [1, 2, 3]));
                        shp = alphaShape(V(1, :)', V(2, :)', V(3, :)', inf);
                        [tri, P] = shp.alphaTriangulation();
                        p.Faces = tri;
                        p.Vertices = P;
                    end
            end
        end
    end
end

