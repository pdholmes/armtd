classdef rotatotope_v3
    %ROTATOTOPE construction and manipulation
    %   Rotatotopes (Definition 4) form the backbone of ARMTD. In short,
    %   they are used to represent a zonotope rotated by sets of rotation
    %   matrices. ARMTD uses rotatotopes to overapproximate the forward
    %   occupany map of a robot link for a set of possible configurations.
    
    properties
        rotation_axes = []; % column array of vectors representing 3D axes of rotation in body-fixed frame
        Jit = {}; % cell of zonotopes representing JRS of i-th joint at time t
        Li; % zonotope representing i-th link (or joint) volume
        
        n_generators = 30; % maximum number of rotatotope generators to keep after reduction
        
        trig_dim = [1, 2]; % cos(q_i) and sin(q_i) dimensions in each Jit
        
        k_dim = {}; % for each JRS, dimensions of trajectory parameters
        k_names = {}; % names of trajectory parameters for each JRS
        k_list = {}; % will store list of all trajectory parameters
        
        Vit; % output rotatotope center and generators
        
        dimension; % 3D or 2D space
        
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

        % used for tracking gen coefficients:
        gen_tracker = [];
        gen_tracker_list = {};
        gen_tracker_JRS_idxs = {};
        gen_tracker_reduction_idxs = {};
        n_total_gen;
    end
    
    methods
        function obj = rotatotope_v3(varargin)
            %ROTATOTOPE Construct an instance of this class
            %   obj = rotatotopes(rotation_axes, Jit, Li) where rotation
            %   axes is an array of column vectors specifying 3D axes of
            %   rotation (in body-fixed frame), Jit is a cell of zonotopes
            %   specifying the JRS of predecessor joints for the current
            %   time step, and Li is a zonotope representing link (or
            %   joint) volume to be rotated.
            if nargin == 8
                % parse the input arguments:
                obj.rotation_axes = varargin{1};
                obj.Jit = varargin{2};
                obj.Li = varargin{3}; % see eqn. 1, lemma 7
                obj.k_dim = varargin{4};
                obj.k_names = varargin{5};
                obj.gen_tracker_list = varargin{6};
                obj.gen_tracker_JRS_idxs = varargin{7};
                obj.gen_tracker_reduction_idxs = varargin{8};
                obj.n_total_gen = length(obj.gen_tracker_list);
                if length(obj.Jit) ~= size(obj.rotation_axes, 2)
                    error('Specify as many JRS zonotopes as rotation axes');
                end
                if ~isa(obj.Li, 'zonotope')
                   error('Specify the link volume Li as a zonotope'); 
                end

                obj = obj.multiply();
            else
                % error('rotatotope requires 5 arguments');
            end

        end

        function obj = multiply(obj)
            
            % infer dimension from link volume
            obj.dimension = size(obj.Li.Z, 1);
            if ~(obj.dimension == 2 || obj.dimension == 3)
                error('Specify the link volume as 2D or 3D');
            end
            
            % initialize K information
            for i = 1:length(obj.Jit)
                % get generator matrix
                G = obj.Jit{i}.Z(:, 2:end);
                G(:, ~any(G)) = []; % delete zero columns of G
                
                n_k = length(obj.k_dim{i});
                for j = 1:n_k
                    if ~any(strcmp(obj.k_names{i}{j}, obj.k_list))
                        % add parameter to list
                        obj.k_list{end+1, 1} = obj.k_names{i}{j};
                        obj.c_k(end+1, 1) = obj.Jit{i}.Z(obj.k_dim{i}(j), 1);
                        k_col = find(G(obj.k_dim{i}(j), :) ~= 0);
                        if length(k_col) == 1
                            obj.delta_k(end+1, 1) = G(obj.k_dim{i}(j), k_col);
                        elseif length(k_col) > 1
                            error('More than one generator for k-dimension');
                        elseif isempty(k_col)
                            error('No generator found for that k-dimension');
                        end
                    else
                        % we've already seen this parameter... check
                        % that the intervals for parameter are the same.
                        k_idx = find(strcmp(obj.k_names{i}{j}, obj.k_list));
                        c_k_tmp = obj.Jit{i}.Z(obj.k_dim{i}(j), 1);
                        k_col = find(G(obj.k_dim{i}(j), :) ~= 0);
                        if length(k_col) == 1
                            delta_k_tmp = G(obj.k_dim{i}(j), k_col);
                        elseif length(k_col) > 1
                            error('More than one generator for k-dimension');
                        elseif isempty(k_col)
                            error('No generator found for that k-dimension');
                        end
                        if (c_k_tmp ~= obj.c_k(k_idx, 1)) || (delta_k_tmp ~= obj.delta_k(k_idx, 1))
                            error('Same parameter used for multiple joints, but parameter range is different');
                        end
                    end
                end
            end
            
            % initialize outputs
            Vit_tmp = obj.Li.Z;
            if obj.dimension == 2
                Vit_tmp = [Vit_tmp; zeros(1, size(Vit_tmp, 2))];
            end
            n_k = length(obj.k_list);
            n_vec = size(Vit_tmp, 2);
            fully_slc_tmp = zeros(1, n_vec);
            fully_slc_tmp(1) = 1;
            k_slc_tmp = zeros(n_k, n_vec);
            gen_tracker_tmp = zeros(obj.n_total_gen, n_vec);
            
            % apply the rotations specified in the JRS:
            for i = length(obj.Jit):-1:1
                
                % get generator matrix
                G = obj.Jit{i}.Z(:, 2:end);
%                 G(:, ~any(G)) = []; % delete zero columns of G
                
                fully_slc_new = [];
                k_slc_new = [];
                Vit_new = [];
                gen_tracker_new = [];
                
                % multiply link volume by rotation matrix created from
                % center of JRS; no indeterminate is added (see eqn. (13))
                % also see line 9 of Alg. 2
                M = obj.make_matrix(obj.rotation_axes(:, i),...
                    obj.Jit{i}.Z(obj.trig_dim, 1), true);
                Vit_new = M*Vit_tmp;
                fully_slc_new = [fully_slc_new, fully_slc_tmp];
                k_slc_new = [k_slc_new, k_slc_tmp];
                gen_tracker_new = [gen_tracker_new, gen_tracker_tmp];
                
                % multiply link volume by rotation matrix created from
                % generators of JRS; indeterminate is added (see eqn. (13)).
                % also see line 9 of Alg. 2
                if size(G, 2) ~= length(obj.gen_tracker_JRS_idxs{i})
                    error('number of JRS generators doesn''t match gen tracker');
                end
                for j = 1:size(G, 2)
                    M = obj.make_matrix(obj.rotation_axes(:, i),...
                        G(obj.trig_dim, j), false);
                    Vit_new(:, (end+1):(end+size(Vit_tmp, 2))) = M*Vit_tmp;
                    if any(G(obj.k_dim{i}, j)) ~= 0
                        % if generator is k-sliceable, then we can still
                        % evaluate indeterminate later on.
                        % add 1 to this row of k_slc, implying an
                        % additional indeterminate
                        k_row = find(G(obj.k_dim{i}, j) ~= 0);
                        k_name = obj.k_names{i}{k_row};
                        k_idx = find(strcmp(k_name, obj.k_list));
                        
                        fully_slc_new = [fully_slc_new, fully_slc_tmp];
                        k_slc_tmp_2 = k_slc_tmp;
                        k_slc_tmp_2(k_idx, :) = k_slc_tmp_2(k_idx, :) + 1;
                        k_slc_new = [k_slc_new, k_slc_tmp_2];
                    else
                        fully_slc_new = [fully_slc_new, zeros(1, n_vec)];
                        k_slc_new = [k_slc_new, k_slc_tmp];
                    end
                    gen_tracker_tmp_2 = gen_tracker_tmp;
                    gen_tracker_tmp_2(obj.gen_tracker_JRS_idxs{i}(j), :) = gen_tracker_tmp_2(obj.gen_tracker_JRS_idxs{i}(j), :) + 1;
                    gen_tracker_new = [gen_tracker_new, gen_tracker_tmp_2];
                end
                
                % reduce number of generators
                % see Appendix D.D 
                [Vit_tmp, fully_slc_tmp, k_slc_tmp, gen_tracker_tmp] = obj.reduce(Vit_new, fully_slc_new, k_slc_new, gen_tracker_new, obj.gen_tracker_reduction_idxs{i});
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
            obj.gen_tracker = gen_tracker_tmp(:, 2:end);

            % combine generators with same coefficients:
            obj = obj.combine_generators();
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
        
        function [Vit_tmp, fully_slc_tmp, k_slc_tmp, gen_tracker_tmp] = reduce(obj, Vit_new, fully_slc_new, k_slc_new, gen_tracker_new, gen_tracker_reduction_idxs)
            % look at Appendix D.D for more details
            % based off of the "reduceGirard.m" function included in CORA
            
            % only reduce if number of gens is greated than desired
            if size(Vit_new, 2)-1 > obj.n_generators
                c = Vit_new(:, 1);
                G = Vit_new(:, 2:end);
                
                fully_slc_new_G = fully_slc_new(:, 2:end);
                k_slc_new_G = k_slc_new(:, 2:end);
                gen_tracker_new_G = gen_tracker_new(:, 2:end);

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
                gen_tracker_reduction_tmp = zeros(size(gen_tracker_new, 1), size(G_box, 2));
                gen_tracker_reduction_tmp(gen_tracker_reduction_idxs, :) = eye(3);
                gen_tracker_tmp = [gen_tracker_new(:, 1), gen_tracker_new_G(:, indices((n_reduced_generators+1):end)), gen_tracker_reduction_tmp];
            else
                Vit_tmp = Vit_new;
                fully_slc_tmp = fully_slc_new;
                k_slc_tmp = k_slc_new;
                gen_tracker_tmp = gen_tracker_new;
            end
        end
        
        % function obj = redefine_k_slc( obj, rot_axes_new, Jit_new )

        %      % initialize K information
        %     for i = 1:length(Jit_new)
        %         % get generator matrix
        %         G = Jit_new{i}.Z(:, 2:end);
        %         G(:, ~any(G)) = []; % delete zero columns of G

        %         n_k = length(obj.k_dim{i});
        %         for j = 1:n_k
        %             if ~any(strcmp(obj.k_names{i}{j}, obj.k_list))
        %                 % add parameter to list
        %                 obj.k_list{end+1, 1} = obj.k_names{i}{j};
        %                 obj.c_k(end+1, 1) = Jit_new{i}.Z(obj.k_dim{i}(j), 1);
        %                 k_col = find(G(obj.k_dim{i}(j), :) ~= 0);
        %                 if length(k_col) == 1
        %                     obj.delta_k(end+1, 1) = G(obj.k_dim{i}(j), k_col);
        %                 elseif length(k_col) > 1
        %                     error('More than one generator for k-dimension');
        %                 elseif isempty(k_col)
        %                     error('No generator found for that k-dimension');
        %                 end
        %             else
        %                 % we've already seen this parameter... check
        %                 % that the intervals for parameter are the same.
        %                 k_idx = find(strcmp(obj.k_names{i}{j}, obj.k_list));
        %                 c_k_tmp = Jit_new{i}.Z(obj.k_dim{i}(j), 1);
        %                 k_col = find(G(obj.k_dim{i}(j), :) ~= 0);
        %                 if length(k_col) == 1
        %                     delta_k_tmp = G(obj.k_dim{i}(j), k_col);
        %                 elseif length(k_col) > 1
        %                     error('More than one generator for k-dimension');
        %                 elseif isempty(k_col)
        %                     error('No generator found for that k-dimension');
        %                 end
        %                 if (c_k_tmp ~= obj.c_k(k_idx, 1)) ||...
        %                         (delta_k_tmp ~= obj.delta_k(k_idx, 1))
        %                     error('Same parameter used for multiple joints, but parameter range is different');
        %                 end
        %             end
        %         end
        %     end

        %     % initialize outputs
        %     Vit_tmp = obj.Li.Z;
        %     if obj.dimension == 2
        %         Vit_tmp = [Vit_tmp; zeros(1, size(Vit_tmp, 2))];
        %     end
        %     n_k = length(obj.k_list);
        %     n_vec = size(Vit_tmp, 2);
        %     fully_slc_tmp = zeros(1, n_vec);
        %     fully_slc_tmp(1) = 1;
        %     k_slc_tmp = zeros(n_k, n_vec);

        %     % apply the rotations specified in the JRS:
        %     for i = length(Jit_new):-1:1

        %         % get generator matrix
        %         G = Jit_new{i}.Z(:, 2:end);
        %         G(:, ~any(G)) = []; % delete zero columns of G

        %         fully_slc_new = [];
        %         k_slc_new = [];
        %         Vit_new = [];

        %         % multiply link volume by rotation matrix created from
        %         % center of JRS; no indeterminate is added (see eqn. (13))
        %         % also see line 9 of Alg. 2
        %         M = obj.make_matrix(rot_axes_new(:, i),...
        %             Jit_new{i}.Z(obj.trig_dim, 1), true);
        %         Vit_new = M*Vit_tmp;
        %         fully_slc_new = [fully_slc_new, fully_slc_tmp];
        %         k_slc_new = [k_slc_new, k_slc_tmp];

        %         % multiply link volume by rotation matrix created from
        %         % generators of JRS; indeterminate is added (see eqn. (13)).
        %         % also see line 9 of Alg. 2
        %         for j = 1:size(G, 2)
        %             M = obj.make_matrix(rot_axes_new(:, i),...
        %                 G(obj.trig_dim, j), false);
        %             Vit_new(:, (end+1):(end+size(Vit_tmp, 2))) = M*Vit_tmp;
        %             if any(G(obj.k_dim{i}, j)) ~= 0
        %                 % if generator is k-sliceable, then we can still
        %                 % evaluate indeterminate later on.
        %                 % add 1 to this row of k_slc, implying an
        %                 % additional indeterminate
        %                 k_row = find(G(obj.k_dim{i}, j) ~= 0);
        %                 k_name = obj.k_names{i}{k_row};
        %                 k_idx = find(strcmp(k_name, obj.k_list));

        %                 fully_slc_new = [fully_slc_new, fully_slc_tmp];
        %                 k_slc_tmp_2 = k_slc_tmp;
        %                 k_slc_tmp_2(k_idx, :) = k_slc_tmp_2(k_idx, :) + 1;
        %                 k_slc_new = [k_slc_new, k_slc_tmp_2];
        %             else
        %                 fully_slc_new = [fully_slc_new, zeros(1, n_vec)];
        %                 k_slc_new = [k_slc_new, k_slc_tmp];
        %             end
        %         end

        %         % reduce number of generators
        %         % see Appendix D.D 
        %         [Vit_tmp, fully_slc_tmp, k_slc_tmp] = obj.reduce(Vit_new,...
        %             fully_slc_new, k_slc_new);
        %         n_vec = size(Vit_tmp, 2);
        %     end 

        %     % disregard first column (rotatotope center has no indeterminates)
        %     obj.fully_slc = fully_slc_tmp(1, 2:end);
        %     obj.k_slc = k_slc_tmp(:, 2:end);

        % end
        
        function obj = stack(obj, prev, stack_coeff)
            % "stack" (a.k.a, take Minkowski sum) of this rotatotope and
            % the rotatotope specified by prev (usually, the rotatotope
            % describing possible positions of the predecessor joint)
            % see eqn. (16) for more details, and line 13 of Alg. 2
            %
            % new argument `stack_coeff` let's you multiply previous
            % rotatotope by a coefficient when stacking.

            if ~exist('stack_coeff', 'var')
                stack_coeff = 1;
            end
            
            c = obj.Vit(:, 1) + stack_coeff*prev.Vit(:, 1);
            G = [obj.Vit(:, 2:end), stack_coeff*prev.Vit(:, 2:end)];
            obj.Vit = [c, G];
            
            obj.fully_slc = [obj.fully_slc, prev.fully_slc];
            obj.gen_tracker = [obj.gen_tracker, prev.gen_tracker];

            % merge lists of parameters and indeterminates
            new_k_list = unique([obj.k_list(:); prev.k_list(:)]);
            new_k_slc = [];
            c_k = [];
            delta_k = [];
            for i = 1:length(new_k_list)
                curr_k_idx = find(strcmp(new_k_list{i}, obj.k_list));
                prev_k_idx = find(strcmp(new_k_list{i}, prev.k_list));
                if isempty(curr_k_idx)
                    k_slc_1 = zeros(1, size(obj.Vit(:, 2:end), 2));
                    k_slc_2 = prev.k_slc(prev_k_idx, :);
                    c_k(i, 1) = prev.c_k(prev_k_idx);
                    delta_k(i, 1) = prev.delta_k(prev_k_idx);
                elseif isempty(prev_k_idx)
                    k_slc_1 = obj.k_slc(curr_k_idx, :);
                    k_slc_2 = zeros(1, size(prev.Vit(:, 2:end), 2));
                    c_k(i, 1) = obj.c_k(curr_k_idx);
                    delta_k(i, 1) = obj.delta_k(curr_k_idx);
                else
                    k_slc_1 = obj.k_slc(curr_k_idx, :);
                    k_slc_2 = prev.k_slc(prev_k_idx, :);
                    if (obj.c_k(curr_k_idx) ~= prev.c_k(prev_k_idx)) || (obj.delta_k(curr_k_idx) ~= prev.delta_k(prev_k_idx))
                        error('shared parameter detected, but ranges are different.')
                    end
                    c_k(i, 1) = obj.c_k(curr_k_idx);
                    delta_k(i, 1) = obj.delta_k(curr_k_idx);
                end
                new_k_slc(i, :) = [k_slc_1, k_slc_2];
            end
            
            obj.k_list = new_k_list;
            obj.k_slc = new_k_slc;
            obj.c_k = c_k;
            obj.delta_k = delta_k;

            % combine generators with same coefficients:
            obj = obj.combine_generators();
        end
        
        function Z = slice(obj, k_names, k_values, fully_slc_only)
            % slice all generators, including the not-fully-k-sliceable
            % generators, unless the fully_slc_only flag is true.
            % this algorithm is described in Alg. 1
            if ~exist('fully_slc_only', 'var')
                fully_slc_only = false;
            end
            c = obj.Vit(:, 1);
            G_sliced = obj.Vit(:, 2:end);
            slice_coefficient = zeros(length(k_names), 1);
            k_slc_tmp = obj.k_slc; % use this too keep track of gens that have been sliced
            for i = 1:length(k_names)
                k_idx = find(strcmp(k_names{i}, obj.k_list));
                if isempty(k_idx)
                    disp([k_names{i} ' not found in k_list.']);
                    continue;
                end
                if abs(k_values(i) - obj.c_k(k_idx)) > obj.delta_k(k_idx)
                    error('Slice point is out of bounds');
                end
                slice_coefficient(i) = (k_values(i) - obj.c_k(k_idx))/obj.delta_k(k_idx);
                % see Alg. 1, line 6:
                if ~fully_slc_only
                    slice_gen_idxs = find(obj.k_slc(k_idx, :) ~= 0);
                else
                    slice_gen_idxs = find(obj.k_slc(k_idx, :) ~= 0 & obj.fully_slc == 1);
                end
                for j = 1:length(slice_gen_idxs)
                    % slice gens... if indeterminate appears twice, then
                    % square coefficient, three times => cube, etc.
                    G_sliced(:, slice_gen_idxs(j)) = G_sliced(:, slice_gen_idxs(j))*slice_coefficient(i)^obj.k_slc(k_idx, slice_gen_idxs(j));
                    k_slc_tmp(k_idx, slice_gen_idxs(j)) = 0; % remove indeterminate coefficients associated with parameter.
                end
            end
            
            % take the fully sliced generators, add to center
            % see Alg. 1, line 11
            c_out = c + sum(G_sliced(:, obj.fully_slc == 1 & all(k_slc_tmp == 0, 1)), 2);
            G_out = G_sliced(:, ~(obj.fully_slc == 1 & all(k_slc_tmp == 0, 1)));
            
            Z = [c_out, G_out];
        end

        function obj = post_slice(obj, k_names, k_values)
            % post_slice... basically like slice, but updates the rotatotope.
            % need to implement so that k_values can be a range.
            % also assuming that the gens don't have to be fully sliceable.

            if size(k_values, 2) == 2
                for i = 1:size(k_values, 1)
                    if k_values(i, 2) < k_values(i, 1)
                        error('If specifying a range, make sure second element is larger than first.');
                    end
                end
            end

            c_sliced = obj.Vit(:, 1);
            G_sliced = obj.Vit(:, 2:end);
            slice_coefficient = zeros(length(k_names), 1);
            k_slc_tmp = obj.k_slc; % use this to keep track of gens that have been sliced
            fully_slc_tmp = obj.fully_slc; % use this to remove full sliceability (lol) for gens that are sliced, but not completely.
            gen_tracker_tmp = obj.gen_tracker;
            for i = 1:length(k_names)
                k_idx = find(strcmp(k_names{i}, obj.k_list));
                if isempty(k_idx)
                    disp([k_names{i} ' not found in k_list.']);
                    continue;
                end
                for j = 1:size(k_values(i, :), 2)
                    if abs(k_values(i, j) - obj.c_k(k_idx)) > obj.delta_k(k_idx)
                        error('Slice point is out of bounds');
                    end
                    slice_coefficient(i, j) = (k_values(i, j) - obj.c_k(k_idx))/obj.delta_k(k_idx);
                end
                % see Alg. 1, line 6:
                slice_gen_idxs = find(obj.k_slc(k_idx, :) ~= 0);

                for j = 1:length(slice_gen_idxs)
                    % slice gens... if indeterminate appears twice, then
                    % square coefficient, three times => cube, etc.
                    if size(slice_coefficient(i, :), 2) == 1 || (size(slice_coefficient(i, :), 2) == 2 && (slice_coefficient(i, 1) == slice_coefficient(i, 2)))
                        G_sliced(:, slice_gen_idxs(j)) = G_sliced(:, slice_gen_idxs(j))*slice_coefficient(i)^obj.k_slc(k_idx, slice_gen_idxs(j));
                        k_slc_tmp(k_idx, slice_gen_idxs(j)) = 0; % remove indeterminate coefficients associated with parameter.
                    elseif size(slice_coefficient(i, :), 2) == 2 % slicing by range
                        G_sliced_lb = G_sliced(:, slice_gen_idxs(j))*slice_coefficient(i, 1)^obj.k_slc(k_idx, slice_gen_idxs(j));
                        G_sliced_ub = G_sliced(:, slice_gen_idxs(j))*slice_coefficient(i, 2)^obj.k_slc(k_idx, slice_gen_idxs(j));
                        if ~fully_slc_tmp(slice_gen_idxs(j))
                            G_sliced(:, slice_gen_idxs(j)) = max(G_sliced_lb, G_sliced_ub);
                            k_slc_tmp(k_idx, slice_gen_idxs(j)) = 0; % remove indeterminate coefficients associated with parameter.
                        else
                            % woo. this gen is fully sliceable, so we know the value's between G_sliced_lb and G_sliced_ub
                            k_slc_tmp_tmp = k_slc_tmp(:, slice_gen_idxs(j));
                            k_slc_tmp_tmp(k_idx) = [];
                            
                            if all(k_slc_tmp_tmp == 0, 1) 
                                % if this gen depends on no other traj. params, 
                                % shift and scale, and
                                % also remove "full sliceability", cuz we already sliced (just not all the way).
                                c_sliced = c_sliced + (G_sliced_lb + G_sliced_ub)/2; % shift center 
                                G_sliced(:, slice_gen_idxs(j)) = (G_sliced_ub - G_sliced_lb)/2; % scale generator
                                k_slc_tmp(k_idx, slice_gen_idxs(j)) = 0; % remove indeterminate coefficients associated with parameter.
                                fully_slc_tmp(slice_gen_idxs(j)) = 0;
                            else
                                % wow these if else statements suck.
                                % this gen still depends on another param.
                                % can't shift, but can scale gen according
                                % to G_sliced_ub.
                                G_sliced(:, slice_gen_idxs(j)) = G_sliced_ub;
                                k_slc_tmp(k_idx, slice_gen_idxs(j)) = 0; % remove indeterminate coefficients associated with parameter.
                            end
                        end
                    end
                end
            end
            
            % take the fully sliced generators, add to center
            % see Alg. 1, line 11
            fully_slcd_idxs = fully_slc_tmp == 1 & all(k_slc_tmp == 0, 1);
            c_out = c_sliced + sum(G_sliced(:, fully_slcd_idxs), 2);
            G_out = G_sliced(:, ~fully_slcd_idxs);

            % remove fully sliced generators
            fully_slc_tmp = fully_slc_tmp(~fully_slcd_idxs);
            k_slc_tmp = k_slc_tmp(:, ~fully_slcd_idxs);
            gen_tracker_tmp = gen_tracker_tmp(:, ~fully_slcd_idxs);

            % remove rows associated with k_names
            k_list_tmp = obj.k_list;
            c_k_tmp = obj.c_k;
            delta_k_tmp = obj.delta_k;
            for i = 1:length(k_names)
                k_idx = find(strcmp(k_names{i}, obj.k_list));
                k_list_tmp(k_idx) = [];
                k_slc_tmp(k_idx, :) = [];
                c_k_tmp(k_idx) = [];
                delta_k_tmp(k_idx) = [];
            end
            
            % save
            obj.Vit = [c_out, G_out];
            obj.fully_slc = fully_slc_tmp;
            obj.k_slc = k_slc_tmp;
            obj.gen_tracker = gen_tracker_tmp;
            obj.k_list = k_list_tmp;
            obj.c_k = c_k_tmp;
            obj.delta_k = delta_k_tmp;


        end

        function obj = combine_generators(obj)
            % combine generators that have the same coefficients!!
            c_tmp = obj.Vit(:, 1);
            G_tmp = obj.Vit(:, 2:end);
            fully_slc_tmp = obj.fully_slc;
            k_slc_tmp = obj.k_slc;
            gen_tracker_tmp = obj.gen_tracker;
            n_vec = size(G_tmp, 2);
            % loop through generators...
            % check for generators whose coefficients are exactly the
            % same... these can be combined into one!! must be fully
            % sliceable and have the same k_slc values.
            
            % old:
%             i = 1;
%             while i < n_vec
%                 if fully_slc_tmp(i) == 0
%                     i = i+1;
%                     continue;
%                 end
%                 curr_k_slc = k_slc_tmp(:, i);
%                 duplicate_idxs = find(fully_slc_tmp((i+1):end) == 1 & all(k_slc_tmp(:, (i+1):end) == curr_k_slc, 1));
%                 if ~isempty(duplicate_idxs)
%                     duplicate_idxs = i + duplicate_idxs;
%                     % found generator(s) whose coefficients match
%                     G_tmp(:, i) = G_tmp(:, i) + sum(G_tmp(:, duplicate_idxs), 2);
%                     G_tmp(:, duplicate_idxs) = [];
%                     fully_slc_tmp(:, duplicate_idxs) = [];
%                     k_slc_tmp(:, duplicate_idxs) = [];
%                     n_vec = size(G_tmp, 2);
%                 end
%                 i = i+1;
%             end
%             obj.Vit = [c_tmp, G_tmp];
%             obj.fully_slc = fully_slc_tmp;
%             obj.k_slc = k_slc_tmp;

            % new:
            i = 1;
            while i < n_vec
                curr_gen_track = gen_tracker_tmp(:, i);
                duplicate_idxs = find(all(gen_tracker_tmp(:, (i+1):end) == curr_gen_track, 1));
                if ~isempty(duplicate_idxs)
                    duplicate_idxs = i + duplicate_idxs;
                    % found generator(s) whose coefficients match
                    G_tmp(:, i) = G_tmp(:, i) + sum(G_tmp(:, duplicate_idxs), 2);
                    G_tmp(:, duplicate_idxs) = [];
                    fully_slc_tmp(:, duplicate_idxs) = [];
                    k_slc_tmp(:, duplicate_idxs) = [];
                    gen_tracker_tmp(:, duplicate_idxs) = [];
                    n_vec = size(G_tmp, 2);
                end
                i = i+1;
            end
            obj.Vit = [c_tmp, G_tmp];
            obj.fully_slc = fully_slc_tmp;
            obj.k_slc = k_slc_tmp;
            obj.gen_tracker = gen_tracker_tmp;
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
            % %Shannon added: back down to 2!!
            % if obj.dimension == 2
            %     Z = zonotope([obj.Vit(1:2, :), buffer/2*eye(obj.dimension)]);
            % else
            %     Z = zonotope([obj.Vit, buffer/2*eye(obj.dimension)]);
            % end

            % check for indices where the coefficient can only be in [0, 1]
            % after slicing...
            c_tmp = obj.Vit(:, 1);
            G_tmp = obj.Vit(:, 2:end);
            %Shannon added: back down to 2!!
            if obj.dimension == 2
                c_tmp = c_tmp(1:2, :);
                G_tmp = G_tmp(1:2, :);
            end
            idxs = (obj.fully_slc == 1) & ~all(obj.k_slc == 0, 1) & all(mod(obj.k_slc, 2) == 0, 1);
            c_tmp = c_tmp + 0.5*sum(G_tmp(:, idxs), 2);
            G_tmp(:, idxs) = 0.5*G_tmp(:, idxs);
            
            Z = zonotope([c_tmp, G_tmp]);

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
        
        function [p] = plot_slice(obj, k_names, k_values, color, buffer, p)
            if length(k_names) ~= length(k_values)
               error('length of k_names doesn''t match length of k_values'); 
            end
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
            
            Z = obj.slice(k_names, k_values, fully_slc_only);
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