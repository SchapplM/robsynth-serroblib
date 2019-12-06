% Convert inverse dynamics regressor matrix to a vector for
% S5PPRPR3
% Use sparsity of the regressor matrix: 33/(5*13) elements are non-zero
%
% Input:
% RM [5x13]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [33x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5PPRPR3_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(3, 3); RM(1, 4); RM(2, 4); RM(3, 4); RM(1, 5); RM(2, 5); RM(3, 5); RM(1, 6); RM(2, 6); RM(3, 6); RM(4, 6); RM(3, 7); RM(5, 7); RM(3, 8); RM(5, 8); RM(3, 9); RM(5, 9); RM(3, 10); RM(5, 10); RM(5, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(4, 12); RM(5, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(4, 13); RM(5, 13);];
RV = t1;
