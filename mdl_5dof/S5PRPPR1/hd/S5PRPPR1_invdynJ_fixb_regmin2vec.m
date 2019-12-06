% Convert inverse dynamics regressor matrix to a vector for
% S5PRPPR1
% Use sparsity of the regressor matrix: 46/(5*19) elements are non-zero
%
% Input:
% RM [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [46x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5PRPPR1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(2, 2); RM(2, 3); RM(2, 4); RM(2, 5); RM(3, 5); RM(2, 6); RM(3, 6); RM(2, 7); RM(3, 7); RM(1, 8); RM(2, 8); RM(3, 8); RM(2, 9); RM(3, 9); RM(4, 9); RM(2, 10); RM(3, 10); RM(4, 10); RM(2, 11); RM(3, 11); RM(4, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(4, 12); RM(2, 13); RM(5, 13); RM(2, 14); RM(5, 14); RM(2, 15); RM(5, 15); RM(2, 16); RM(5, 16); RM(2, 17); RM(5, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(5, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(5, 19);];
RV = t1;
