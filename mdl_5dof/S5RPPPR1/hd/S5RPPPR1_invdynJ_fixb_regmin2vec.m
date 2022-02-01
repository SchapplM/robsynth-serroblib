% Convert inverse dynamics regressor matrix to a vector for
% S5RPPPR1
% Use sparsity of the regressor matrix: 42/(5*17) elements are non-zero
%
% Input:
% RM [5x17]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [42x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5RPPPR1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(3, 5); RM(1, 6); RM(3, 6); RM(1, 7); RM(2, 7); RM(3, 7); RM(1, 8); RM(3, 8); RM(4, 8); RM(1, 9); RM(3, 9); RM(4, 9); RM(1, 10); RM(2, 10); RM(3, 10); RM(4, 10); RM(1, 11); RM(5, 11); RM(1, 12); RM(5, 12); RM(1, 13); RM(5, 13); RM(1, 14); RM(5, 14); RM(1, 15); RM(5, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(4, 16); RM(5, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(4, 17); RM(5, 17);];
RV = t1;
