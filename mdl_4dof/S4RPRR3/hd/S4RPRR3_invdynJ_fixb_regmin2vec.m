% Convert inverse dynamics regressor matrix to a vector for
% S4RPRR3
% Use sparsity of the regressor matrix: 42/(4*18) elements are non-zero
%
% Input:
% RM [4x18]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [42x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:49
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S4RPRR3_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(3, 5); RM(1, 6); RM(3, 6); RM(1, 7); RM(3, 7); RM(1, 8); RM(3, 8); RM(3, 9); RM(1, 10); RM(2, 10); RM(3, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(1, 12); RM(3, 12); RM(4, 12); RM(1, 13); RM(3, 13); RM(4, 13); RM(1, 14); RM(3, 14); RM(4, 14); RM(1, 15); RM(3, 15); RM(4, 15); RM(3, 16); RM(4, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(4, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18);];
RV = t1;
