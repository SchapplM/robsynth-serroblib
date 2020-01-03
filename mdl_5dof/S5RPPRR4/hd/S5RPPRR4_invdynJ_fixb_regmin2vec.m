% Convert inverse dynamics regressor matrix to a vector for
% S5RPPRR4
% Use sparsity of the regressor matrix: 66/(5*25) elements are non-zero
%
% Input:
% RM [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [66x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5RPPRR4_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(1, 8); RM(2, 8); RM(3, 8); RM(1, 9); RM(2, 9); RM(3, 9); RM(1, 10); RM(2, 10); RM(3, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(1, 12); RM(4, 12); RM(1, 13); RM(4, 13); RM(1, 14); RM(4, 14); RM(1, 15); RM(4, 15); RM(1, 16); RM(4, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(4, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(1, 19); RM(4, 19); RM(5, 19); RM(1, 20); RM(4, 20); RM(5, 20); RM(1, 21); RM(4, 21); RM(5, 21); RM(1, 22); RM(4, 22); RM(5, 22); RM(1, 23); RM(4, 23); RM(5, 23); RM(1, 24); RM(2, 24); RM(3, 24); RM(4, 24); RM(5, 24); RM(1, 25); RM(2, 25); RM(3, 25); RM(4, 25); RM(5, 25);];
RV = t1;
