% Convert inverse dynamics regressor matrix to a vector for
% S5PRRPR8
% Use sparsity of the regressor matrix: 41/(5*15) elements are non-zero
%
% Input:
% RM [5x15]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [41x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5PRRPR8_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(2, 2); RM(1, 3); RM(2, 3); RM(1, 4); RM(2, 4); RM(2, 5); RM(3, 5); RM(1, 6); RM(2, 6); RM(3, 6); RM(1, 7); RM(2, 7); RM(3, 7); RM(1, 8); RM(2, 8); RM(3, 8); RM(4, 8); RM(2, 9); RM(3, 9); RM(5, 9); RM(2, 10); RM(3, 10); RM(5, 10); RM(2, 11); RM(3, 11); RM(5, 11); RM(2, 12); RM(3, 12); RM(5, 12); RM(5, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(4, 14); RM(5, 14); RM(1, 15); RM(2, 15); RM(3, 15); RM(4, 15); RM(5, 15);];
RV = t1;
