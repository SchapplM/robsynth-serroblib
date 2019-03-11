% Convert inverse dynamics regressor matrix to a vector for
% S6PRPPRR2
% Use sparsity of the regressor matrix: 65/(6*22) elements are non-zero
%
% Input:
% RM [6x22]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [65x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S6PRPPRR2_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(2, 2); RM(1, 3); RM(2, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(3, 5); RM(1, 6); RM(2, 6); RM(4, 6); RM(1, 7); RM(2, 7); RM(4, 7); RM(1, 8); RM(2, 8); RM(3, 8); RM(4, 8); RM(2, 9); RM(5, 9); RM(2, 10); RM(5, 10); RM(2, 11); RM(5, 11); RM(2, 12); RM(5, 12); RM(5, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(4, 14); RM(5, 14); RM(1, 15); RM(2, 15); RM(3, 15); RM(4, 15); RM(5, 15); RM(2, 16); RM(5, 16); RM(6, 16); RM(2, 17); RM(5, 17); RM(6, 17); RM(2, 18); RM(5, 18); RM(6, 18); RM(2, 19); RM(5, 19); RM(6, 19); RM(2, 20); RM(5, 20); RM(6, 20); RM(1, 21); RM(2, 21); RM(3, 21); RM(4, 21); RM(5, 21); RM(6, 21); RM(1, 22); RM(2, 22); RM(3, 22); RM(4, 22); RM(5, 22); RM(6, 22);];
RV  = t1;
