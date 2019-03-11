% Convert inverse dynamics regressor matrix to a vector for
% S6PRPRRR2
% Use sparsity of the regressor matrix: 83/(6*26) elements are non-zero
%
% Input:
% RM [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [83x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S6PRPRRR2_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(2, 2); RM(1, 3); RM(2, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(3, 5); RM(2, 6); RM(4, 6); RM(2, 7); RM(4, 7); RM(2, 8); RM(4, 8); RM(2, 9); RM(4, 9); RM(4, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(4, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(4, 12); RM(2, 13); RM(4, 13); RM(5, 13); RM(2, 14); RM(4, 14); RM(5, 14); RM(2, 15); RM(4, 15); RM(5, 15); RM(2, 16); RM(4, 16); RM(5, 16); RM(2, 17); RM(4, 17); RM(5, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(5, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(5, 19); RM(2, 20); RM(4, 20); RM(5, 20); RM(6, 20); RM(2, 21); RM(4, 21); RM(5, 21); RM(6, 21); RM(2, 22); RM(4, 22); RM(5, 22); RM(6, 22); RM(2, 23); RM(4, 23); RM(5, 23); RM(6, 23); RM(2, 24); RM(4, 24); RM(5, 24); RM(6, 24); RM(1, 25); RM(2, 25); RM(3, 25); RM(4, 25); RM(5, 25); RM(6, 25); RM(1, 26); RM(2, 26); RM(3, 26); RM(4, 26); RM(5, 26); RM(6, 26);];
RV  = t1;
