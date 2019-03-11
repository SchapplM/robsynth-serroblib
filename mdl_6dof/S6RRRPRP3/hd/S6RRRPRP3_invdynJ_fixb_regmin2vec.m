% Convert inverse dynamics regressor matrix to a vector for
% S6RRRPRP3
% Use sparsity of the regressor matrix: 106/(6*32) elements are non-zero
%
% Input:
% RM [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [106x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S6RRRPRP3_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(2, 8); RM(1, 9); RM(2, 9); RM(1, 10); RM(2, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(2, 15); RM(3, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(1, 17); RM(2, 17); RM(3, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(1, 20); RM(2, 20); RM(3, 20); RM(4, 20); RM(1, 21); RM(2, 21); RM(3, 21); RM(4, 21); RM(1, 22); RM(2, 22); RM(3, 22); RM(5, 22); RM(1, 23); RM(2, 23); RM(3, 23); RM(5, 23); RM(1, 24); RM(2, 24); RM(3, 24); RM(5, 24); RM(1, 25); RM(2, 25); RM(3, 25); RM(5, 25); RM(1, 26); RM(2, 26); RM(3, 26); RM(5, 26); RM(1, 27); RM(2, 27); RM(3, 27); RM(4, 27); RM(5, 27); RM(1, 28); RM(2, 28); RM(3, 28); RM(4, 28); RM(5, 28); RM(1, 29); RM(2, 29); RM(3, 29); RM(4, 29); RM(5, 29); RM(6, 29); RM(1, 30); RM(2, 30); RM(3, 30); RM(4, 30); RM(5, 30); RM(6, 30); RM(1, 31); RM(2, 31); RM(3, 31); RM(4, 31); RM(5, 31); RM(6, 31); RM(1, 32); RM(2, 32); RM(3, 32); RM(4, 32); RM(5, 32); RM(6, 32);];
RV  = t1;
