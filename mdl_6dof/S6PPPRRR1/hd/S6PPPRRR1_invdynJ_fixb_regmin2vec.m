% Convert inverse dynamics regressor matrix to a vector for
% S6PPPRRR1
% Use sparsity of the regressor matrix: 61/(6*20) elements are non-zero
%
% Input:
% RM [6x20]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [61x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S6PPPRRR1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(1, 3); RM(2, 3); RM(3, 3); RM(4, 4); RM(1, 5); RM(2, 5); RM(3, 5); RM(4, 5); RM(1, 6); RM(2, 6); RM(3, 6); RM(4, 6); RM(4, 7); RM(5, 7); RM(4, 8); RM(5, 8); RM(4, 9); RM(5, 9); RM(4, 10); RM(5, 10); RM(5, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(4, 12); RM(5, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(4, 13); RM(5, 13); RM(4, 14); RM(5, 14); RM(6, 14); RM(4, 15); RM(5, 15); RM(6, 15); RM(4, 16); RM(5, 16); RM(6, 16); RM(4, 17); RM(5, 17); RM(6, 17); RM(4, 18); RM(5, 18); RM(6, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(5, 19); RM(6, 19); RM(1, 20); RM(2, 20); RM(3, 20); RM(4, 20); RM(5, 20); RM(6, 20);];
RV  = t1;
