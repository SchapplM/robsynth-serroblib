% Convert inverse dynamics regressor matrix to a vector for
% S5RPRRP12
% Use sparsity of the regressor matrix: 67/(5*24) elements are non-zero
%
% Input:
% RM [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [67x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5RPRRP12_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(3, 7); RM(1, 8); RM(3, 8); RM(1, 9); RM(3, 9); RM(1, 10); RM(3, 10); RM(3, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(1, 14); RM(3, 14); RM(4, 14); RM(1, 15); RM(3, 15); RM(4, 15); RM(1, 16); RM(3, 16); RM(4, 16); RM(1, 17); RM(3, 17); RM(4, 17); RM(1, 18); RM(3, 18); RM(4, 18); RM(1, 19); RM(2, 19); RM(3, 19); RM(4, 19); RM(1, 20); RM(2, 20); RM(3, 20); RM(4, 20); RM(1, 21); RM(2, 21); RM(3, 21); RM(4, 21); RM(5, 21); RM(1, 22); RM(2, 22); RM(3, 22); RM(4, 22); RM(5, 22); RM(1, 23); RM(2, 23); RM(3, 23); RM(4, 23); RM(5, 23); RM(1, 24); RM(2, 24); RM(3, 24); RM(4, 24); RM(5, 24);];
RV = t1;
