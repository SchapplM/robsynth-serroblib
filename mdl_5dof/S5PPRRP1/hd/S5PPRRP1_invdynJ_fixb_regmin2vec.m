% Convert inverse dynamics regressor matrix to a vector for
% S5PPRRP1
% Use sparsity of the regressor matrix: 44/(5*16) elements are non-zero
%
% Input:
% RM [5x16]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [44x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 14:48
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5PPRRP1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(3, 3); RM(1, 4); RM(3, 4); RM(1, 5); RM(3, 5); RM(3, 6); RM(4, 6); RM(3, 7); RM(4, 7); RM(3, 8); RM(4, 8); RM(3, 9); RM(4, 9); RM(4, 10); RM(1, 11); RM(2, 11); RM(3, 11); RM(4, 11); RM(1, 12); RM(2, 12); RM(3, 12); RM(4, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(4, 13); RM(5, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(4, 14); RM(5, 14); RM(1, 15); RM(3, 15); RM(4, 15); RM(5, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(4, 16); RM(5, 16);];
RV = t1;
