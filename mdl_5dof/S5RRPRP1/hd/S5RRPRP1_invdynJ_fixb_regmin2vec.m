% Convert inverse dynamics regressor matrix to a vector for
% S5RRPRP1
% Use sparsity of the regressor matrix: 52/(5*18) elements are non-zero
%
% Input:
% RM [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [52x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S5RRPRP1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3); RM(1, 4); RM(2, 4); RM(1, 5); RM(2, 5); RM(1, 6); RM(2, 6); RM(1, 7); RM(2, 7); RM(3, 7); RM(1, 8); RM(2, 8); RM(4, 8); RM(1, 9); RM(2, 9); RM(4, 9); RM(1, 10); RM(2, 10); RM(4, 10); RM(1, 11); RM(2, 11); RM(4, 11); RM(4, 12); RM(1, 13); RM(2, 13); RM(3, 13); RM(4, 13); RM(1, 14); RM(2, 14); RM(3, 14); RM(4, 14); RM(1, 15); RM(2, 15); RM(3, 15); RM(4, 15); RM(5, 15); RM(1, 16); RM(2, 16); RM(3, 16); RM(4, 16); RM(5, 16); RM(1, 17); RM(2, 17); RM(4, 17); RM(5, 17); RM(1, 18); RM(2, 18); RM(3, 18); RM(4, 18); RM(5, 18);];
RV = t1;
