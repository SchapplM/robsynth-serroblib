% Convert inverse dynamics regressor matrix to a vector for
% S3PPP1
% Use sparsity of the regressor matrix: 6/(3*3) elements are non-zero
%
% Input:
% RM [3x3]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [6x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S3PPP1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(1, 3); RM(2, 3); RM(3, 3);];
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(1, 3); RM(2, 3); RM(3, 3);];
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(1, 3); RM(2, 3); RM(3, 3);];
t1 = [RM(1, 1); RM(1, 2); RM(2, 2); RM(1, 3); RM(2, 3); RM(3, 3);];
RV  = t1;
