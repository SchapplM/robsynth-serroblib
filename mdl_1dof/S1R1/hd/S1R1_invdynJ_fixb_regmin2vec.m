% Convert inverse dynamics regressor matrix to a vector for
% S1R1
% Use sparsity of the regressor matrix: 3/(1*3) elements are non-zero
%
% Input:
% RM [1x3]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [3x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2020-06-19 09:13
% Revision: caa0dbda1e8a16d11faaa29ba3bbef6afcd619f7 (2020-05-25)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S1R1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(1, 3);];
RV = t1;
