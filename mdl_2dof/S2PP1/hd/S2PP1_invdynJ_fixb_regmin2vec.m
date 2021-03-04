% Convert inverse dynamics regressor matrix to a vector for
% S2PP1
% Use sparsity of the regressor matrix: 3/(2*2) elements are non-zero
%
% Input:
% RM [2x2]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [3x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2021-03-03 18:41
% Revision: 33b345ae0dd6ec4aa15499ab3d43edbbded0bea5 (2021-02-02)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S2PP1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1); RM(1, 2); RM(2, 2);];
RV = t1;
