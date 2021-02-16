% Convert inverse dynamics regressor matrix to a vector for
% S1P1
% Use sparsity of the regressor matrix: 1/(1*1) elements are non-zero
%
% Input:
% RM [1x1]
%   minimal parameter regressor of inverse dynamics joint torque vector
%
% Output:
% RV [1x1]
%   vector of non-Null entries of the input matrix. (columns, then rows).

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RV = S1P1_invdynJ_fixb_regmin2vec(RM)

% From invdyn_joint_fixb_regressor_minpar_occupancy_vector_matlab.m
t1 = [RM(1, 1);];
RV = t1;
