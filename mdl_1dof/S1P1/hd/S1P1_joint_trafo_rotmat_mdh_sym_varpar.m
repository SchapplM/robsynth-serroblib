% Calculate homogenous joint transformation matrices for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% T_mdh [4x4x1]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(1+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S1P1_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [1x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:46
% EndTime: 2021-01-14 12:21:46
% DurationCPUTime: 0.02s
% Computational Cost: add. (0->0), mult. (0->0), div. (0->0), fcn. (0->0), ass. (0->1)
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, qJ(1);];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,1);             % numerisch
else,                         T_mdh = sym('xx', [4,4,1]); end % symbolisch

for i = 1:1
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
