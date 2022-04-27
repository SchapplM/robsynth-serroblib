% Calculate homogenous joint transformation matrices for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)
% T_stack [(5+1)*3 x 4]
%   stacked matrices from T_mdh into one 2D array, last row left out.
%   Last row only contains [0 0 0 1].

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-31 21:42
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [T_mdh, T_stack] = S5PPPRR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-31 21:42:26
% EndTime: 2022-01-31 21:42:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (9->9), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t62 = cos(qJ(4));
t61 = cos(qJ(5));
t60 = sin(qJ(4));
t59 = sin(qJ(5));
t58 = cos(pkin(7));
t57 = cos(pkin(8));
t56 = cos(pkin(9));
t55 = sin(pkin(7));
t54 = sin(pkin(8));
t53 = sin(pkin(9));
t1 = [t58, -t55, 0, 0; t55, t58, 0, 0; 0, 0, 1, qJ(1); t57, -t54, 0, pkin(1); 0, 0, -1, -qJ(2); t54, t57, 0, 0; t56, -t53, 0, pkin(2); 0, 0, -1, -qJ(3); t53, t56, 0, 0; t62, -t60, 0, pkin(3); 0, 0, -1, -pkin(5); t60, t62, 0, 0; t61, -t59, 0, pkin(4); 0, 0, -1, -pkin(6); t59, t61, 0, 0;];
T_stack = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = [T_stack((i-1)*3+1 : 3*i, :);[0 0 0 1]];
end
