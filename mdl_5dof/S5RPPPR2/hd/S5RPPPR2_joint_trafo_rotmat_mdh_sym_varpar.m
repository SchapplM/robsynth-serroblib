% Calculate homogenous joint transformation matrices for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5RPPPR2_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:11
% EndTime: 2020-01-03 11:22:11
% DurationCPUTime: 0.04s
% Computational Cost: add. (9->9), mult. (0->0), div. (0->0), fcn. (20->10), ass. (0->11)
t56 = cos(qJ(1));
t55 = cos(qJ(5));
t54 = sin(qJ(1));
t53 = sin(qJ(5));
t52 = cos(pkin(7));
t51 = cos(pkin(8));
t50 = cos(pkin(9));
t49 = sin(pkin(7));
t48 = sin(pkin(8));
t47 = sin(pkin(9));
t1 = [0, 0, 1, pkin(5); t54, t56, 0, 0; -t56, t54, 0, 0; 0, 0, 0, 1; t52, -t49, 0, pkin(1); 0, 0, -1, -qJ(2); t49, t52, 0, 0; 0, 0, 0, 1; t51, -t48, 0, pkin(2); 0, 0, -1, -qJ(3); t48, t51, 0, 0; 0, 0, 0, 1; t50, -t47, 0, pkin(3); 0, 0, -1, -qJ(4); t47, t50, 0, 0; 0, 0, 0, 1; t55, -t53, 0, pkin(4); 0, 0, -1, -pkin(6); t53, t55, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
