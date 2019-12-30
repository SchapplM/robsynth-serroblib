% Calculate homogenous joint transformation matrices for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% T_mdh [4x4x4]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-29 12:35
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S4PRRR7_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-29 12:33:44
% EndTime: 2019-12-29 12:33:44
% DurationCPUTime: 0.06s
% Computational Cost: add. (8->8), mult. (6->6), div. (0->0), fcn. (26->10), ass. (0->11)
t54 = cos(qJ(2));
t53 = cos(qJ(3));
t52 = cos(qJ(4));
t51 = sin(qJ(2));
t50 = sin(qJ(3));
t49 = sin(qJ(4));
t48 = cos(pkin(4));
t47 = cos(pkin(8));
t46 = sin(pkin(4));
t45 = sin(pkin(8));
t1 = [t47, -t45, 0, 0; t45, t47, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t54, -t51, 0, pkin(1); t48 * t51, t48 * t54, -t46, -t46 * pkin(5); t46 * t51, t46 * t54, t48, t48 * pkin(5); 0, 0, 0, 1; t53, -t50, 0, pkin(2); 0, 0, -1, -pkin(6); t50, t53, 0, 0; 0, 0, 0, 1; t52, -t49, 0, pkin(3); 0, 0, -1, -pkin(7); t49, t52, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,4);             % numerisch
else,                         T_mdh = sym('xx', [4,4,4]); end % symbolisch

for i = 1:4
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
