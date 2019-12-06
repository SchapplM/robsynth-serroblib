% Calculate homogenous joint transformation matrices for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% T_mdh [4x4x5]
%   homogenous transformation matrices for joint transformation (MDH)
%   Transformation matrices from one joint to the next (not: from base to joints)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T_mdh = S5PRPRR6_joint_trafo_rotmat_mdh_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_joint_trafo_rotmat_mdh_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_joint_trafo_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From joint_transformation_mdh_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:56:18
% EndTime: 2019-12-05 15:56:18
% DurationCPUTime: 0.03s
% Computational Cost: add. (9->9), mult. (6->6), div. (0->0), fcn. (30->12), ass. (0->13)
t68 = cos(qJ(2));
t67 = cos(qJ(4));
t66 = cos(qJ(5));
t65 = sin(qJ(2));
t64 = sin(qJ(4));
t63 = sin(qJ(5));
t62 = cos(pkin(5));
t61 = cos(pkin(9));
t60 = cos(pkin(10));
t59 = sin(pkin(5));
t58 = sin(pkin(9));
t57 = sin(pkin(10));
t1 = [t61, -t58, 0, 0; t58, t61, 0, 0; 0, 0, 1, qJ(1); 0, 0, 0, 1; t68, -t65, 0, pkin(1); t62 * t65, t62 * t68, -t59, -t59 * pkin(6); t59 * t65, t59 * t68, t62, t62 * pkin(6); 0, 0, 0, 1; t60, -t57, 0, pkin(2); 0, 0, -1, -qJ(3); t57, t60, 0, 0; 0, 0, 0, 1; t67, -t64, 0, pkin(3); t64, t67, 0, 0; 0, 0, 1, pkin(7); 0, 0, 0, 1; t66, -t63, 0, pkin(4); 0, 0, -1, -pkin(8); t63, t66, 0, 0; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_mdh = NaN(4,4,5);             % numerisch
else,                         T_mdh = sym('xx', [4,4,5]); end % symbolisch

for i = 1:5
  T_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
