% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% T_c_mdh [4x4x(6+1)]
%   homogenous transformation matrices for each (body) frame (MDH)
%   1:  mdh base (link 0) -> mdh base link 0 (unit matrix, no information)
%   ...
%   7:  mdh base (link 0) -> mdh frame (7-1), link (7-1)
%   ...
%   6+1:  mdh base (link 0) -> mdh frame (6)

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:09
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:09:27
% EndTime: 2018-11-23 15:09:28
% DurationCPUTime: 0.16s
% Computational Cost: add. (552->74), mult. (621->70), div. (0->0), fcn. (697->14), ass. (0->55)
t70 = cos(qJ(3));
t33 = sin(pkin(10));
t34 = sin(pkin(6));
t69 = t33 * t34;
t35 = cos(pkin(10));
t68 = t35 * t34;
t60 = pkin(6) + qJ(2);
t50 = sin(t60) / 0.2e1;
t61 = pkin(6) - qJ(2);
t54 = sin(t61);
t22 = t50 - t54 / 0.2e1;
t40 = cos(qJ(2));
t14 = t35 * t22 + t33 * t40;
t37 = sin(qJ(3));
t58 = t34 * t70;
t5 = t14 * t37 + t35 * t58;
t67 = t5 * qJ(4);
t16 = -t33 * t22 + t35 * t40;
t7 = t16 * t37 - t33 * t58;
t66 = t7 * qJ(4);
t65 = pkin(5) + qJ(4);
t64 = pkin(8) - qJ(5);
t51 = cos(t61) / 0.2e1;
t55 = cos(t60);
t23 = t51 - t55 / 0.2e1;
t62 = cos(pkin(6));
t17 = t23 * t37 - t62 * t70;
t63 = t17 * qJ(4);
t59 = qJ(1) + 0;
t57 = t35 * pkin(1) + pkin(7) * t69 + 0;
t56 = t62 * pkin(7) + t59;
t53 = t16 * pkin(2) + t57;
t52 = t23 * pkin(2) + t56;
t49 = t33 * pkin(1) - pkin(7) * t68 + 0;
t48 = t14 * pkin(2) + t49;
t38 = sin(qJ(2));
t44 = t51 + t55 / 0.2e1;
t15 = t33 * t44 + t35 * t38;
t47 = t15 * pkin(8) + t53;
t21 = t50 + t54 / 0.2e1;
t46 = -t21 * pkin(8) + t52;
t13 = t33 * t38 - t35 * t44;
t45 = t13 * pkin(8) + t48;
t8 = t16 * t70 + t37 * t69;
t4 = t8 * pkin(3);
t43 = t8 * pkin(4) + t64 * t15 + t4 + t53;
t18 = t23 * t70 + t62 * t37;
t12 = t18 * pkin(3);
t42 = t18 * pkin(4) - t64 * t21 + t12 + t52;
t6 = t14 * t70 - t37 * t68;
t2 = t6 * pkin(3);
t41 = t6 * pkin(4) + t64 * t13 + t2 + t48;
t39 = cos(qJ(6));
t36 = sin(qJ(6));
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t35, -t33, 0, 0; t33, t35, 0, 0; 0, 0, 1, t59; 0, 0, 0, 1; t16, -t15, t69, t57; t14, -t13, -t68, t49; t23, t21, t62, t56; 0, 0, 0, 1; t8, -t7, t15, t47; t6, -t5, t13, t45; t18, -t17, -t21, t46; 0, 0, 0, 1; t8, t15, t7, t4 + t47 + t66; t6, t13, t5, t2 + t45 + t67; t18, -t21, t17, t12 + t46 + t63; 0, 0, 0, 1; t7, -t8, -t15, t43 + t66; t5, -t6, -t13, t41 + t67; t17, -t18, t21, t42 + t63; 0, 0, 0, 1; -t15 * t36 + t7 * t39, -t15 * t39 - t7 * t36, t8, t8 * pkin(9) + t65 * t7 + t43; -t13 * t36 + t5 * t39, -t13 * t39 - t5 * t36, t6, t6 * pkin(9) + t65 * t5 + t41; t17 * t39 + t21 * t36, -t17 * t36 + t21 * t39, t18, t18 * pkin(9) + t65 * t17 + t42; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
