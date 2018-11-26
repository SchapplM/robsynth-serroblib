% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:04
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:03:36
% EndTime: 2018-11-23 17:03:36
% DurationCPUTime: 0.24s
% Computational Cost: add. (677->83), mult. (509->86), div. (0->0), fcn. (549->20), ass. (0->62)
t47 = pkin(6) - qJ(2);
t35 = cos(t47) / 0.2e1;
t46 = pkin(6) + qJ(2);
t43 = cos(t46);
t83 = t35 - t43 / 0.2e1;
t34 = sin(t46) / 0.2e1;
t41 = sin(t47);
t23 = t34 - t41 / 0.2e1;
t45 = qJ(2) + pkin(11);
t39 = sin(t45);
t54 = sin(qJ(1));
t58 = cos(qJ(1));
t73 = pkin(6) - t45;
t64 = cos(t73) / 0.2e1;
t72 = pkin(6) + t45;
t69 = cos(t72);
t59 = t69 / 0.2e1 + t64;
t9 = t54 * t39 - t58 * t59;
t82 = t9 * pkin(9);
t79 = pkin(5) + pkin(9);
t11 = t58 * t39 + t54 * t59;
t78 = t11 * pkin(9);
t63 = sin(t72) / 0.2e1;
t68 = sin(t73);
t21 = t68 / 0.2e1 + t63;
t77 = t21 * pkin(9);
t48 = sin(pkin(6));
t36 = t54 * t48;
t76 = t58 * t48;
t75 = pkin(7) + 0;
t50 = pkin(8) + qJ(3);
t17 = t23 * pkin(2) - t48 * t50;
t57 = cos(qJ(2));
t38 = t57 * pkin(2) + pkin(1);
t74 = t58 * t17 + t54 * t38 + 0;
t20 = t63 - t68 / 0.2e1;
t42 = cos(t45);
t10 = t58 * t20 + t54 * t42;
t71 = t10 * pkin(3) + t74;
t70 = -t54 * t17 + t58 * t38 + 0;
t49 = cos(pkin(6));
t67 = t83 * pkin(2) + t49 * t50 + t75;
t12 = -t54 * t20 + t58 * t42;
t66 = t12 * pkin(3) + t70;
t22 = t64 - t69 / 0.2e1;
t65 = t22 * pkin(3) + t67;
t52 = sin(qJ(4));
t56 = cos(qJ(4));
t3 = t10 * t52 + t56 * t76;
t4 = t10 * t56 - t52 * t76;
t62 = t4 * pkin(4) + t3 * qJ(5) + t71;
t5 = t12 * t52 - t56 * t36;
t6 = t12 * t56 + t52 * t36;
t61 = t6 * pkin(4) + t5 * qJ(5) + t66;
t14 = t22 * t52 - t49 * t56;
t15 = t22 * t56 + t49 * t52;
t60 = t15 * pkin(4) + t14 * qJ(5) + t65;
t55 = cos(qJ(6));
t53 = sin(qJ(2));
t51 = sin(qJ(6));
t24 = t35 + t43 / 0.2e1;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t75; 0, 0, 0, 1; -t54 * t23 + t58 * t57, -t54 * t24 - t58 * t53, t36, t58 * pkin(1) + pkin(8) * t36 + 0; t58 * t23 + t54 * t57, t58 * t24 - t54 * t53, -t76, t54 * pkin(1) - pkin(8) * t76 + 0; t83, t34 + t41 / 0.2e1, t49, t49 * pkin(8) + t75; 0, 0, 0, 1; t12, -t11, t36, t70; t10, -t9, -t76, t74; t22, t21, t49, t67; 0, 0, 0, 1; t6, -t5, t11, t66 + t78; t4, -t3, t9, t71 + t82; t15, -t14, -t21, t65 - t77; 0, 0, 0, 1; t11, -t6, t5, t61 + t78; t9, -t4, t3, t62 + t82; -t21, -t15, t14, t60 - t77; 0, 0, 0, 1; t11 * t55 + t5 * t51, -t11 * t51 + t5 * t55, t6, t6 * pkin(10) + t79 * t11 + t61; t3 * t51 + t9 * t55, t3 * t55 - t9 * t51, t4, t4 * pkin(10) + t79 * t9 + t62; t14 * t51 - t21 * t55, t14 * t55 + t21 * t51, t15, t15 * pkin(10) - t79 * t21 + t60; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
