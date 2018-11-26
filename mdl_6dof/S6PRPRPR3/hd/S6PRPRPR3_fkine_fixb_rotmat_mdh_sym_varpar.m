% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2018-11-23 14:56
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6PRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:56:18
% EndTime: 2018-11-23 14:56:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (677->83), mult. (509->86), div. (0->0), fcn. (549->20), ass. (0->62)
t47 = pkin(6) - qJ(2);
t36 = cos(t47) / 0.2e1;
t46 = pkin(6) + qJ(2);
t43 = cos(t46);
t83 = t36 - t43 / 0.2e1;
t35 = sin(t46) / 0.2e1;
t41 = sin(t47);
t23 = t35 - t41 / 0.2e1;
t45 = qJ(2) + pkin(11);
t39 = sin(t45);
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t73 = pkin(6) - t45;
t65 = cos(t73) / 0.2e1;
t72 = pkin(6) + t45;
t69 = cos(t72);
t59 = t69 / 0.2e1 + t65;
t9 = t48 * t39 - t50 * t59;
t82 = t9 * pkin(8);
t79 = pkin(5) + pkin(8);
t11 = t50 * t39 + t48 * t59;
t78 = t11 * pkin(8);
t64 = sin(t72) / 0.2e1;
t68 = sin(t73);
t21 = t68 / 0.2e1 + t64;
t77 = t21 * pkin(8);
t49 = sin(pkin(6));
t33 = t48 * t49;
t76 = t50 * t49;
t75 = qJ(1) + 0;
t52 = pkin(7) + qJ(3);
t17 = t23 * pkin(2) - t49 * t52;
t58 = cos(qJ(2));
t38 = t58 * pkin(2) + pkin(1);
t74 = t50 * t17 + t48 * t38 + 0;
t20 = t64 - t68 / 0.2e1;
t42 = cos(t45);
t10 = t50 * t20 + t48 * t42;
t71 = t10 * pkin(3) + t74;
t70 = -t48 * t17 + t50 * t38 + 0;
t51 = cos(pkin(6));
t67 = t83 * pkin(2) + t51 * t52 + t75;
t12 = -t48 * t20 + t50 * t42;
t66 = t12 * pkin(3) + t70;
t22 = t65 - t69 / 0.2e1;
t63 = t22 * pkin(3) + t67;
t54 = sin(qJ(4));
t57 = cos(qJ(4));
t3 = t10 * t54 + t57 * t76;
t4 = t10 * t57 - t54 * t76;
t62 = t4 * pkin(4) + t3 * qJ(5) + t71;
t5 = t12 * t54 - t57 * t33;
t6 = t12 * t57 + t54 * t33;
t61 = t6 * pkin(4) + t5 * qJ(5) + t66;
t14 = t22 * t54 - t51 * t57;
t15 = t22 * t57 + t51 * t54;
t60 = t15 * pkin(4) + t14 * qJ(5) + t63;
t56 = cos(qJ(6));
t55 = sin(qJ(2));
t53 = sin(qJ(6));
t24 = t36 + t43 / 0.2e1;
t1 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t50, -t48, 0, 0; t48, t50, 0, 0; 0, 0, 1, t75; 0, 0, 0, 1; -t48 * t23 + t50 * t58, -t48 * t24 - t50 * t55, t33, t50 * pkin(1) + pkin(7) * t33 + 0; t50 * t23 + t48 * t58, t50 * t24 - t48 * t55, -t76, t48 * pkin(1) - pkin(7) * t76 + 0; t83, t35 + t41 / 0.2e1, t51, t51 * pkin(7) + t75; 0, 0, 0, 1; t12, -t11, t33, t70; t10, -t9, -t76, t74; t22, t21, t51, t67; 0, 0, 0, 1; t6, -t5, t11, t66 + t78; t4, -t3, t9, t71 + t82; t15, -t14, -t21, t63 - t77; 0, 0, 0, 1; t11, -t6, t5, t61 + t78; t9, -t4, t3, t62 + t82; -t21, -t15, t14, t60 - t77; 0, 0, 0, 1; t11 * t56 + t5 * t53, -t11 * t53 + t5 * t56, t6, t6 * pkin(9) + t79 * t11 + t61; t3 * t53 + t9 * t56, t3 * t56 - t9 * t53, t4, t4 * pkin(9) + t79 * t9 + t62; t14 * t53 - t21 * t56, t14 * t56 + t21 * t53, t15, t15 * pkin(9) - t79 * t21 + t60; 0, 0, 0, 1;];
T_ges = t1;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
