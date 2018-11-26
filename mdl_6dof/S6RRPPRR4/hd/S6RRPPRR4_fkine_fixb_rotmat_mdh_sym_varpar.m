% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2018-11-23 16:51
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T_c_mdh = S6RRPPRR4_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 16:50:43
% EndTime: 2018-11-23 16:50:43
% DurationCPUTime: 0.22s
% Computational Cost: add. (606->79), mult. (445->86), div. (0->0), fcn. (473->20), ass. (0->58)
t47 = pkin(6) - qJ(2);
t35 = cos(t47) / 0.2e1;
t46 = pkin(6) + qJ(2);
t42 = cos(t46);
t80 = t35 - t42 / 0.2e1;
t34 = sin(t46) / 0.2e1;
t40 = sin(t47);
t22 = t34 - t40 / 0.2e1;
t48 = sin(pkin(6));
t54 = sin(qJ(1));
t36 = t54 * t48;
t58 = cos(qJ(1));
t77 = t58 * t48;
t76 = pkin(7) + 0;
t45 = qJ(2) + pkin(11);
t50 = pkin(8) + qJ(3);
t15 = t22 * pkin(2) - t48 * t50;
t57 = cos(qJ(2));
t37 = t57 * pkin(2) + pkin(1);
t75 = t58 * t15 + t54 * t37 + 0;
t74 = pkin(6) - t45;
t73 = pkin(6) + t45;
t72 = -t54 * t15 + t58 * t37 + 0;
t71 = cos(t73);
t70 = sin(t74);
t49 = cos(pkin(6));
t69 = t80 * pkin(2) + t49 * t50 + t76;
t68 = cos(t74) / 0.2e1;
t67 = sin(t73) / 0.2e1;
t38 = sin(t45);
t65 = t71 / 0.2e1 + t68;
t8 = t54 * t38 - t58 * t65;
t41 = cos(t45);
t64 = t67 - t70 / 0.2e1;
t9 = t54 * t41 + t58 * t64;
t66 = t9 * pkin(3) + t8 * qJ(4) + t75;
t10 = t58 * t38 + t54 * t65;
t11 = t58 * t41 - t54 * t64;
t63 = t11 * pkin(3) + t10 * qJ(4) + t72;
t20 = t70 / 0.2e1 + t67;
t21 = t68 - t71 / 0.2e1;
t62 = t21 * pkin(3) - t20 * qJ(4) + t69;
t61 = pkin(4) * t36 + t11 * pkin(9) + t63;
t60 = t49 * pkin(4) + t21 * pkin(9) + t62;
t59 = -pkin(4) * t77 + t9 * pkin(9) + t66;
t56 = cos(qJ(5));
t55 = cos(qJ(6));
t53 = sin(qJ(2));
t52 = sin(qJ(5));
t51 = sin(qJ(6));
t23 = t35 + t42 / 0.2e1;
t13 = -t20 * t52 + t49 * t56;
t12 = t20 * t56 + t49 * t52;
t4 = t8 * t52 - t56 * t77;
t3 = t52 * t77 + t8 * t56;
t2 = t10 * t52 + t56 * t36;
t1 = -t10 * t56 + t52 * t36;
t5 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t58, -t54, 0, 0; t54, t58, 0, 0; 0, 0, 1, t76; 0, 0, 0, 1; -t54 * t22 + t58 * t57, -t54 * t23 - t58 * t53, t36, t58 * pkin(1) + pkin(8) * t36 + 0; t58 * t22 + t54 * t57, t58 * t23 - t54 * t53, -t77, t54 * pkin(1) - pkin(8) * t77 + 0; t80, t34 + t40 / 0.2e1, t49, t49 * pkin(8) + t76; 0, 0, 0, 1; t11, -t10, t36, t72; t9, -t8, -t77, t75; t21, t20, t49, t69; 0, 0, 0, 1; t36, -t11, t10, t63; -t77, -t9, t8, t66; t49, -t21, -t20, t62; 0, 0, 0, 1; t2, -t1, t11, t61; t4, t3, t9, t59; t13, -t12, t21, t60; 0, 0, 0, 1; t11 * t51 + t2 * t55, t11 * t55 - t2 * t51, t1, t2 * pkin(5) + t1 * pkin(10) + t61; t4 * t55 + t9 * t51, -t4 * t51 + t9 * t55, -t3, t4 * pkin(5) - t3 * pkin(10) + t59; t13 * t55 + t21 * t51, -t13 * t51 + t21 * t55, t12, t13 * pkin(5) + t12 * pkin(10) + t60; 0, 0, 0, 1;];
T_ges = t5;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
