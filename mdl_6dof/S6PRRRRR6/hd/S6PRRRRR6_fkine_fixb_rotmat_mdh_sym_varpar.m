% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2018-11-23 15:36
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:36:06
% EndTime: 2018-11-23 15:36:06
% DurationCPUTime: 0.44s
% Computational Cost: add. (3259->102), mult. (3274->129), div. (0->0), fcn. (3324->30), ass. (0->92)
t78 = sin(pkin(14));
t81 = sin(pkin(6));
t120 = t78 * t81;
t77 = pkin(6) - qJ(2);
t66 = cos(t77) / 0.2e1;
t76 = pkin(6) + qJ(2);
t70 = cos(t76);
t58 = t66 + t70 / 0.2e1;
t82 = cos(pkin(14));
t90 = sin(qJ(2));
t47 = -t78 * t58 - t82 * t90;
t80 = sin(pkin(7));
t84 = cos(pkin(7));
t39 = t84 * t120 - t47 * t80;
t64 = sin(t76) / 0.2e1;
t68 = sin(t77);
t53 = t64 + t68 / 0.2e1;
t85 = cos(pkin(6));
t44 = -t53 * t80 + t85 * t84;
t74 = pkin(7) + qJ(3);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - qJ(3);
t67 = sin(t75);
t51 = t63 + t67 / 0.2e1;
t65 = cos(t75) / 0.2e1;
t69 = cos(t74);
t56 = t65 + t69 / 0.2e1;
t59 = t66 - t70 / 0.2e1;
t89 = sin(qJ(3));
t32 = t85 * t51 + t53 * t56 - t59 * t89;
t79 = sin(pkin(8));
t83 = cos(pkin(8));
t23 = -t32 * t79 + t44 * t83;
t54 = t64 - t68 / 0.2e1;
t95 = cos(qJ(2));
t48 = -t78 * t54 + t82 * t95;
t28 = t120 * t51 + t47 * t56 - t48 * t89;
t19 = -t28 * t79 + t39 * t83;
t119 = t82 * t81;
t45 = t82 * t58 - t78 * t90;
t46 = t82 * t54 + t78 * t95;
t26 = -t119 * t51 + t45 * t56 - t46 * t89;
t38 = -t119 * t84 - t45 * t80;
t18 = -t26 * t79 + t38 * t83;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t114 = qJ(1) + 0;
t113 = t82 * pkin(1) + pkin(9) * t120 + 0;
t112 = t85 * pkin(9) + t114;
t111 = cos(t116);
t110 = sin(t117);
t109 = cos(t117) / 0.2e1;
t108 = sin(t116) / 0.2e1;
t107 = t78 * pkin(1) - pkin(9) * t119 + 0;
t106 = t48 * pkin(2) + t39 * pkin(10) + t113;
t105 = t59 * pkin(2) + t44 * pkin(10) + t112;
t104 = t109 + t111 / 0.2e1;
t103 = t108 + t110 / 0.2e1;
t52 = t63 - t67 / 0.2e1;
t57 = t65 - t69 / 0.2e1;
t94 = cos(qJ(3));
t29 = t120 * t57 + t47 * t52 + t48 * t94;
t102 = t29 * pkin(3) + t19 * pkin(11) + t106;
t33 = t53 * t52 + t85 * t57 + t59 * t94;
t101 = t33 * pkin(3) + t23 * pkin(11) + t105;
t100 = t46 * pkin(2) + pkin(10) * t38 + t107;
t88 = sin(qJ(4));
t11 = -t103 * t39 - t104 * t28 + t29 * t88;
t50 = t108 - t110 / 0.2e1;
t55 = t109 - t111 / 0.2e1;
t93 = cos(qJ(4));
t12 = t28 * t50 + t29 * t93 + t39 * t55;
t99 = t12 * pkin(4) + t11 * pkin(12) + t102;
t14 = -t103 * t44 - t104 * t32 + t33 * t88;
t15 = t32 * t50 + t33 * t93 + t44 * t55;
t98 = t15 * pkin(4) + t14 * pkin(12) + t101;
t27 = -t119 * t57 + t45 * t52 + t46 * t94;
t97 = t27 * pkin(3) + t18 * pkin(11) + t100;
t10 = t26 * t50 + t27 * t93 + t38 * t55;
t9 = -t103 * t38 - t104 * t26 + t27 * t88;
t96 = t10 * pkin(4) + t9 * pkin(12) + t97;
t92 = cos(qJ(5));
t91 = cos(qJ(6));
t87 = sin(qJ(5));
t86 = sin(qJ(6));
t6 = t15 * t92 + t23 * t87;
t5 = t15 * t87 - t23 * t92;
t4 = t12 * t92 + t19 * t87;
t3 = t12 * t87 - t19 * t92;
t2 = t10 * t92 + t18 * t87;
t1 = t10 * t87 - t18 * t92;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t82, -t78, 0, 0; t78, t82, 0, 0; 0, 0, 1, t114; 0, 0, 0, 1; t48, t47, t120, t113; t46, t45, -t119, t107; t59, t53, t85, t112; 0, 0, 0, 1; t29, t28, t39, t106; t27, t26, t38, t100; t33, t32, t44, t105; 0, 0, 0, 1; t12, -t11, t19, t102; t10, -t9, t18, t97; t15, -t14, t23, t101; 0, 0, 0, 1; t4, -t3, t11, t99; t2, -t1, t9, t96; t6, -t5, t14, t98; 0, 0, 0, 1; t11 * t86 + t4 * t91, t11 * t91 - t4 * t86, t3, t4 * pkin(5) + t3 * pkin(13) + t99; t2 * t91 + t9 * t86, -t2 * t86 + t9 * t91, t1, t2 * pkin(5) + t1 * pkin(13) + t96; t14 * t86 + t6 * t91, t14 * t91 - t6 * t86, t5, t6 * pkin(5) + t5 * pkin(13) + t98; 0, 0, 0, 1;];
T_ges = t7;
%% Postprocessing: Reshape Output
% Convert Maple format (2-dimensional tensor) to Matlab format (3-dimensional tensor)
% Fallunterscheidung der Initialisierung für symbolische Eingabe
if isa([qJ; pkin], 'double'), T_c_mdh = NaN(4,4,6+1);               % numerisch
else,                         T_c_mdh = sym('xx', [4,4,6+1]); end % symbolisch
for i = 1:6+1
  T_c_mdh(:,:,i) = T_ges((i-1)*4+1 : 4*i, :);
end
Tc_stack = NaN(3*size(T_c_mdh,3),4);
% Zusätzliche Ausgabe: Als 2D-array gestapelt, ohne Zeile mit 0001
for i = 1:size(T_c_mdh,3), Tc_stack((i-1)*3+1:3*i,1:4) = T_c_mdh(1:3,1:4,i); end
