% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PRPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2018-11-23 15:07
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:07:05
% EndTime: 2018-11-23 15:07:05
% DurationCPUTime: 0.45s
% Computational Cost: add. (3259->104), mult. (3274->132), div. (0->0), fcn. (3324->30), ass. (0->95)
t74 = pkin(7) + pkin(14);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - pkin(14);
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t64 = cos(t75) / 0.2e1;
t68 = cos(t74);
t52 = t64 + t68 / 0.2e1;
t76 = pkin(6) + qJ(2);
t65 = sin(t76) / 0.2e1;
t77 = pkin(6) - qJ(2);
t69 = sin(t77);
t55 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t59 = t66 - t70 / 0.2e1;
t78 = sin(pkin(14));
t87 = cos(pkin(6));
t32 = t87 * t50 + t55 * t52 - t59 * t78;
t86 = cos(pkin(7));
t119 = t87 * t86;
t81 = sin(pkin(7));
t44 = -t55 * t81 + t119;
t80 = sin(pkin(8));
t85 = cos(pkin(8));
t23 = -t32 * t80 + t44 * t85;
t79 = sin(pkin(13));
t82 = sin(pkin(6));
t121 = t79 * t82;
t58 = t66 + t70 / 0.2e1;
t84 = cos(pkin(13));
t91 = sin(qJ(2));
t47 = -t79 * t58 - t84 * t91;
t56 = t65 - t69 / 0.2e1;
t95 = cos(qJ(2));
t48 = -t79 * t56 + t84 * t95;
t28 = t121 * t50 + t47 * t52 - t48 * t78;
t115 = t86 * t121;
t39 = -t47 * t81 + t115;
t19 = -t28 * t80 + t39 * t85;
t120 = t84 * t82;
t45 = t84 * t58 - t79 * t91;
t46 = t84 * t56 + t79 * t95;
t26 = -t120 * t50 + t45 * t52 - t46 * t78;
t38 = -t120 * t86 - t45 * t81;
t18 = -t26 * t80 + t38 * t85;
t118 = qJ(3) * t81;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t114 = qJ(1) + 0;
t113 = t84 * pkin(1) + pkin(9) * t121 + 0;
t112 = t87 * pkin(9) + t114;
t111 = cos(t116);
t110 = sin(t117);
t109 = cos(t117) / 0.2e1;
t108 = sin(t116) / 0.2e1;
t107 = t79 * pkin(1) - pkin(9) * t120 + 0;
t106 = t48 * pkin(2) + qJ(3) * t115 - t118 * t47 + t113;
t105 = t59 * pkin(2) + qJ(3) * t119 - t118 * t55 + t112;
t104 = t109 + t111 / 0.2e1;
t103 = t108 + t110 / 0.2e1;
t51 = t63 - t67 / 0.2e1;
t53 = t64 - t68 / 0.2e1;
t83 = cos(pkin(14));
t29 = t121 * t53 + t47 * t51 + t48 * t83;
t102 = t29 * pkin(3) + t19 * pkin(10) + t106;
t33 = t55 * t51 + t87 * t53 + t59 * t83;
t101 = t33 * pkin(3) + t23 * pkin(10) + t105;
t100 = t46 * pkin(2) + qJ(3) * t38 + t107;
t90 = sin(qJ(4));
t11 = -t103 * t39 - t104 * t28 + t29 * t90;
t54 = t108 - t110 / 0.2e1;
t57 = t109 - t111 / 0.2e1;
t94 = cos(qJ(4));
t12 = t28 * t54 + t29 * t94 + t39 * t57;
t99 = t12 * pkin(4) + t11 * pkin(11) + t102;
t14 = -t103 * t44 - t104 * t32 + t33 * t90;
t15 = t32 * t54 + t33 * t94 + t44 * t57;
t98 = t15 * pkin(4) + t14 * pkin(11) + t101;
t27 = -t120 * t53 + t45 * t51 + t46 * t83;
t97 = t27 * pkin(3) + t18 * pkin(10) + t100;
t10 = t26 * t54 + t27 * t94 + t38 * t57;
t9 = -t103 * t38 - t104 * t26 + t27 * t90;
t96 = t10 * pkin(4) + t9 * pkin(11) + t97;
t93 = cos(qJ(5));
t92 = cos(qJ(6));
t89 = sin(qJ(5));
t88 = sin(qJ(6));
t6 = t15 * t93 + t23 * t89;
t5 = t15 * t89 - t23 * t93;
t4 = t12 * t93 + t19 * t89;
t3 = t12 * t89 - t19 * t93;
t2 = t10 * t93 + t18 * t89;
t1 = t10 * t89 - t18 * t93;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t84, -t79, 0, 0; t79, t84, 0, 0; 0, 0, 1, t114; 0, 0, 0, 1; t48, t47, t121, t113; t46, t45, -t120, t107; t59, t55, t87, t112; 0, 0, 0, 1; t29, t28, t39, t106; t27, t26, t38, t100; t33, t32, t44, t105; 0, 0, 0, 1; t12, -t11, t19, t102; t10, -t9, t18, t97; t15, -t14, t23, t101; 0, 0, 0, 1; t4, -t3, t11, t99; t2, -t1, t9, t96; t6, -t5, t14, t98; 0, 0, 0, 1; t11 * t88 + t4 * t92, t11 * t92 - t4 * t88, t3, t4 * pkin(5) + t3 * pkin(12) + t99; t2 * t92 + t9 * t88, -t2 * t88 + t9 * t92, t1, t2 * pkin(5) + t1 * pkin(12) + t96; t14 * t88 + t6 * t92, t14 * t92 - t6 * t88, t5, t6 * pkin(5) + t5 * pkin(12) + t98; 0, 0, 0, 1;];
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
