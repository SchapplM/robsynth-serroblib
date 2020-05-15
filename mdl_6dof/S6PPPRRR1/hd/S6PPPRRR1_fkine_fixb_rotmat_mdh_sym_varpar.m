% Calculate forward kinematics (homogenous transformation matrices) for fixed-base
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2018-11-23 14:48
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function [T_c_mdh, Tc_stack] = S6PPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_fkine_fixb_rotmat_mdh_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From fkine_mdh_floatb_twist_rotmat_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 14:48:37
% EndTime: 2018-11-23 14:48:38
% DurationCPUTime: 0.46s
% Computational Cost: add. (3259->104), mult. (3274->133), div. (0->0), fcn. (3324->30), ass. (0->96)
t74 = pkin(7) + pkin(14);
t63 = sin(t74) / 0.2e1;
t75 = pkin(7) - pkin(14);
t67 = sin(t75);
t50 = t63 + t67 / 0.2e1;
t76 = pkin(6) + pkin(13);
t64 = sin(t76) / 0.2e1;
t77 = pkin(6) - pkin(13);
t68 = sin(t77);
t52 = t64 + t68 / 0.2e1;
t65 = cos(t75) / 0.2e1;
t69 = cos(t74);
t54 = t65 + t69 / 0.2e1;
t66 = cos(t77) / 0.2e1;
t70 = cos(t76);
t57 = t66 - t70 / 0.2e1;
t78 = sin(pkin(14));
t89 = cos(pkin(6));
t32 = t89 * t50 + t52 * t54 - t57 * t78;
t88 = cos(pkin(7));
t120 = t89 * t88;
t82 = sin(pkin(7));
t48 = -t52 * t82 + t120;
t81 = sin(pkin(8));
t87 = cos(pkin(8));
t23 = -t32 * t81 + t48 * t87;
t80 = sin(pkin(12));
t83 = sin(pkin(6));
t122 = t80 * t83;
t56 = t66 + t70 / 0.2e1;
t79 = sin(pkin(13));
t86 = cos(pkin(12));
t46 = -t80 * t56 - t86 * t79;
t53 = t64 - t68 / 0.2e1;
t85 = cos(pkin(13));
t47 = -t80 * t53 + t86 * t85;
t28 = t122 * t50 + t46 * t54 - t47 * t78;
t115 = t88 * t122;
t39 = -t46 * t82 + t115;
t19 = -t28 * t81 + t39 * t87;
t121 = t86 * t83;
t44 = t86 * t56 - t80 * t79;
t45 = t86 * t53 + t80 * t85;
t26 = -t121 * t50 + t44 * t54 - t45 * t78;
t38 = -t121 * t88 - t44 * t82;
t18 = -t26 * t81 + t38 * t87;
t119 = qJ(2) * t83;
t118 = qJ(3) * t82;
t117 = pkin(8) - qJ(4);
t116 = pkin(8) + qJ(4);
t114 = qJ(1) + 0;
t113 = t86 * pkin(1) + t80 * t119 + 0;
t112 = t89 * qJ(2) + t114;
t111 = cos(t116);
t110 = sin(t117);
t109 = cos(t117) / 0.2e1;
t108 = sin(t116) / 0.2e1;
t107 = t80 * pkin(1) - t119 * t86 + 0;
t106 = t47 * pkin(2) + qJ(3) * t115 - t118 * t46 + t113;
t105 = t57 * pkin(2) + qJ(3) * t120 - t118 * t52 + t112;
t104 = t109 + t111 / 0.2e1;
t103 = t108 + t110 / 0.2e1;
t51 = t63 - t67 / 0.2e1;
t55 = t65 - t69 / 0.2e1;
t84 = cos(pkin(14));
t29 = t122 * t55 + t46 * t51 + t47 * t84;
t102 = t29 * pkin(3) + t19 * pkin(9) + t106;
t33 = t52 * t51 + t89 * t55 + t57 * t84;
t101 = t33 * pkin(3) + t23 * pkin(9) + t105;
t100 = t45 * pkin(2) + qJ(3) * t38 + t107;
t92 = sin(qJ(4));
t11 = -t103 * t39 - t104 * t28 + t29 * t92;
t58 = t108 - t110 / 0.2e1;
t59 = t109 - t111 / 0.2e1;
t95 = cos(qJ(4));
t12 = t28 * t58 + t29 * t95 + t39 * t59;
t99 = t12 * pkin(4) + t11 * pkin(10) + t102;
t14 = -t103 * t48 - t104 * t32 + t33 * t92;
t15 = t32 * t58 + t33 * t95 + t48 * t59;
t98 = t15 * pkin(4) + t14 * pkin(10) + t101;
t27 = -t121 * t55 + t44 * t51 + t45 * t84;
t97 = t27 * pkin(3) + t18 * pkin(9) + t100;
t10 = t26 * t58 + t27 * t95 + t38 * t59;
t9 = -t103 * t38 - t104 * t26 + t27 * t92;
t96 = t10 * pkin(4) + t9 * pkin(10) + t97;
t94 = cos(qJ(5));
t93 = cos(qJ(6));
t91 = sin(qJ(5));
t90 = sin(qJ(6));
t6 = t15 * t94 + t23 * t91;
t5 = t15 * t91 - t23 * t94;
t4 = t12 * t94 + t19 * t91;
t3 = t12 * t91 - t19 * t94;
t2 = t10 * t94 + t18 * t91;
t1 = t10 * t91 - t18 * t94;
t7 = [1, 0, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0; 0, 0, 0, 1; t86, -t80, 0, 0; t80, t86, 0, 0; 0, 0, 1, t114; 0, 0, 0, 1; t47, t46, t122, t113; t45, t44, -t121, t107; t57, t52, t89, t112; 0, 0, 0, 1; t29, t28, t39, t106; t27, t26, t38, t100; t33, t32, t48, t105; 0, 0, 0, 1; t12, -t11, t19, t102; t10, -t9, t18, t97; t15, -t14, t23, t101; 0, 0, 0, 1; t4, -t3, t11, t99; t2, -t1, t9, t96; t6, -t5, t14, t98; 0, 0, 0, 1; t11 * t90 + t4 * t93, t11 * t93 - t4 * t90, t3, t4 * pkin(5) + t3 * pkin(11) + t99; t2 * t93 + t9 * t90, -t2 * t90 + t9 * t93, t1, t2 * pkin(5) + t1 * pkin(11) + t96; t14 * t90 + t6 * t93, t14 * t93 - t6 * t90, t5, t6 * pkin(5) + t5 * pkin(11) + t98; 0, 0, 0, 1;];
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
