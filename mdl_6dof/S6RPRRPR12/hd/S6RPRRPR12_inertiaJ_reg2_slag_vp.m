% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t93 = sin(qJ(4));
t82 = t93 ^ 2;
t96 = cos(qJ(4));
t84 = t96 ^ 2;
t162 = t82 + t84;
t87 = sin(pkin(7));
t94 = sin(qJ(3));
t142 = t87 * t94;
t90 = cos(pkin(7));
t56 = t96 * t142 + t93 * t90;
t51 = t56 ^ 2;
t54 = t93 * t142 - t96 * t90;
t79 = t87 ^ 2;
t97 = cos(qJ(3));
t72 = t79 * t97 ^ 2;
t161 = t54 ^ 2 + t51 + t72;
t86 = sin(pkin(12));
t88 = sin(pkin(6));
t143 = t86 * t88;
t91 = cos(pkin(6));
t152 = pkin(1) * t91;
t89 = cos(pkin(12));
t69 = t89 * t152;
t40 = t91 * pkin(2) + t69 + (-pkin(9) * t90 - qJ(2)) * t143;
t47 = (-pkin(9) * t86 * t87 - pkin(2) * t89 - pkin(1)) * t88;
t108 = t40 * t90 + t47 * t87;
t139 = t89 * t90;
t117 = t88 * t139;
t128 = qJ(2) * t88;
t53 = t89 * t128 + t86 * t152;
t34 = (t87 * t91 + t117) * pkin(9) + t53;
t17 = t108 * t97 - t94 * t34;
t39 = t91 * t142 + (t94 * t139 + t86 * t97) * t88;
t140 = t88 * t89;
t50 = -t87 * t140 + t91 * t90;
t29 = t39 * t93 - t50 * t96;
t160 = t29 ^ 2;
t31 = t39 * t96 + t50 * t93;
t28 = t31 ^ 2;
t141 = t87 * t97;
t37 = -t97 * t117 - t91 * t141 + t94 * t143;
t36 = t37 ^ 2;
t159 = -0.2e1 * t31;
t158 = -0.2e1 * t37;
t157 = 0.2e1 * t88;
t156 = -0.2e1 * t93;
t155 = 0.2e1 * t96;
t154 = 2 * qJ(5);
t153 = pkin(4) + pkin(11);
t151 = t37 * pkin(4);
t92 = sin(qJ(6));
t95 = cos(qJ(6));
t19 = -t29 * t95 + t37 * t92;
t150 = t19 * t92;
t21 = t29 * t92 + t37 * t95;
t149 = t21 * t95;
t148 = t29 * t96;
t147 = t31 * t37;
t146 = t31 * t92;
t26 = t31 * t93;
t145 = t37 * t29;
t80 = t88 ^ 2;
t144 = t80 * t89;
t138 = t92 * t93;
t137 = t92 * t96;
t136 = t92 * t153;
t135 = t93 * t37;
t134 = t93 * t96;
t133 = t95 * t92;
t132 = t95 * t96;
t131 = t95 * t153;
t130 = t96 * t37;
t24 = -t87 * t40 + t90 * t47;
t12 = t37 * pkin(3) - t39 * pkin(10) + t24;
t18 = t108 * t94 + t97 * t34;
t16 = t50 * pkin(10) + t18;
t9 = t93 * t12 + t96 * t16;
t129 = t162 * pkin(10) ^ 2;
t81 = t92 ^ 2;
t83 = t95 ^ 2;
t70 = t81 + t83;
t127 = qJ(5) * t96;
t35 = t37 * qJ(5);
t126 = t56 * qJ(5);
t125 = t29 * t159;
t124 = t91 * t157;
t123 = -0.2e1 * t134;
t121 = pkin(10) * t135;
t120 = pkin(10) * t130;
t119 = t93 * t141;
t118 = t96 * t141;
t116 = t92 * t132;
t6 = -t35 - t9;
t115 = -t93 * qJ(5) - pkin(3);
t8 = t96 * t12 - t93 * t16;
t3 = t31 * pkin(5) - t153 * t37 - t8;
t15 = -t50 * pkin(3) - t17;
t101 = -t31 * qJ(5) + t15;
t5 = t153 * t29 + t101;
t1 = t95 * t3 - t92 * t5;
t2 = t92 * t3 + t95 * t5;
t114 = t1 * t95 + t2 * t92;
t7 = -t8 - t151;
t113 = -t6 * t96 + t7 * t93;
t112 = -t8 * t93 + t9 * t96;
t111 = -pkin(4) * t93 + t127;
t110 = -t56 * t29 + t54 * t31;
t109 = -t93 * t29 + t31 * t96;
t58 = -t153 * t96 + t115;
t76 = t93 * pkin(10);
t65 = t93 * pkin(5) + t76;
t41 = -t92 * t58 + t95 * t65;
t42 = t95 * t58 + t92 * t65;
t22 = t41 * t95 + t42 * t92;
t43 = t92 * t141 + t95 * t54;
t44 = t95 * t141 - t92 * t54;
t23 = t43 * t95 - t44 * t92;
t107 = t54 * t93 + t56 * t96;
t106 = -t153 * t93 + t127;
t105 = t29 * t141 + t54 * t37;
t104 = t31 * t141 + t56 * t37;
t103 = (t26 - t148) * pkin(10);
t102 = t107 * pkin(10);
t99 = qJ(5) ^ 2;
t78 = t96 * pkin(10);
t74 = t95 * t93;
t71 = 0.2e1 * t134;
t66 = t96 * pkin(5) + t78;
t62 = -t96 * pkin(4) + t115;
t61 = 0.2e1 * t162 * pkin(10);
t59 = t70 * t153;
t52 = -t86 * t128 + t69;
t27 = t31 * t95;
t10 = t29 * pkin(4) + t101;
t4 = -t29 * pkin(5) - t6;
t11 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t80 * t86 ^ 2, 0.2e1 * t86 * t144, t86 * t124, t80 * t89 ^ 2, t89 * t124, t91 ^ 2, 0.2e1 * pkin(1) * t144 + 0.2e1 * t52 * t91, -0.2e1 * t80 * pkin(1) * t86 - 0.2e1 * t53 * t91 (-t52 * t86 + t53 * t89) * t157, t80 * pkin(1) ^ 2 + t52 ^ 2 + t53 ^ 2, t39 ^ 2, t39 * t158, 0.2e1 * t39 * t50, t36, t50 * t158, t50 ^ 2, 0.2e1 * t17 * t50 + 0.2e1 * t24 * t37, -0.2e1 * t18 * t50 + 0.2e1 * t24 * t39, -0.2e1 * t17 * t39 - 0.2e1 * t18 * t37, t17 ^ 2 + t18 ^ 2 + t24 ^ 2, t28, t125, 0.2e1 * t147, t160, -0.2e1 * t145, t36, 0.2e1 * t15 * t29 + 0.2e1 * t8 * t37, 0.2e1 * t15 * t31 - 0.2e1 * t9 * t37, -0.2e1 * t9 * t29 - 0.2e1 * t8 * t31, t15 ^ 2 + t8 ^ 2 + t9 ^ 2, t36, -0.2e1 * t147, 0.2e1 * t145, t28, t125, t160, 0.2e1 * t6 * t29 + 0.2e1 * t7 * t31, -0.2e1 * t10 * t29 + 0.2e1 * t7 * t37, -0.2e1 * t10 * t31 - 0.2e1 * t6 * t37, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, 0.2e1 * t21 * t31, t19 ^ 2, t19 * t159, t28, 0.2e1 * t1 * t31 + 0.2e1 * t4 * t19, -0.2e1 * t2 * t31 + 0.2e1 * t4 * t21, -0.2e1 * t1 * t21 - 0.2e1 * t2 * t19, t1 ^ 2 + t2 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t143, 0, -t88 * pkin(1), 0, 0, 0, 0, 0, 0, t50 * t141 + t90 * t37, -t50 * t142 + t90 * t39 (-t37 * t94 - t39 * t97) * t87, t24 * t90 + (t17 * t97 + t18 * t94) * t87, 0, 0, 0, 0, 0, 0, -t105, -t104, t110, -t15 * t141 - t8 * t54 + t9 * t56, 0, 0, 0, 0, 0, 0, t110, t105, t104, -t10 * t141 + t7 * t54 - t6 * t56, 0, 0, 0, 0, 0, 0, t56 * t19 + t43 * t31, t56 * t21 + t44 * t31, t44 * t19 - t43 * t21, t1 * t43 - t2 * t44 + t4 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t94 ^ 2 + t90 ^ 2 + t72, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43 ^ 2 + t44 ^ 2 + t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, -t37, t50, t17, -t18, 0, 0, t26, t109, t135, -t148, t130, 0, -pkin(3) * t29 - t15 * t96 - t121, -pkin(3) * t31 + t15 * t93 - t120, t103 + t112, -t15 * pkin(3) + pkin(10) * t112, 0, -t135, -t130, t26, t109, -t148, t103 + t113, t10 * t96 - t62 * t29 + t121, -t10 * t93 - t62 * t31 + t120, pkin(10) * t113 + t10 * t62, -t21 * t137 (-t149 + t150) * t96, -t137 * t31 + t21 * t93, t19 * t132, -t132 * t31 - t19 * t93, t26, t1 * t93 + t132 * t4 + t66 * t19 + t41 * t31, -t137 * t4 - t2 * t93 + t66 * t21 - t42 * t31, -t42 * t19 - t41 * t21 + (t1 * t92 - t2 * t95) * t96, t1 * t41 + t2 * t42 + t4 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t142, 0, 0, 0, 0, 0, 0, 0, 0, t118, -t119, t107, pkin(3) * t141 + t102, 0, 0, 0, 0, 0, 0, t107, -t118, t119, -t62 * t141 + t102, 0, 0, 0, 0, 0, 0, t132 * t56 + t43 * t93, -t137 * t56 + t44 * t93 (t43 * t92 + t44 * t95) * t96, t43 * t41 - t44 * t42 + t56 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t82, t71, 0, t84, 0, 0, pkin(3) * t155, pkin(3) * t156, t61, pkin(3) ^ 2 + t129, 0, 0, 0, t82, t71, t84, t61, t62 * t155, t62 * t156, t62 ^ 2 + t129, t81 * t84, 0.2e1 * t84 * t133, t92 * t123, t83 * t84, t95 * t123, t82, 0.2e1 * t132 * t66 + 0.2e1 * t41 * t93, -0.2e1 * t137 * t66 - 0.2e1 * t42 * t93 (t41 * t92 - t42 * t95) * t155, t41 ^ 2 + t42 ^ 2 + t66 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, t37, t8, -t9, 0, 0, t37, -t31, t29, 0, 0, 0, -t31 * pkin(4) - qJ(5) * t29, -t8 - 0.2e1 * t151, -t6 + t35, -t7 * pkin(4) - t6 * qJ(5), t149, -t95 * t19 - t21 * t92, t27, t150, -t146, 0, qJ(5) * t19 - t131 * t31 + t4 * t92, qJ(5) * t21 + t136 * t31 + t4 * t95 (t153 * t21 - t1) * t95 + (t153 * t19 - t2) * t92, t4 * qJ(5) - t114 * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, -t56, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t56, -t54 * pkin(4) + t126, 0, 0, 0, 0, 0, 0, t56 * t92, t56 * t95, -t23, -t153 * t23 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, t96, 0, -t76, -t78, 0, 0, 0, -t93, -t96, 0, 0, 0, t111, t76, t78, t111 * pkin(10), -t116 (t81 - t83) * t96, t74, t116, -t138, 0, t106 * t95 + t66 * t92, -t106 * t92 + t66 * t95, -t22, t66 * qJ(5) - t153 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(4), t154, pkin(4) ^ 2 + t99, t83, -0.2e1 * t133, 0, t81, 0, 0, t92 * t154, t95 * t154, 0.2e1 * t59, t153 ^ 2 * t70 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t37, 0, t7, 0, 0, 0, 0, 0, 0, t27, -t146, -t149 - t150, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, 0, t76, 0, 0, 0, 0, 0, 0, t74, -t138, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, t31, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t44, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t137, 0, -t132, t93, t41, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, 0, -t92, 0, -t131, t136, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t95, -t92, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t11;
