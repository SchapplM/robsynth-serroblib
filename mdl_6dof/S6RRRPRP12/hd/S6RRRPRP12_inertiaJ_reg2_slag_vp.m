% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP12_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP12_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP12_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t93 = sin(qJ(5));
t86 = t93 ^ 2;
t96 = cos(qJ(5));
t88 = t96 ^ 2;
t67 = t86 + t88;
t94 = sin(qJ(3));
t87 = t94 ^ 2;
t97 = cos(qJ(3));
t89 = t97 ^ 2;
t167 = t87 + t89;
t91 = sin(pkin(6));
t98 = cos(qJ(2));
t144 = t91 * t98;
t95 = sin(qJ(2));
t145 = t91 * t95;
t92 = cos(pkin(6));
t46 = t94 * t145 - t92 * t97;
t29 = t93 * t144 + t46 * t96;
t24 = t29 * t93;
t31 = t96 * t144 - t46 * t93;
t25 = t96 * t31;
t166 = t24 + t25;
t165 = (-t25 + t24) * t97;
t164 = t29 ^ 2;
t163 = t46 ^ 2;
t48 = t97 * t145 + t92 * t94;
t45 = t48 ^ 2;
t162 = -0.2e1 * t48;
t161 = 0.2e1 * t91;
t160 = -0.2e1 * t94;
t159 = 0.2e1 * t97;
t158 = 2 * qJ(4);
t157 = pkin(3) + pkin(10);
t156 = pkin(1) * t95;
t155 = pkin(1) * t98;
t134 = qJ(6) * t48;
t68 = pkin(8) * t145;
t37 = t68 + (-pkin(2) - t155) * t92;
t102 = -t48 * qJ(4) + t37;
t11 = t157 * t46 + t102;
t130 = pkin(8) * t144;
t38 = t130 + (pkin(9) + t156) * t92;
t39 = (-pkin(2) * t98 - pkin(9) * t95 - pkin(1)) * t91;
t18 = -t94 * t38 + t97 * t39;
t69 = pkin(3) * t144;
t17 = -t18 + t69;
t8 = t48 * pkin(4) + pkin(10) * t144 + t17;
t4 = t96 * t11 + t93 * t8;
t1 = t134 + t4;
t154 = t1 * t93;
t153 = t4 * t93;
t152 = t48 * pkin(5);
t151 = t94 * pkin(5);
t150 = t29 * t48;
t149 = t31 * t29;
t148 = t31 * t157;
t147 = t46 * t97;
t41 = t48 * t93;
t42 = t48 * t94;
t43 = t48 * t96;
t85 = t91 ^ 2;
t146 = t85 * t98;
t143 = t93 * t97;
t142 = t93 * t157;
t141 = t94 * t97;
t140 = t94 * t157;
t139 = t96 * t93;
t138 = t96 * t97;
t79 = t96 * t157;
t19 = t97 * t38 + t94 * t39;
t117 = -t94 * qJ(4) - pkin(2);
t53 = -t157 * t97 + t117;
t81 = t94 * pkin(9);
t63 = t94 * pkin(4) + t81;
t27 = t96 * t53 + t93 * t63;
t137 = t67 * t157 ^ 2;
t136 = t167 * pkin(9) ^ 2;
t83 = t97 * pkin(9);
t64 = t97 * pkin(4) + t83;
t135 = qJ(4) * t97;
t133 = t94 * qJ(6);
t132 = t46 * t162;
t131 = t92 * t161;
t129 = t29 * t138;
t128 = t31 * t143;
t127 = t48 * t144;
t126 = t46 * t144;
t125 = t94 * t144;
t124 = t97 * t144;
t123 = t89 * t139;
t122 = t48 * t142;
t121 = t94 * t138;
t120 = t48 * t79;
t119 = qJ(4) * t144;
t118 = t93 * t11 - t96 * t8;
t116 = t93 * t53 - t96 * t63;
t115 = pkin(9) * t125;
t114 = pkin(9) * t124;
t2 = t118 - t152;
t113 = -t2 * t96 + t154;
t112 = -t118 * t96 + t153;
t111 = -pkin(3) * t94 + t135;
t61 = pkin(5) * t96 + t93 * qJ(6);
t110 = t93 * pkin(5) - t96 * qJ(6);
t16 = t119 - t19;
t109 = -t16 * t97 + t17 * t94;
t108 = -t18 * t94 + t19 * t97;
t21 = t133 + t27;
t22 = t116 - t151;
t9 = t21 * t93 - t22 * t96;
t13 = -t116 * t96 + t27 * t93;
t107 = t96 * t29 + t31 * t93;
t105 = -t94 * t46 + t48 * t97;
t104 = t48 * t138 - t29 * t94;
t103 = (t42 - t147) * pkin(9);
t12 = -t46 * pkin(4) - t16;
t100 = qJ(4) ^ 2;
t76 = t96 * t94;
t75 = t88 * t89;
t74 = t93 * t94;
t73 = t86 * t89;
t71 = t85 * t98 ^ 2;
t70 = 0.2e1 * t141;
t66 = t94 * t79;
t65 = t93 * t138;
t62 = -0.2e1 * t93 * t141;
t60 = -t97 * pkin(3) + t117;
t59 = qJ(4) + t110;
t58 = 0.2e1 * t167 * pkin(9);
t57 = t67 * t157;
t54 = (-t86 + t88) * t97;
t52 = -0.2e1 * t57;
t51 = t92 * t156 + t130;
t50 = t92 * t155 - t68;
t33 = t61 * t97 + t64;
t28 = t31 ^ 2;
t23 = t29 * t142;
t20 = t31 * t162;
t15 = -t48 * t143 - t31 * t94;
t14 = t46 * pkin(3) + t102;
t5 = -t29 * pkin(5) + t31 * qJ(6) + t12;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t85 * t95 ^ 2, 0.2e1 * t95 * t146, t95 * t131, t71, t98 * t131, t92 ^ 2, 0.2e1 * pkin(1) * t146 + 0.2e1 * t50 * t92, -0.2e1 * t85 * t156 - 0.2e1 * t51 * t92 (-t50 * t95 + t51 * t98) * t161, t85 * pkin(1) ^ 2 + t50 ^ 2 + t51 ^ 2, t45, t132, -0.2e1 * t127, t163, 0.2e1 * t126, t71, -0.2e1 * t18 * t144 + 0.2e1 * t37 * t46, 0.2e1 * t19 * t144 + 0.2e1 * t37 * t48, -0.2e1 * t18 * t48 - 0.2e1 * t19 * t46, t18 ^ 2 + t19 ^ 2 + t37 ^ 2, t71, 0.2e1 * t127, -0.2e1 * t126, t45, t132, t163, 0.2e1 * t16 * t46 + 0.2e1 * t17 * t48, -0.2e1 * t14 * t46 - 0.2e1 * t17 * t144, -0.2e1 * t14 * t48 + 0.2e1 * t16 * t144, t14 ^ 2 + t16 ^ 2 + t17 ^ 2, t28, -0.2e1 * t149, t20, t164, 0.2e1 * t150, t45, -0.2e1 * t118 * t48 - 0.2e1 * t12 * t29, -0.2e1 * t12 * t31 - 0.2e1 * t4 * t48, -0.2e1 * t118 * t31 + 0.2e1 * t29 * t4, t118 ^ 2 + t12 ^ 2 + t4 ^ 2, t28, t20, 0.2e1 * t149, t45, -0.2e1 * t150, t164, -0.2e1 * t2 * t48 - 0.2e1 * t29 * t5, 0.2e1 * t1 * t29 - 0.2e1 * t2 * t31, 0.2e1 * t1 * t48 + 0.2e1 * t31 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t145, 0, t144, t92, t50, -t51, 0, 0, t42, t105, -t125, -t147, -t124, 0, -pkin(2) * t46 - t37 * t97 + t115, -pkin(2) * t48 + t37 * t94 + t114, t103 + t108, -t37 * pkin(2) + pkin(9) * t108, 0, t125, t124, t42, t105, -t147, t103 + t109, t14 * t97 - t60 * t46 - t115, -t14 * t94 - t60 * t48 - t114, pkin(9) * t109 + t14 * t60, t128, -t165, t15, -t129, -t104, t42, -t116 * t48 - t118 * t94 + t12 * t138 - t64 * t29, -t12 * t143 - t27 * t48 - t64 * t31 - t4 * t94, -t116 * t31 + t27 * t29 + (-t118 * t93 - t4 * t96) * t97, t116 * t118 + t12 * t64 + t4 * t27, t128, t15, t165, t42, t104, -t129, t138 * t5 - t2 * t94 - t22 * t48 - t33 * t29, t21 * t29 - t22 * t31 + (-t1 * t96 - t2 * t93) * t97, t1 * t94 + t143 * t5 + t21 * t48 + t33 * t31, t1 * t21 + t2 * t22 + t33 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t87, t70, 0, t89, 0, 0, pkin(2) * t159, pkin(2) * t160, t58, pkin(2) ^ 2 + t136, 0, 0, 0, t87, t70, t89, t58, t60 * t159, t60 * t160, t60 ^ 2 + t136, t73, 0.2e1 * t123, t62, t75, -0.2e1 * t121, t87, -0.2e1 * t116 * t94 + 0.2e1 * t138 * t64, -0.2e1 * t143 * t64 - 0.2e1 * t27 * t94 (-t116 * t93 - t27 * t96) * t159, t116 ^ 2 + t27 ^ 2 + t64 ^ 2, t73, t62, -0.2e1 * t123, t87, 0.2e1 * t121, t75, 0.2e1 * t138 * t33 - 0.2e1 * t22 * t94 (-t21 * t96 - t22 * t93) * t159, 0.2e1 * t143 * t33 + 0.2e1 * t21 * t94, t21 ^ 2 + t22 ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t46, -t144, t18, -t19, 0, 0, -t144, -t48, t46, 0, 0, 0, -pkin(3) * t48 - qJ(4) * t46, -t18 + 0.2e1 * t69, -0.2e1 * t119 + t19, -pkin(3) * t17 - qJ(4) * t16, -t25, t107, t43, -t24, -t41, 0, -qJ(4) * t29 + t12 * t93 - t120, -qJ(4) * t31 + t12 * t96 + t122, -t153 - t23 + (t118 - t148) * t96, t12 * qJ(4) - t112 * t157, -t25, t43, -t107, 0, t41, -t24, -t59 * t29 + t5 * t93 - t120, -t154 - t23 + (t2 - t148) * t96, t59 * t31 - t5 * t96 - t122, -t113 * t157 + t5 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, t97, 0, -t81, -t83, 0, 0, 0, -t94, -t97, 0, 0, 0, t111, t81, t83, t111 * pkin(9), -t65, -t54, t76, t65, -t74, 0, t135 * t96 + t64 * t93 - t66, t64 * t96 + (-t135 + t140) * t93, -t13, t64 * qJ(4) - t13 * t157, -t65, t76, t54, 0, t74, t65, t138 * t59 + t33 * t93 - t66, -t9, -t33 * t96 + (t59 * t97 - t140) * t93, -t157 * t9 + t33 * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3), t158, pkin(3) ^ 2 + t100, t88, -0.2e1 * t139, 0, t86, 0, 0, t93 * t158, t96 * t158, -t52, t100 + t137, t88, 0, 0.2e1 * t139, 0, 0, t86, 0.2e1 * t59 * t93, -t52, -0.2e1 * t59 * t96, t59 ^ 2 + t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t144, 0, t17, 0, 0, 0, 0, 0, 0, t43, -t41, t166, t112, 0, 0, 0, 0, 0, 0, t43, t166, t41, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, 0, 0, t81, 0, 0, 0, 0, 0, 0, t76, -t74, 0, t13, 0, 0, 0, 0, 0, 0, t76, 0, t74, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t67, -t57, 0, 0, 0, 0, 0, 0, 0, -t67, 0, -t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, t29, t48, -t118, -t4, 0, 0, 0, -t31, 0, t48, -t29, 0, -t118 + 0.2e1 * t152, pkin(5) * t31 + qJ(6) * t29, 0.2e1 * t134 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143, 0, -t138, t94, -t116, -t27, 0, 0, 0, -t143, 0, t94, t138, 0, -t116 + 0.2e1 * t151, t110 * t97, 0.2e1 * t133 + t27, -pkin(5) * t22 + qJ(6) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, -t93, 0, -t79, t142, 0, 0, 0, t96, 0, 0, t93, 0, -t79, -t61, -t142, -t61 * t157; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, -t93, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, t93, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t31, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94, -t143, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, 0, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
