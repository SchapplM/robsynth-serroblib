% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t105 = sin(pkin(11));
t108 = cos(pkin(11));
t111 = sin(qJ(3));
t113 = cos(qJ(3));
t83 = t105 * t113 + t108 * t111;
t172 = -0.2e1 * t83;
t104 = sin(pkin(12));
t107 = cos(pkin(12));
t110 = sin(qJ(6));
t161 = cos(qJ(6));
t171 = -t110 * t104 + t161 * t107;
t109 = cos(pkin(6));
t106 = sin(pkin(6));
t112 = sin(qJ(2));
t136 = t106 * t112;
t69 = t109 * t113 - t111 * t136;
t70 = t109 * t111 + t113 * t136;
t48 = t105 * t70 - t108 * t69;
t44 = t48 ^ 2;
t152 = -qJ(4) - pkin(9);
t122 = t152 * t111;
t87 = t152 * t113;
t59 = -t105 * t87 - t108 * t122;
t170 = t59 ^ 2;
t79 = t105 * t111 - t108 * t113;
t75 = t79 ^ 2;
t169 = -0.2e1 * t48;
t168 = 0.2e1 * t48;
t167 = 0.2e1 * t79;
t157 = t108 * pkin(3);
t97 = -pkin(4) - t157;
t86 = -t107 * pkin(5) + t97;
t166 = 0.2e1 * t86;
t98 = -t113 * pkin(3) - pkin(2);
t165 = 0.2e1 * t98;
t164 = 0.2e1 * t106;
t163 = 0.2e1 * t113;
t158 = t105 * pkin(3);
t93 = qJ(5) + t158;
t162 = pkin(10) + t93;
t160 = pkin(1) * t112;
t114 = cos(qJ(2));
t159 = pkin(1) * t114;
t135 = t106 * t114;
t50 = t105 * t69 + t108 * t70;
t37 = t104 * t50 + t107 * t135;
t39 = -t104 * t135 + t107 * t50;
t21 = -t110 * t37 + t161 * t39;
t156 = t21 * t171;
t47 = t171 * t83;
t155 = t47 * t171;
t34 = t48 * t79;
t84 = t161 * t104 + t110 * t107;
t154 = t84 * t48;
t153 = t84 * t79;
t128 = pkin(3) * t135;
t127 = pkin(8) * t135;
t64 = t127 + (pkin(9) + t160) * t109;
t65 = (-pkin(2) * t114 - pkin(9) * t112 - pkin(1)) * t106;
t40 = -t111 * t64 + t113 * t65;
t26 = -qJ(4) * t70 - t128 + t40;
t41 = t111 * t65 + t113 * t64;
t31 = qJ(4) * t69 + t41;
t14 = t105 * t26 + t108 * t31;
t11 = -qJ(5) * t135 + t14;
t90 = pkin(8) * t136;
t63 = t90 + (-pkin(2) - t159) * t109;
t55 = -t69 * pkin(3) + t63;
t18 = t48 * pkin(4) - t50 * qJ(5) + t55;
t6 = t104 * t18 + t107 * t11;
t151 = t105 * t31 - t108 * t26;
t56 = t79 * pkin(4) - t83 * qJ(5) + t98;
t61 = t105 * t122 - t108 * t87;
t28 = t104 * t56 + t107 * t61;
t150 = t104 * t48;
t149 = t104 * t79;
t148 = t104 * t83;
t147 = t107 * t39;
t43 = t107 * t48;
t146 = t107 * t83;
t145 = t37 * t107;
t144 = t39 * t104;
t143 = t59 * t104;
t142 = t59 * t107;
t141 = t69 * t113;
t140 = t70 * t111;
t101 = t107 ^ 2;
t99 = t104 ^ 2;
t139 = t99 + t101;
t100 = t106 ^ 2;
t138 = t100 * t114;
t137 = t104 * t107;
t134 = t109 * t112;
t102 = t111 ^ 2;
t103 = t113 ^ 2;
t132 = t102 + t103;
t131 = t79 * t172;
t130 = -0.2e1 * t135;
t129 = 0.2e1 * t135;
t12 = pkin(4) * t135 + t151;
t126 = t83 * t137;
t125 = t111 * t135;
t124 = t113 * t135;
t5 = -t104 * t11 + t107 * t18;
t27 = -t104 * t61 + t107 * t56;
t121 = t104 * t6 + t107 * t5;
t120 = -t79 * t93 + t83 * t97;
t119 = t104 * t28 + t107 * t27;
t118 = -t104 * t27 + t107 * t28;
t117 = -t40 * t111 + t41 * t113;
t94 = t100 * t114 ^ 2;
t78 = t84 ^ 2;
t77 = t83 ^ 2;
t76 = t171 ^ 2;
t74 = pkin(1) * t134 + t127;
t73 = t109 * t159 - t90;
t72 = t162 * t107;
t71 = t162 * t104;
t68 = t107 * t79;
t58 = t171 * t79;
t54 = -t110 * t71 + t161 * t72;
t53 = -t110 * t72 - t161 * t71;
t45 = t84 * t83;
t42 = pkin(5) * t148 + t59;
t36 = t84 * t45;
t35 = t171 * t48;
t32 = t104 * t37;
t23 = -pkin(10) * t148 + t28;
t22 = pkin(5) * t79 - pkin(10) * t146 + t27;
t19 = t110 * t39 + t161 * t37;
t15 = t84 * t19;
t9 = t110 * t22 + t161 * t23;
t8 = -t110 * t23 + t161 * t22;
t7 = pkin(5) * t37 + t12;
t4 = -pkin(10) * t37 + t6;
t3 = pkin(5) * t48 - pkin(10) * t39 + t5;
t2 = t110 * t3 + t161 * t4;
t1 = -t110 * t4 + t161 * t3;
t10 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t100 * t112 ^ 2, 0.2e1 * t112 * t138, t134 * t164, t94, t109 * t129, t109 ^ 2, 0.2e1 * pkin(1) * t138 + 0.2e1 * t109 * t73, -0.2e1 * t100 * t160 - 0.2e1 * t109 * t74 (-t112 * t73 + t114 * t74) * t164, pkin(1) ^ 2 * t100 + t73 ^ 2 + t74 ^ 2, t70 ^ 2, 0.2e1 * t70 * t69, t70 * t130, t69 ^ 2, t69 * t130, t94, -0.2e1 * t40 * t135 - 0.2e1 * t63 * t69, 0.2e1 * t41 * t135 + 0.2e1 * t63 * t70, -0.2e1 * t40 * t70 + 0.2e1 * t41 * t69, t40 ^ 2 + t41 ^ 2 + t63 ^ 2, t50 ^ 2, t50 * t169, t50 * t130, t44, t48 * t129, t94, 0.2e1 * t135 * t151 + 0.2e1 * t48 * t55, 0.2e1 * t14 * t135 + 0.2e1 * t50 * t55, -0.2e1 * t14 * t48 + 0.2e1 * t151 * t50, t14 ^ 2 + t151 ^ 2 + t55 ^ 2, t39 ^ 2, -0.2e1 * t39 * t37, t39 * t168, t37 ^ 2, t37 * t169, t44, 0.2e1 * t12 * t37 + 0.2e1 * t48 * t5, 0.2e1 * t12 * t39 - 0.2e1 * t48 * t6, -0.2e1 * t37 * t6 - 0.2e1 * t39 * t5, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, t21 * t168, t19 ^ 2, t19 * t169, t44, 0.2e1 * t1 * t48 + 0.2e1 * t19 * t7, -0.2e1 * t2 * t48 + 0.2e1 * t21 * t7, -0.2e1 * t1 * t21 - 0.2e1 * t19 * t2, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, t135, t109, t73, -t74, 0, 0, t140, t111 * t69 + t113 * t70, -t125, t141, -t124, 0, pkin(2) * t69 + pkin(9) * t125 - t113 * t63, -pkin(2) * t70 + pkin(9) * t124 + t111 * t63 (t140 + t141) * pkin(9) + t117, -t63 * pkin(2) + t117 * pkin(9), t50 * t83, -t48 * t83 - t50 * t79, -t83 * t135, t34, t79 * t135, 0, t59 * t135 + t48 * t98 + t55 * t79, t61 * t135 + t50 * t98 + t55 * t83, -t14 * t79 + t151 * t83 - t48 * t61 + t50 * t59, t14 * t61 + t151 * t59 + t55 * t98, t39 * t146 (-t144 - t145) * t83, t48 * t146 + t39 * t79, t37 * t148, -t48 * t148 - t37 * t79, t34, t12 * t148 + t27 * t48 + t37 * t59 + t5 * t79, t12 * t146 - t28 * t48 + t39 * t59 - t6 * t79, -t121 * t83 - t27 * t39 - t28 * t37, t12 * t59 + t27 * t5 + t28 * t6, t21 * t47, -t19 * t47 - t21 * t45, t21 * t79 + t47 * t48, t19 * t45, -t19 * t79 - t45 * t48, t34, t1 * t79 + t19 * t42 + t45 * t7 + t48 * t8, -t2 * t79 + t21 * t42 + t47 * t7 - t48 * t9, -t1 * t47 - t19 * t9 - t2 * t45 - t21 * t8, t1 * t8 + t2 * t9 + t42 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t102, t111 * t163, 0, t103, 0, 0, pkin(2) * t163, -0.2e1 * pkin(2) * t111, 0.2e1 * t132 * pkin(9), t132 * pkin(9) ^ 2 + pkin(2) ^ 2, t77, t131, 0, t75, 0, 0, t79 * t165, t83 * t165, 0.2e1 * t59 * t83 - 0.2e1 * t61 * t79, t61 ^ 2 + t98 ^ 2 + t170, t101 * t77, -0.2e1 * t77 * t137, t146 * t167, t99 * t77, t104 * t131, t75, 0.2e1 * t83 * t143 + 0.2e1 * t27 * t79, 0.2e1 * t83 * t142 - 0.2e1 * t28 * t79, t119 * t172, t27 ^ 2 + t28 ^ 2 + t170, t47 ^ 2, -0.2e1 * t47 * t45, t47 * t167, t45 ^ 2, -t45 * t167, t75, 0.2e1 * t42 * t45 + 0.2e1 * t79 * t8, 0.2e1 * t42 * t47 - 0.2e1 * t79 * t9, -0.2e1 * t45 * t9 - 0.2e1 * t47 * t8, t42 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, t69, -t135, t40, -t41, 0, 0, 0, 0, t50, 0, -t48, -t135, -t108 * t128 - t151, t105 * t128 - t14 (-t105 * t48 - t108 * t50) * pkin(3) (t105 * t14 - t108 * t151) * pkin(3), t144, -t32 + t147, t150, -t145, t43, 0, -t107 * t12 - t93 * t150 + t37 * t97, t104 * t12 + t39 * t97 - t93 * t43 (-t37 * t93 + t6) * t107 + (t39 * t93 - t5) * t104, t12 * t97 + (-t5 * t104 + t6 * t107) * t93, t21 * t84, -t15 + t156, t154, -t19 * t171, t35, 0, -t171 * t7 + t19 * t86 + t48 * t53, t21 * t86 - t48 * t54 + t7 * t84, -t1 * t84 + t171 * t2 - t19 * t54 - t21 * t53, t1 * t53 + t2 * t54 + t7 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, t113, 0, -t111 * pkin(9), -t113 * pkin(9), 0, 0, 0, 0, t83, 0, -t79, 0, -t59, -t61 (-t105 * t79 - t108 * t83) * pkin(3) (t105 * t61 - t108 * t59) * pkin(3), t126 (t101 - t99) * t83, t149, -t126, t68, 0, t104 * t120 - t142, t120 * t107 + t143, t118, t118 * t93 + t59 * t97, t47 * t84, -t36 + t155, t153, -t45 * t171, t58, 0, -t171 * t42 + t45 * t86 + t53 * t79, t42 * t84 + t47 * t86 - t54 * t79, t171 * t9 - t45 * t54 - t47 * t53 - t8 * t84, t42 * t86 + t53 * t8 + t54 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t157, -0.2e1 * t158, 0 (t105 ^ 2 + t108 ^ 2) * pkin(3) ^ 2, t99, 0.2e1 * t137, 0, t101, 0, 0, -0.2e1 * t97 * t107, 0.2e1 * t97 * t104, 0.2e1 * t139 * t93, t139 * t93 ^ 2 + t97 ^ 2, t78, 0.2e1 * t84 * t171, 0, t76, 0, 0, -t171 * t166, t84 * t166, 0.2e1 * t171 * t54 - 0.2e1 * t53 * t84, t53 ^ 2 + t54 ^ 2 + t86 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t50, 0, t55, 0, 0, 0, 0, 0, 0, t43, -t150, -t32 - t147, t121, 0, 0, 0, 0, 0, 0, t35, -t154, -t15 - t156, t1 * t171 + t2 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t83, 0, t98, 0, 0, 0, 0, 0, 0, t68, -t149, -t139 * t83, t119, 0, 0, 0, 0, 0, 0, t58, -t153, -t36 - t155, t171 * t8 + t84 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171 * t53 + t54 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 + t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, t39, 0, t12, 0, 0, 0, 0, 0, 0, t19, t21, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t148, t146, 0, t59, 0, 0, 0, 0, 0, 0, t45, t47, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, t104, 0, t97, 0, 0, 0, 0, 0, 0, -t171, t84, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, t48, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, -t45, t79, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, t171, 0, t53, -t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, -t84, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t10;
