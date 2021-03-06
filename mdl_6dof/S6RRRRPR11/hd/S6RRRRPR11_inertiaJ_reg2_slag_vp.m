% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 23:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 23:37:10
% EndTime: 2019-05-07 23:37:26
% DurationCPUTime: 3.36s
% Computational Cost: add. (5617->299), mult. (12826->596), div. (0->0), fcn. (14802->12), ass. (0->135)
t108 = cos(pkin(6));
t111 = sin(qJ(3));
t114 = cos(qJ(3));
t106 = sin(pkin(6));
t112 = sin(qJ(2));
t137 = t106 * t112;
t74 = -t108 * t114 + t111 * t137;
t73 = t74 ^ 2;
t105 = sin(pkin(12));
t107 = cos(pkin(12));
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t79 = t105 * t110 - t107 * t113;
t97 = -t113 * pkin(4) - pkin(3);
t63 = t79 * pkin(5) + t97;
t165 = 0.2e1 * t63;
t164 = -0.2e1 * t74;
t163 = 0.2e1 * t74;
t76 = t108 * t111 + t114 * t137;
t162 = -0.2e1 * t76;
t161 = 0.2e1 * t97;
t160 = 0.2e1 * t106;
t159 = -0.2e1 * t111;
t158 = -0.2e1 * t114;
t157 = 0.2e1 * t114;
t156 = cos(qJ(6));
t155 = pkin(1) * t112;
t115 = cos(qJ(2));
t154 = pkin(1) * t115;
t153 = pkin(3) * t113;
t152 = pkin(9) * t110;
t102 = t111 ^ 2;
t151 = t102 * pkin(9);
t150 = t105 * pkin(4);
t149 = t107 * pkin(4);
t99 = t111 * pkin(9);
t148 = t114 * pkin(4);
t147 = -qJ(5) - pkin(10);
t90 = pkin(8) * t137;
t64 = t90 + (-pkin(2) - t154) * t108;
t34 = t74 * pkin(3) - t76 * pkin(10) + t64;
t136 = t106 * t115;
t125 = pkin(8) * t136;
t65 = t125 + (pkin(9) + t155) * t108;
t66 = (-pkin(2) * t115 - pkin(9) * t112 - pkin(1)) * t106;
t38 = t111 * t66 + t114 * t65;
t36 = -pkin(10) * t136 + t38;
t19 = -t110 * t36 + t113 * t34;
t55 = -t110 * t136 + t76 * t113;
t12 = t74 * pkin(4) - t55 * qJ(5) + t19;
t20 = t110 * t34 + t113 * t36;
t53 = t76 * t110 + t113 * t136;
t15 = -t53 * qJ(5) + t20;
t7 = t105 * t12 + t107 * t15;
t131 = t113 * t111;
t85 = -t114 * pkin(3) - t111 * pkin(10) - pkin(2);
t82 = t113 * t85;
t50 = -qJ(5) * t131 + t82 + (-pkin(4) - t152) * t114;
t130 = t113 * t114;
t126 = pkin(9) * t130;
t56 = t126 + (-qJ(5) * t111 + t85) * t110;
t29 = t105 * t50 + t107 * t56;
t134 = t110 * t111;
t84 = pkin(4) * t134 + t99;
t146 = t110 * t74;
t145 = t113 * t74;
t37 = -t111 * t65 + t114 * t66;
t35 = pkin(3) * t136 - t37;
t144 = t35 * t110;
t143 = t35 * t113;
t142 = t53 * t113;
t141 = t55 * t110;
t140 = t74 * t114;
t139 = t76 * t111;
t100 = t106 ^ 2;
t138 = t100 * t115;
t135 = t108 * t112;
t133 = t110 * t113;
t132 = t110 * t114;
t101 = t110 ^ 2;
t103 = t113 ^ 2;
t129 = t101 + t103;
t128 = 0.2e1 * t136;
t127 = t111 * t157;
t124 = t111 * t136;
t123 = t114 * t136;
t122 = t110 * t131;
t109 = sin(qJ(6));
t32 = -t105 * t53 + t107 * t55;
t6 = -t105 * t15 + t107 * t12;
t4 = t74 * pkin(5) - t32 * pkin(11) + t6;
t30 = t105 * t55 + t107 * t53;
t5 = -t30 * pkin(11) + t7;
t1 = -t109 * t5 + t156 * t4;
t28 = -t105 * t56 + t107 * t50;
t86 = t147 * t110;
t87 = t147 * t113;
t57 = t105 * t87 + t107 * t86;
t69 = -t105 * t134 + t107 * t131;
t22 = -t114 * pkin(5) - t69 * pkin(11) + t28;
t81 = t105 * t113 + t107 * t110;
t67 = t81 * t111;
t25 = -t67 * pkin(11) + t29;
t8 = -t109 * t25 + t156 * t22;
t58 = t105 * t86 - t107 * t87;
t121 = -t19 * t110 + t20 * t113;
t61 = -pkin(9) * t132 + t82;
t62 = t110 * t85 + t126;
t120 = -t61 * t110 + t62 * t113;
t119 = -t37 * t111 + t38 * t114;
t2 = t109 * t4 + t156 * t5;
t9 = t109 * t22 + t156 * t25;
t26 = t53 * pkin(4) + t35;
t117 = pkin(9) ^ 2;
t104 = t114 ^ 2;
t98 = t102 * t117;
t96 = pkin(5) + t149;
t93 = t100 * t115 ^ 2;
t78 = pkin(1) * t135 + t125;
t77 = t108 * t154 - t90;
t71 = t109 * t96 + t150 * t156;
t70 = -t109 * t150 + t156 * t96;
t51 = t67 * pkin(5) + t84;
t48 = -t109 * t79 + t156 * t81;
t46 = t109 * t81 + t156 * t79;
t43 = -t79 * pkin(11) + t58;
t42 = -t81 * pkin(11) + t57;
t41 = -t109 * t67 + t156 * t69;
t39 = t109 * t69 + t156 * t67;
t24 = t109 * t42 + t156 * t43;
t23 = -t109 * t43 + t156 * t42;
t18 = -t109 * t30 + t156 * t32;
t16 = t109 * t32 + t156 * t30;
t13 = t30 * pkin(5) + t26;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t100 * t112 ^ 2, 0.2e1 * t112 * t138, t135 * t160, t93, t108 * t128, t108 ^ 2, 0.2e1 * pkin(1) * t138 + 0.2e1 * t77 * t108, -0.2e1 * t100 * t155 - 0.2e1 * t78 * t108 (-t112 * t77 + t115 * t78) * t160, t100 * pkin(1) ^ 2 + t77 ^ 2 + t78 ^ 2, t76 ^ 2, t74 * t162, t136 * t162, t73, t74 * t128, t93, -0.2e1 * t136 * t37 + 0.2e1 * t64 * t74, 0.2e1 * t136 * t38 + 0.2e1 * t64 * t76, -0.2e1 * t37 * t76 - 0.2e1 * t38 * t74, t37 ^ 2 + t38 ^ 2 + t64 ^ 2, t55 ^ 2, -0.2e1 * t55 * t53, t55 * t163, t53 ^ 2, t53 * t164, t73, 0.2e1 * t19 * t74 + 0.2e1 * t35 * t53, -0.2e1 * t20 * t74 + 0.2e1 * t35 * t55, -0.2e1 * t19 * t55 - 0.2e1 * t20 * t53, t19 ^ 2 + t20 ^ 2 + t35 ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, t32 * t163, t30 ^ 2, t30 * t164, t73, 0.2e1 * t26 * t30 + 0.2e1 * t6 * t74, 0.2e1 * t26 * t32 - 0.2e1 * t7 * t74, -0.2e1 * t7 * t30 - 0.2e1 * t6 * t32, t26 ^ 2 + t6 ^ 2 + t7 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, t18 * t163, t16 ^ 2, t16 * t164, t73, 0.2e1 * t1 * t74 + 0.2e1 * t13 * t16, 0.2e1 * t13 * t18 - 0.2e1 * t2 * t74, -0.2e1 * t1 * t18 - 0.2e1 * t2 * t16, t1 ^ 2 + t13 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, t136, t108, t77, -t78, 0, 0, t139, -t111 * t74 + t76 * t114, -t124, -t140, -t123, 0, -pkin(2) * t74 + pkin(9) * t124 - t64 * t114, -pkin(2) * t76 + pkin(9) * t123 + t64 * t111 (t139 - t140) * pkin(9) + t119, -t64 * pkin(2) + pkin(9) * t119, t55 * t131 (-t141 - t142) * t111, -t55 * t114 + t131 * t74, t53 * t134, t53 * t114 - t134 * t74, -t140, -t19 * t114 + t61 * t74 + (pkin(9) * t53 + t144) * t111, t20 * t114 - t62 * t74 + (pkin(9) * t55 + t143) * t111, -t62 * t53 - t61 * t55 + (-t110 * t20 - t113 * t19) * t111, t19 * t61 + t20 * t62 + t35 * t99, t32 * t69, -t69 * t30 - t32 * t67, -t32 * t114 + t69 * t74, t30 * t67, t30 * t114 - t67 * t74, -t140, -t6 * t114 + t26 * t67 + t28 * t74 + t84 * t30, t7 * t114 + t26 * t69 - t29 * t74 + t84 * t32, -t28 * t32 - t29 * t30 - t6 * t69 - t7 * t67, t26 * t84 + t6 * t28 + t7 * t29, t18 * t41, -t41 * t16 - t18 * t39, -t18 * t114 + t41 * t74, t16 * t39, t16 * t114 - t39 * t74, -t140, -t1 * t114 + t13 * t39 + t51 * t16 + t8 * t74, t2 * t114 + t13 * t41 + t51 * t18 - t9 * t74, -t1 * t41 - t9 * t16 - t8 * t18 - t2 * t39, t1 * t8 + t13 * t51 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t102, t127, 0, t104, 0, 0, pkin(2) * t157, pkin(2) * t159, 0.2e1 * (t102 + t104) * pkin(9), pkin(2) ^ 2 + t104 * t117 + t98, t103 * t102, -0.2e1 * t102 * t133, t130 * t159, t101 * t102, t110 * t127, t104, 0.2e1 * t110 * t151 - 0.2e1 * t61 * t114, 0.2e1 * t113 * t151 + 0.2e1 * t62 * t114, 0.2e1 * (-t110 * t62 - t113 * t61) * t111, t61 ^ 2 + t62 ^ 2 + t98, t69 ^ 2, -0.2e1 * t69 * t67, t69 * t158, t67 ^ 2, -t67 * t158, t104, -0.2e1 * t28 * t114 + 0.2e1 * t84 * t67, 0.2e1 * t29 * t114 + 0.2e1 * t84 * t69, -0.2e1 * t28 * t69 - 0.2e1 * t29 * t67, t28 ^ 2 + t29 ^ 2 + t84 ^ 2, t41 ^ 2, -0.2e1 * t41 * t39, t41 * t158, t39 ^ 2, t39 * t157, t104, -0.2e1 * t8 * t114 + 0.2e1 * t51 * t39, 0.2e1 * t9 * t114 + 0.2e1 * t51 * t41, -0.2e1 * t9 * t39 - 0.2e1 * t8 * t41, t51 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, -t74, -t136, t37, -t38, 0, 0, t141, -t110 * t53 + t55 * t113, t146, -t142, t145, 0, -pkin(3) * t53 - pkin(10) * t146 - t143, -pkin(3) * t55 - pkin(10) * t145 + t144 (t141 - t142) * pkin(10) + t121, -t35 * pkin(3) + pkin(10) * t121, t32 * t81, -t81 * t30 - t32 * t79, t81 * t74, t30 * t79, -t79 * t74, 0, t26 * t79 + t97 * t30 + t57 * t74, t26 * t81 + t97 * t32 - t58 * t74, -t58 * t30 - t57 * t32 - t6 * t81 - t7 * t79, t26 * t97 + t6 * t57 + t7 * t58, t18 * t48, -t48 * t16 - t18 * t46, t48 * t74, t16 * t46, -t46 * t74, 0, t13 * t46 + t63 * t16 + t23 * t74, t13 * t48 + t63 * t18 - t24 * t74, -t1 * t48 - t24 * t16 - t23 * t18 - t2 * t46, t1 * t23 + t13 * t63 + t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, t114, 0, -t99, -t114 * pkin(9), 0, 0, t122 (-t101 + t103) * t111, -t132, -t122, -t130, 0, -pkin(9) * t131 + (-pkin(3) * t111 + pkin(10) * t114) * t110, pkin(10) * t130 + (t152 - t153) * t111, t120, -pkin(3) * t99 + pkin(10) * t120, t69 * t81, -t81 * t67 - t69 * t79, -t81 * t114, t67 * t79, t79 * t114, 0, -t57 * t114 + t97 * t67 + t84 * t79, t58 * t114 + t97 * t69 + t84 * t81, -t28 * t81 - t29 * t79 - t57 * t69 - t58 * t67, t28 * t57 + t29 * t58 + t84 * t97, t41 * t48, -t48 * t39 - t41 * t46, -t48 * t114, t39 * t46, t46 * t114, 0, -t23 * t114 + t63 * t39 + t51 * t46, t24 * t114 + t63 * t41 + t51 * t48, -t23 * t41 - t24 * t39 - t9 * t46 - t8 * t48, t8 * t23 + t9 * t24 + t51 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t101, 0.2e1 * t133, 0, t103, 0, 0, 0.2e1 * t153, -0.2e1 * pkin(3) * t110, 0.2e1 * t129 * pkin(10), pkin(10) ^ 2 * t129 + pkin(3) ^ 2, t81 ^ 2, -0.2e1 * t81 * t79, 0, t79 ^ 2, 0, 0, t79 * t161, t81 * t161, -0.2e1 * t57 * t81 - 0.2e1 * t58 * t79, t57 ^ 2 + t58 ^ 2 + t97 ^ 2, t48 ^ 2, -0.2e1 * t48 * t46, 0, t46 ^ 2, 0, 0, t46 * t165, t48 * t165, -0.2e1 * t23 * t48 - 0.2e1 * t24 * t46, t23 ^ 2 + t24 ^ 2 + t63 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, 0, -t53, t74, t19, -t20, 0, 0, 0, 0, t32, 0, -t30, t74, t149 * t74 + t6, -t150 * t74 - t7 (-t105 * t30 - t107 * t32) * pkin(4) (t105 * t7 + t107 * t6) * pkin(4), 0, 0, t18, 0, -t16, t74, t70 * t74 + t1, -t71 * t74 - t2, -t71 * t16 - t70 * t18, t1 * t70 + t2 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, 0, -t134, -t114, t61, -t62, 0, 0, 0, 0, t69, 0, -t67, -t114, -t107 * t148 + t28, t105 * t148 - t29 (-t105 * t67 - t107 * t69) * pkin(4) (t105 * t29 + t107 * t28) * pkin(4), 0, 0, t41, 0, -t39, -t114, -t70 * t114 + t8, t71 * t114 - t9, -t71 * t39 - t70 * t41, t8 * t70 + t9 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110, 0, t113, 0, -t110 * pkin(10), -t113 * pkin(10), 0, 0, 0, 0, t81, 0, -t79, 0, t57, -t58 (-t105 * t79 - t107 * t81) * pkin(4) (t105 * t58 + t107 * t57) * pkin(4), 0, 0, t48, 0, -t46, 0, t23, -t24, -t71 * t46 - t70 * t48, t23 * t70 + t24 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t149, -0.2e1 * t150, 0 (t105 ^ 2 + t107 ^ 2) * pkin(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t71, 0, t70 ^ 2 + t71 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t32, 0, t26, 0, 0, 0, 0, 0, 0, t16, t18, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t69, 0, t84, 0, 0, 0, 0, 0, 0, t39, t41, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t81, 0, t97, 0, 0, 0, 0, 0, 0, t46, t48, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, t74, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, -t39, -t114, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t46, 0, t23, -t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
