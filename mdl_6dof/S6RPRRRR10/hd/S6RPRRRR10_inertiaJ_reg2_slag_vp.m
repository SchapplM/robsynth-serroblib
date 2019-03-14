% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRR10_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_inertiaJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t113 = cos(qJ(6));
t101 = t113 ^ 2;
t109 = sin(qJ(6));
t99 = t109 ^ 2;
t138 = t101 + t99;
t110 = sin(qJ(5));
t159 = t110 * pkin(4);
t92 = pkin(12) + t159;
t153 = t138 * t92;
t108 = cos(pkin(6));
t112 = sin(qJ(3));
t116 = cos(qJ(3));
t105 = sin(pkin(6));
t106 = cos(pkin(13));
t107 = cos(pkin(7));
t133 = t106 * t107;
t127 = t105 * t133;
t104 = sin(pkin(7));
t135 = t104 * t116;
t103 = sin(pkin(13));
t137 = t103 * t105;
t54 = -t108 * t135 + t112 * t137 - t116 * t127;
t111 = sin(qJ(4));
t115 = cos(qJ(4));
t160 = pkin(1) * t108;
t86 = t106 * t160;
t57 = t108 * pkin(2) + t86 + (-pkin(9) * t107 - qJ(2)) * t137;
t63 = (-pkin(9) * t103 * t104 - pkin(2) * t106 - pkin(1)) * t105;
t37 = -t104 * t57 + t107 * t63;
t136 = t104 * t112;
t56 = t108 * t136 + (t103 * t116 + t112 * t133) * t105;
t24 = t54 * pkin(3) - t56 * pkin(10) + t37;
t122 = t104 * t63 + t107 * t57;
t134 = t105 * t106;
t69 = qJ(2) * t134 + t103 * t160;
t50 = (t104 * t108 + t127) * pkin(9) + t69;
t34 = t122 * t112 + t116 * t50;
t66 = t104 * t134 - t108 * t107;
t26 = -t66 * pkin(10) + t34;
t14 = -t111 * t26 + t115 * t24;
t162 = t54 * pkin(4);
t40 = -t66 * t111 + t56 * t115;
t11 = -t40 * pkin(11) + t14 + t162;
t114 = cos(qJ(5));
t15 = t111 * t24 + t115 * t26;
t39 = -t56 * t111 - t66 * t115;
t13 = t39 * pkin(11) + t15;
t145 = t114 * t13;
t8 = t110 * t11 + t145;
t6 = t54 * pkin(12) + t8;
t33 = -t112 * t50 + t122 * t116;
t25 = t66 * pkin(3) - t33;
t16 = -t39 * pkin(4) + t25;
t29 = t110 * t40 - t114 * t39;
t31 = t110 * t39 + t114 * t40;
t9 = t29 * pkin(5) - t31 * pkin(12) + t16;
t3 = t109 * t9 + t113 * t6;
t1 = t3 * t113;
t2 = -t109 * t6 + t113 * t9;
t124 = -t2 * t109 + t1;
t173 = t29 ^ 2;
t70 = t115 * t107 - t111 * t136;
t71 = t111 * t107 + t115 * t136;
t44 = t110 * t71 - t114 * t70;
t172 = t44 ^ 2;
t53 = t54 ^ 2;
t163 = -pkin(11) - pkin(10);
t128 = t163 * t111;
t81 = t163 * t115;
t60 = -t110 * t81 - t114 * t128;
t171 = t60 ^ 2;
t75 = t110 * t111 - t114 * t115;
t170 = t75 ^ 2;
t169 = -0.2e1 * t29;
t168 = 0.2e1 * t54;
t167 = -0.2e1 * t56;
t94 = -t115 * pkin(4) - pkin(3);
t166 = 0.2e1 * t94;
t165 = 0.2e1 * t105;
t164 = 0.2e1 * t115;
t158 = t114 * pkin(4);
t93 = -pkin(5) - t158;
t161 = pkin(5) - t93;
t156 = t29 * t75;
t155 = t44 * t60;
t126 = -t114 * t11 + t110 * t13;
t5 = -t54 * pkin(5) + t126;
t154 = t5 * t113;
t152 = t138 * pkin(12);
t98 = t105 ^ 2;
t151 = t106 * t98;
t27 = t109 * t29;
t77 = t110 * t115 + t114 * t111;
t150 = t109 * t77;
t149 = t109 * t92;
t148 = t111 * t54;
t28 = t113 * t29;
t147 = t113 * t77;
t146 = t113 * t92;
t144 = t115 * t54;
t19 = t109 * t31 - t54 * t113;
t17 = t19 * t113;
t21 = t54 * t109 + t113 * t31;
t18 = t21 * t109;
t143 = t21 * t113;
t142 = t39 * t115;
t141 = t40 * t111;
t43 = t44 * t109;
t140 = t44 * t113;
t58 = t60 * t109;
t139 = t60 * t113;
t132 = t109 * t113;
t100 = t111 ^ 2;
t102 = t115 ^ 2;
t131 = t100 + t102;
t130 = -0.2e1 * t77 * t75;
t129 = t108 * t165;
t125 = -pkin(5) * t77 - pkin(12) * t75;
t123 = -t75 * t92 + t77 * t93;
t51 = t75 * pkin(5) - t77 * pkin(12) + t94;
t62 = t110 * t128 - t114 * t81;
t35 = -t109 * t62 + t113 * t51;
t36 = t109 * t51 + t113 * t62;
t22 = -t35 * t109 + t36 * t113;
t46 = t110 * t70 + t114 * t71;
t41 = -t109 * t46 - t113 * t135;
t42 = -t109 * t135 + t113 * t46;
t32 = -t41 * t109 + t42 * t113;
t121 = -t14 * t111 + t15 * t115;
t120 = -t70 * t111 + t71 * t115;
t97 = t104 ^ 2;
t90 = t97 * t116 ^ 2;
t87 = 0.2e1 * t132;
t74 = t77 ^ 2;
t73 = t113 * t75;
t72 = t109 * t75;
t68 = -qJ(2) * t137 + t86;
t65 = t77 * t132;
t52 = (t101 - t99) * t77;
t12 = -t109 * t19 + t143;
t4 = t5 * t109;
t7 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t98 * t103 ^ 2, 0.2e1 * t103 * t151, t103 * t129, t98 * t106 ^ 2, t106 * t129, t108 ^ 2, 0.2e1 * pkin(1) * t151 + 0.2e1 * t68 * t108, -0.2e1 * t98 * pkin(1) * t103 - 0.2e1 * t69 * t108 (-t103 * t68 + t106 * t69) * t165, t98 * pkin(1) ^ 2 + t68 ^ 2 + t69 ^ 2, t56 ^ 2, t54 * t167, t66 * t167, t53, t66 * t168, t66 ^ 2, -0.2e1 * t33 * t66 + 0.2e1 * t37 * t54, 0.2e1 * t34 * t66 + 0.2e1 * t37 * t56, -0.2e1 * t33 * t56 - 0.2e1 * t34 * t54, t33 ^ 2 + t34 ^ 2 + t37 ^ 2, t40 ^ 2, 0.2e1 * t40 * t39, t40 * t168, t39 ^ 2, t39 * t168, t53, 0.2e1 * t14 * t54 - 0.2e1 * t25 * t39, -0.2e1 * t15 * t54 + 0.2e1 * t25 * t40, -0.2e1 * t14 * t40 + 0.2e1 * t15 * t39, t14 ^ 2 + t15 ^ 2 + t25 ^ 2, t31 ^ 2, t31 * t169, t31 * t168, t173, t54 * t169, t53, -0.2e1 * t126 * t54 + 0.2e1 * t16 * t29, 0.2e1 * t16 * t31 - 0.2e1 * t8 * t54, 0.2e1 * t126 * t31 - 0.2e1 * t8 * t29, t126 ^ 2 + t16 ^ 2 + t8 ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, 0.2e1 * t21 * t29, t19 ^ 2, t19 * t169, t173, 0.2e1 * t5 * t19 + 0.2e1 * t2 * t29, 0.2e1 * t5 * t21 - 0.2e1 * t3 * t29, -0.2e1 * t3 * t19 - 0.2e1 * t2 * t21, t2 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, t137, 0, -t105 * pkin(1), 0, 0, 0, 0, 0, 0, t107 * t54 - t66 * t135, t107 * t56 + t66 * t136 (-t112 * t54 - t116 * t56) * t104, t37 * t107 + (t112 * t34 + t116 * t33) * t104, 0, 0, 0, 0, 0, 0, t135 * t39 + t70 * t54, -t135 * t40 - t71 * t54, t71 * t39 - t70 * t40, -t135 * t25 + t14 * t70 + t15 * t71, 0, 0, 0, 0, 0, 0, -t135 * t29 - t44 * t54, -t135 * t31 - t46 * t54, -t46 * t29 + t44 * t31, t126 * t44 - t135 * t16 + t8 * t46, 0, 0, 0, 0, 0, 0, t44 * t19 + t41 * t29, t44 * t21 - t42 * t29, -t42 * t19 - t41 * t21, t2 * t41 + t3 * t42 + t5 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t112 ^ 2 + t107 ^ 2 + t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 ^ 2 + t71 ^ 2 + t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 + t172 + t90, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 ^ 2 + t42 ^ 2 + t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t54, -t66, t33, -t34, 0, 0, t141, t111 * t39 + t40 * t115, t148, t142, t144, 0, pkin(3) * t39 - pkin(10) * t148 - t25 * t115, -pkin(3) * t40 - pkin(10) * t144 + t25 * t111 (t141 + t142) * pkin(10) + t121, -t25 * pkin(3) + pkin(10) * t121, t31 * t77, -t77 * t29 - t31 * t75, t77 * t54, t156, -t75 * t54, 0, t16 * t75 + t94 * t29 - t60 * t54, t16 * t77 + t94 * t31 - t62 * t54, t126 * t77 - t62 * t29 + t60 * t31 - t8 * t75, t126 * t60 + t16 * t94 + t8 * t62, t77 * t143 (-t18 - t17) * t77, t29 * t147 + t21 * t75, t19 * t150, -t29 * t150 - t19 * t75, t156, t5 * t150 + t60 * t19 + t2 * t75 + t35 * t29, t147 * t5 + t60 * t21 - t36 * t29 - t3 * t75, -t36 * t19 - t35 * t21 + (-t109 * t3 - t113 * t2) * t77, t2 * t35 + t3 * t36 + t5 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, -t136, 0, 0, 0, 0, 0, 0, 0, 0, t115 * t135, -t111 * t135, t120, pkin(3) * t135 + pkin(10) * t120, 0, 0, 0, 0, 0, 0, -t75 * t135, -t77 * t135, t44 * t77 - t46 * t75, -t135 * t94 + t46 * t62 + t155, 0, 0, 0, 0, 0, 0, t41 * t75 + t43 * t77, t140 * t77 - t42 * t75 (-t109 * t42 - t113 * t41) * t77, t41 * t35 + t42 * t36 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t100, t111 * t164, 0, t102, 0, 0, pkin(3) * t164, -0.2e1 * pkin(3) * t111, 0.2e1 * t131 * pkin(10), pkin(10) ^ 2 * t131 + pkin(3) ^ 2, t74, t130, 0, t170, 0, 0, t75 * t166, t77 * t166, 0.2e1 * t60 * t77 - 0.2e1 * t62 * t75, t62 ^ 2 + t94 ^ 2 + t171, t101 * t74, -0.2e1 * t74 * t132, 0.2e1 * t75 * t147, t99 * t74, t109 * t130, t170, 0.2e1 * t35 * t75 + 0.2e1 * t58 * t77, 0.2e1 * t139 * t77 - 0.2e1 * t36 * t75, 0.2e1 * (-t109 * t36 - t113 * t35) * t77, t35 ^ 2 + t36 ^ 2 + t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, t39, t54, t14, -t15, 0, 0, 0, 0, t31, 0, -t29, t54, t54 * t158 - t126, -t145 + (-t11 - t162) * t110 (-t110 * t29 - t114 * t31) * pkin(4) (t110 * t8 - t114 * t126) * pkin(4), t18, t12, t27, -t17, t28, 0, -t149 * t29 + t93 * t19 - t154, -t146 * t29 + t93 * t21 + t4, -t92 * t17 + t1 + (t21 * t92 - t2) * t109, t124 * t92 + t5 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, -t71, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t46, 0 (t110 * t46 - t114 * t44) * pkin(4), 0, 0, 0, 0, 0, 0, -t140, t43, t32, t32 * t92 + t44 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, 0, t115, 0, -t111 * pkin(10), -t115 * pkin(10), 0, 0, 0, 0, t77, 0, -t75, 0, -t60, -t62 (-t110 * t75 - t114 * t77) * pkin(4) (t110 * t62 - t114 * t60) * pkin(4), t65, t52, t72, -t65, t73, 0, t109 * t123 - t139, t113 * t123 + t58, t22, t22 * t92 + t60 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t158, -0.2e1 * t159, 0 (t110 ^ 2 + t114 ^ 2) * pkin(4) ^ 2, t99, t87, 0, t101, 0, 0, -0.2e1 * t93 * t113, 0.2e1 * t93 * t109, 0.2e1 * t153, t138 * t92 ^ 2 + t93 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t29, t54, -t126, -t8, 0, 0, t18, t12, t27, -t17, t28, 0, -pkin(5) * t19 - pkin(12) * t27 - t154, -pkin(5) * t21 - pkin(12) * t28 + t4 (t18 - t17) * pkin(12) + t124, -t5 * pkin(5) + pkin(12) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t46, 0, 0, 0, 0, 0, 0, 0, 0, -t140, t43, t32, -t44 * pkin(5) + pkin(12) * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, -t75, 0, -t60, -t62, 0, 0, t65, t52, t72, -t65, t73, 0, t109 * t125 - t139, t113 * t125 + t58, t22, -t60 * pkin(5) + pkin(12) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t158, -t159, 0, 0, t99, t87, 0, t101, 0, 0, t161 * t113, -t161 * t109, t152 + t153, -t93 * pkin(5) + pkin(12) * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t99, t87, 0, t101, 0, 0, 0.2e1 * pkin(5) * t113, -0.2e1 * pkin(5) * t109, 0.2e1 * t152, pkin(12) ^ 2 * t138 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, -t19, t29, t2, -t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, -t150, t75, t35, -t36, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, t113, 0, -t149, -t146, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, 0, t113, 0, -t109 * pkin(12), -t113 * pkin(12), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t7;