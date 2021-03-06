% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRP6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 10:16:51
% EndTime: 2019-05-05 10:17:00
% DurationCPUTime: 2.85s
% Computational Cost: add. (1947->243), mult. (4981->430), div. (0->0), fcn. (5800->12), ass. (0->148)
t92 = sin(qJ(5));
t84 = t92 ^ 2;
t96 = cos(qJ(5));
t86 = t96 ^ 2;
t178 = t84 + t86;
t90 = cos(pkin(7));
t99 = cos(qJ(2));
t139 = t90 * t99;
t88 = sin(pkin(7));
t94 = sin(qJ(3));
t144 = t88 * t94;
t89 = sin(pkin(6));
t91 = cos(pkin(6));
t95 = sin(qJ(2));
t98 = cos(qJ(3));
t36 = t91 * t144 + (t94 * t139 + t95 * t98) * t89;
t141 = t89 * t99;
t52 = -t88 * t141 + t91 * t90;
t93 = sin(qJ(4));
t97 = cos(qJ(4));
t23 = t36 * t97 + t52 * t93;
t142 = t89 * t95;
t143 = t88 * t98;
t34 = -t89 * t98 * t139 + t94 * t142 - t91 * t143;
t11 = t23 * t92 - t34 * t96;
t13 = t23 * t96 + t34 * t92;
t118 = t11 * t92 + t13 * t96;
t56 = t97 * t144 + t93 * t90;
t38 = t96 * t143 + t92 * t56;
t150 = t38 * t96;
t40 = -t92 * t143 + t96 * t56;
t33 = t40 * t92;
t177 = (t33 + t150) * t93;
t21 = t36 * t93 - t52 * t97;
t19 = t21 ^ 2;
t176 = t34 ^ 2;
t175 = t38 ^ 2;
t54 = t93 * t144 - t97 * t90;
t53 = t54 ^ 2;
t174 = -0.2e1 * t56;
t112 = -t96 * pkin(5) - t92 * qJ(6);
t64 = -pkin(4) + t112;
t173 = -0.2e1 * t64;
t172 = 0.2e1 * t88;
t171 = -0.2e1 * t93;
t170 = 0.2e1 * t93;
t169 = pkin(2) * t94;
t168 = pkin(2) * t98;
t167 = pkin(4) * t92;
t166 = pkin(4) * t96;
t165 = pkin(10) * t92;
t164 = pkin(11) * t40;
t130 = t54 * qJ(6);
t70 = pkin(9) * t144;
t47 = t70 + (-pkin(3) - t168) * t90;
t24 = t54 * pkin(4) - t56 * pkin(11) + t47;
t125 = pkin(9) * t143;
t48 = t125 + (pkin(10) + t169) * t90;
t49 = (-pkin(3) * t98 - pkin(10) * t94 - pkin(2)) * t88;
t28 = t97 * t48 + t93 * t49;
t26 = -pkin(11) * t143 + t28;
t6 = t92 * t24 + t96 * t26;
t2 = t130 + t6;
t163 = t2 * t96;
t162 = t54 * pkin(5);
t161 = t6 * t96;
t85 = t93 ^ 2;
t160 = t85 * pkin(10);
t159 = t92 * pkin(11);
t158 = t93 * pkin(10);
t157 = t96 * pkin(11);
t155 = t21 * t92;
t154 = t21 * t93;
t153 = t21 * t96;
t27 = -t93 * t48 + t97 * t49;
t25 = pkin(4) * t143 - t27;
t152 = t25 * t92;
t151 = t25 * t96;
t149 = t40 * t38;
t148 = t54 * t38;
t147 = t54 * t97;
t146 = t56 * t93;
t82 = t88 ^ 2;
t145 = t82 * t98;
t140 = t90 * t94;
t50 = t92 * t54;
t138 = t92 * t93;
t137 = t92 * t96;
t136 = t92 * t97;
t135 = t93 * t54;
t134 = t93 * t97;
t133 = t96 * t54;
t78 = t96 * t93;
t132 = t96 * t97;
t65 = -t97 * pkin(4) - t93 * pkin(11) - pkin(3);
t45 = pkin(10) * t132 + t92 * t65;
t131 = t178 * pkin(11) ^ 2;
t129 = t97 * qJ(6);
t128 = 0.2e1 * t143;
t127 = pkin(11) * t50;
t126 = pkin(11) * t133;
t124 = t38 * t138;
t123 = t92 * t134;
t122 = t93 * t143;
t121 = t85 * t137;
t120 = t97 * t143;
t119 = t11 * t40 - t13 * t38;
t117 = -t11 * t54 + t21 * t38;
t116 = t11 * t97 + t21 * t138;
t115 = -t96 * t24 + t92 * t26;
t114 = t11 ^ 2 + t13 ^ 2 + t19;
t113 = t118 * pkin(11);
t111 = -pkin(5) * t92 + t96 * qJ(6);
t110 = t13 * t54 - t21 * t40;
t109 = t23 * t97 + t154;
t108 = -t27 * t93 + t28 * t97;
t106 = t92 * t38 - t40 * t96;
t41 = -t129 + t45;
t60 = t96 * t65;
t42 = -t60 + (pkin(5) + t165) * t97;
t105 = t41 * t96 + t42 * t92;
t44 = -pkin(10) * t136 + t60;
t104 = -t44 * t92 + t45 * t96;
t103 = t92 * t135 - t97 * t38;
t102 = (t11 * t96 - t13 * t92) * t93;
t101 = pkin(10) ^ 2;
t87 = t97 ^ 2;
t80 = t85 * t101;
t77 = t86 * t85;
t76 = t84 * t85;
t74 = t82 * t98 ^ 2;
t72 = pkin(11) * t136;
t69 = t92 * t78;
t66 = t132 * t171;
t63 = 0.2e1 * t178 * pkin(11);
t61 = (t84 - t86) * t93;
t58 = pkin(2) * t140 + t125;
t57 = t90 * t168 - t70;
t51 = (pkin(10) - t111) * t93;
t37 = t40 ^ 2;
t32 = pkin(11) * t150;
t30 = t40 * t78;
t29 = 0.2e1 * t40 * t54;
t20 = -t40 * t97 + t54 * t78;
t7 = t38 * pkin(5) - t40 * qJ(6) + t25;
t3 = t115 - t162;
t1 = t13 * t97 + t21 * t78;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 ^ 2 + (t95 ^ 2 + t99 ^ 2) * t89 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36 ^ 2 + t52 ^ 2 + t176, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 + t176 + t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, -t142, 0, 0, 0, 0, 0, 0, 0, 0, -t52 * t143 - t34 * t90, t52 * t144 - t36 * t90 (t34 * t94 + t36 * t98) * t88, -t52 * t88 * pkin(2) - t34 * t57 + t36 * t58, 0, 0, 0, 0, 0, 0, t143 * t21 + t34 * t54, t143 * t23 + t34 * t56, t21 * t56 - t23 * t54, -t21 * t27 + t23 * t28 + t34 * t47, 0, 0, 0, 0, 0, 0, t117, -t110, t119, t11 * t115 + t13 * t6 + t21 * t25, 0, 0, 0, 0, 0, 0, t117, t119, t110, t11 * t3 + t13 * t2 + t21 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t82 * t94 ^ 2, 0.2e1 * t94 * t145, t140 * t172, t74, t90 * t128, t90 ^ 2, 0.2e1 * pkin(2) * t145 + 0.2e1 * t57 * t90, -0.2e1 * t82 * t169 - 0.2e1 * t58 * t90 (-t57 * t94 + t58 * t98) * t172, t82 * pkin(2) ^ 2 + t57 ^ 2 + t58 ^ 2, t56 ^ 2, t54 * t174, t143 * t174, t53, t54 * t128, t74, -0.2e1 * t143 * t27 + 0.2e1 * t47 * t54, 0.2e1 * t143 * t28 + 0.2e1 * t47 * t56, -0.2e1 * t27 * t56 - 0.2e1 * t28 * t54, t27 ^ 2 + t28 ^ 2 + t47 ^ 2, t37, -0.2e1 * t149, t29, t175, -0.2e1 * t148, t53, -0.2e1 * t115 * t54 + 0.2e1 * t25 * t38, 0.2e1 * t25 * t40 - 0.2e1 * t6 * t54, 0.2e1 * t115 * t40 - 0.2e1 * t6 * t38, t115 ^ 2 + t25 ^ 2 + t6 ^ 2, t37, t29, 0.2e1 * t149, t53, 0.2e1 * t148, t175, -0.2e1 * t3 * t54 + 0.2e1 * t7 * t38, -0.2e1 * t2 * t38 + 0.2e1 * t3 * t40, 0.2e1 * t2 * t54 - 0.2e1 * t7 * t40, t2 ^ 2 + t3 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t36, 0, 0, 0, 0, 0, 0, 0, 0, -t34 * t97, t34 * t93, t109, -t34 * pkin(3) + pkin(10) * t109, 0, 0, 0, 0, 0, 0, t116, t1, t102, pkin(10) * t154 - t11 * t44 + t13 * t45, 0, 0, 0, 0, 0, 0, t116, t102, -t1, t11 * t42 + t13 * t41 + t21 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t144, 0, t143, t90, t57, -t58, 0, 0, t146, t56 * t97 - t135, -t122, -t147, -t120, 0, -pkin(3) * t54 + pkin(10) * t122 - t47 * t97, -pkin(3) * t56 + pkin(10) * t120 + t47 * t93 (t146 - t147) * pkin(10) + t108, -t47 * pkin(3) + pkin(10) * t108, t30, -t177, t20, t124, -t103, -t147, t44 * t54 + t115 * t97 + (pkin(10) * t38 + t152) * t93, -t45 * t54 + t6 * t97 + (pkin(10) * t40 + t151) * t93, -t45 * t38 - t44 * t40 + (t115 * t96 - t6 * t92) * t93, -t115 * t44 + t158 * t25 + t6 * t45, t30, t20, t177, -t147, t103, t124, t138 * t7 + t3 * t97 + t51 * t38 - t42 * t54, -t41 * t38 + t42 * t40 + (-t2 * t92 + t3 * t96) * t93, -t2 * t97 - t51 * t40 + t41 * t54 - t7 * t78, t2 * t41 + t3 * t42 + t7 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t85, 0.2e1 * t134, 0, t87, 0, 0, 0.2e1 * pkin(3) * t97, pkin(3) * t171, 0.2e1 * (t85 + t87) * pkin(10), pkin(3) ^ 2 + t87 * t101 + t80, t77, -0.2e1 * t121, t66, t76, 0.2e1 * t123, t87, 0.2e1 * t92 * t160 - 0.2e1 * t44 * t97, 0.2e1 * t96 * t160 + 0.2e1 * t45 * t97 (-t44 * t96 - t45 * t92) * t170, t44 ^ 2 + t45 ^ 2 + t80, t77, t66, 0.2e1 * t121, t87, -0.2e1 * t123, t76, 0.2e1 * t138 * t51 + 0.2e1 * t42 * t97 (-t41 * t92 + t42 * t96) * t170, -0.2e1 * t41 * t97 - 0.2e1 * t51 * t78, t41 ^ 2 + t42 ^ 2 + t51 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, -t23, 0, 0, 0, 0, 0, 0, 0, 0, -t153, t155, t118, -t21 * pkin(4) + t113, 0, 0, 0, 0, 0, 0, -t153, t118, -t155, t21 * t64 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, 0, -t54, -t143, t27, -t28, 0, 0, t33, -t106, t50, -t150, t133, 0, -pkin(4) * t38 - t127 - t151, -pkin(4) * t40 - t126 + t152, t161 - t32 + (t115 + t164) * t92, -t25 * pkin(4) + (t115 * t92 + t161) * pkin(11), t33, t50, t106, 0, -t133, -t150, t64 * t38 - t7 * t96 - t127, t163 - t32 + (t3 + t164) * t92, -t64 * t40 - t7 * t92 + t126, t7 * t64 + (t3 * t92 + t163) * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, 0, t97, 0, -t158, -t97 * pkin(10), 0, 0, t69, -t61, -t136, -t69, -t132, 0, t72 + (-pkin(10) * t96 - t167) * t93, pkin(11) * t132 + (t165 - t166) * t93, t104, -pkin(4) * t158 + pkin(11) * t104, t69, -t136, t61, 0, t132, -t69, t138 * t64 - t51 * t96 + t72, t105, -t51 * t92 + (-pkin(11) * t97 - t64 * t93) * t96, pkin(11) * t105 + t51 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t84, 0.2e1 * t137, 0, t86, 0, 0, 0.2e1 * t166, -0.2e1 * t167, t63, pkin(4) ^ 2 + t131, t84, 0, -0.2e1 * t137, 0, 0, t86, t96 * t173, t63, t92 * t173, t64 ^ 2 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t13, -t11 * pkin(5) + t13 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t38, t54, -t115, -t6, 0, 0, 0, t40, 0, t54, t38, 0, -t115 + 0.2e1 * t162, -pkin(5) * t40 - t38 * qJ(6), 0.2e1 * t130 + t6, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, -t138, -t97, t44, -t45, 0, 0, 0, t78, 0, -t97, t138, 0, t60 + (-0.2e1 * pkin(5) - t165) * t97, t112 * t93, -0.2e1 * t129 + t45, -t42 * pkin(5) + t41 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, t96, 0, -t159, -t157, 0, 0, 0, t92, 0, 0, -t96, 0, -t159, t111, t157, t111 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t40, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t78, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, 0, t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
