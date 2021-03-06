% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRPRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 13:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRRR5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR5_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR5_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 21:26:12
% EndTime: 2019-05-06 21:26:21
% DurationCPUTime: 3.19s
% Computational Cost: add. (4431->261), mult. (10673->493), div. (0->0), fcn. (12207->12), ass. (0->151)
t102 = cos(qJ(4));
t103 = cos(qJ(2));
t92 = sin(pkin(12));
t93 = sin(pkin(6));
t94 = cos(pkin(12));
t99 = sin(qJ(2));
t56 = (t103 * t92 + t94 * t99) * t93;
t95 = cos(pkin(6));
t98 = sin(qJ(4));
t43 = -t95 * t102 + t56 * t98;
t42 = t43 ^ 2;
t123 = t93 * t103;
t147 = t93 * t99;
t54 = -t94 * t123 + t92 * t147;
t174 = t54 ^ 2;
t173 = -0.2e1 * t43;
t172 = 0.2e1 * t43;
t171 = -0.2e1 * t54;
t101 = cos(qJ(5));
t85 = -t101 * pkin(5) - pkin(4);
t170 = 0.2e1 * t85;
t169 = 0.2e1 * t95;
t168 = 0.2e1 * t98;
t167 = -0.2e1 * t102;
t166 = -pkin(11) - pkin(10);
t165 = pkin(1) * t99;
t164 = t43 * pkin(5);
t142 = pkin(8) + qJ(3);
t160 = t95 * pkin(2);
t76 = t95 * t103 * pkin(1);
t46 = -t142 * t147 + t160 + t76;
t117 = t95 * t165;
t51 = t142 * t123 + t117;
t27 = t92 * t46 + t94 * t51;
t23 = t95 * pkin(9) + t27;
t70 = (-pkin(2) * t103 - pkin(1)) * t93;
t31 = t54 * pkin(3) - t56 * pkin(9) + t70;
t14 = t102 * t31 - t98 * t23;
t9 = -t54 * pkin(4) - t14;
t97 = sin(qJ(5));
t163 = t9 * t97;
t162 = t92 * pkin(2);
t161 = t94 * pkin(2);
t96 = sin(qJ(6));
t159 = t96 * pkin(5);
t100 = cos(qJ(6));
t158 = t100 * pkin(5);
t45 = t56 * t102 + t95 * t98;
t28 = -t54 * t101 + t45 * t97;
t15 = t102 * t23 + t98 * t31;
t10 = t54 * pkin(10) + t15;
t141 = -t94 * t46 + t92 * t51;
t22 = -t95 * pkin(3) + t141;
t13 = t43 * pkin(4) - t45 * pkin(10) + t22;
t7 = t101 * t10 + t97 * t13;
t5 = -t28 * pkin(11) + t7;
t157 = t100 * t5;
t156 = t102 * pkin(4);
t155 = t102 * pkin(5);
t30 = t45 * t101 + t54 * t97;
t18 = t100 * t30 - t96 * t28;
t69 = t100 * t97 + t96 * t101;
t59 = t69 * t98;
t154 = t18 * t59;
t153 = t30 * t97;
t120 = t100 * t101;
t145 = t97 * t98;
t61 = t98 * t120 - t96 * t145;
t152 = t61 * t43;
t151 = t69 * t59;
t79 = pkin(9) + t162;
t150 = t79 * t97;
t88 = t97 ^ 2;
t149 = t88 * t98;
t148 = t9 * t101;
t146 = t97 * t43;
t36 = t98 * t43;
t144 = t98 * t54;
t143 = t98 * t79;
t90 = t101 ^ 2;
t140 = t88 + t90;
t89 = t98 ^ 2;
t91 = t102 ^ 2;
t139 = t89 + t91;
t119 = t101 * t102;
t111 = t79 * t119;
t80 = -pkin(3) - t161;
t66 = -t98 * pkin(10) - t156 + t80;
t35 = t111 + (-pkin(11) * t98 + t66) * t97;
t138 = t100 * t35;
t137 = t101 * t43;
t136 = t101 * t98;
t16 = t100 * t28 + t96 * t30;
t135 = t102 * t16;
t134 = t102 * t28;
t67 = t96 * t97 - t120;
t133 = t102 * t67;
t132 = t102 * t79;
t86 = t93 ^ 2;
t131 = t103 * t86;
t130 = t18 * t102;
t129 = t28 * t101;
t128 = t30 * t101;
t127 = t30 * t102;
t126 = t43 * t102;
t125 = t45 * t102;
t124 = t69 * t102;
t122 = t97 * t101;
t121 = t97 * t102;
t118 = t93 * t169;
t116 = t102 * t168;
t115 = t30 * t145;
t114 = t97 * t36;
t113 = t43 * t136;
t112 = t98 * t122;
t6 = -t97 * t10 + t101 * t13;
t4 = -t30 * pkin(11) + t164 + t6;
t1 = t100 * t4 - t96 * t5;
t110 = t140 * pkin(10);
t62 = t101 * t66;
t34 = -pkin(11) * t136 + t62 + (-pkin(5) - t150) * t102;
t19 = t100 * t34 - t96 * t35;
t109 = -pkin(4) * t98 + pkin(10) * t102;
t108 = t7 * t101 - t6 * t97;
t40 = -t79 * t121 + t62;
t41 = t97 * t66 + t111;
t107 = t41 * t101 - t40 * t97;
t87 = t95 ^ 2;
t84 = t90 * t98;
t83 = t90 * t89;
t81 = t88 * t89;
t77 = t79 ^ 2;
t74 = t166 * t101;
t73 = t166 * t97;
t71 = t89 * t77;
t65 = pkin(8) * t123 + t117;
t64 = -pkin(8) * t147 + t76;
t63 = (pkin(5) * t97 + t79) * t98;
t58 = t61 ^ 2;
t57 = t59 ^ 2;
t52 = t102 * t54;
t50 = -t100 * t74 + t96 * t73;
t49 = t100 * t73 + t96 * t74;
t37 = t61 * t67;
t32 = t59 * t43;
t24 = t98 * t129;
t20 = t96 * t34 + t138;
t11 = t61 * t16;
t8 = t28 * pkin(5) + t9;
t2 = t96 * t4 + t157;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t86 * t99 ^ 2, 0.2e1 * t99 * t131, t99 * t118, t86 * t103 ^ 2, t103 * t118, t87, 0.2e1 * pkin(1) * t131 + 0.2e1 * t64 * t95, -0.2e1 * t86 * t165 - 0.2e1 * t65 * t95, 0.2e1 * (t103 * t65 - t64 * t99) * t93, t86 * pkin(1) ^ 2 + t64 ^ 2 + t65 ^ 2, t56 ^ 2, t56 * t171, t56 * t169, t174, t95 * t171, t87, -0.2e1 * t141 * t95 + 0.2e1 * t70 * t54, -0.2e1 * t27 * t95 + 0.2e1 * t70 * t56, 0.2e1 * t141 * t56 - 0.2e1 * t27 * t54, t141 ^ 2 + t27 ^ 2 + t70 ^ 2, t45 ^ 2, t45 * t173, 0.2e1 * t45 * t54, t42, t43 * t171, t174, 0.2e1 * t14 * t54 + 0.2e1 * t22 * t43, -0.2e1 * t15 * t54 + 0.2e1 * t22 * t45, -0.2e1 * t14 * t45 - 0.2e1 * t15 * t43, t14 ^ 2 + t15 ^ 2 + t22 ^ 2, t30 ^ 2, -0.2e1 * t30 * t28, t30 * t172, t28 ^ 2, t28 * t173, t42, 0.2e1 * t9 * t28 + 0.2e1 * t6 * t43, 0.2e1 * t9 * t30 - 0.2e1 * t7 * t43, -0.2e1 * t7 * t28 - 0.2e1 * t6 * t30, t6 ^ 2 + t7 ^ 2 + t9 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, t18 * t172, t16 ^ 2, t16 * t173, t42, 0.2e1 * t1 * t43 + 0.2e1 * t8 * t16, 0.2e1 * t8 * t18 - 0.2e1 * t2 * t43, -0.2e1 * t1 * t18 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, 0, t123, t95, t64, -t65, 0, 0, 0, 0, t56, 0, -t54, t95, t94 * t160 - t141, -t92 * t160 - t27 (-t54 * t92 - t56 * t94) * pkin(2) (-t141 * t94 + t27 * t92) * pkin(2), t45 * t98, -t36 + t125, t144, -t126, t52, 0, -t22 * t102 - t143 * t54 + t80 * t43, -t132 * t54 + t22 * t98 + t80 * t45 (t45 * t79 - t14) * t98 + (-t43 * t79 + t15) * t102, t22 * t80 + (t15 * t102 - t14 * t98) * t79, t98 * t128, -t24 - t115, t113 - t127, t28 * t145, -t114 + t134, -t126, -t6 * t102 + t40 * t43 + (t28 * t79 + t163) * t98, t7 * t102 - t41 * t43 + (t30 * t79 + t148) * t98, -t41 * t28 - t40 * t30 + (-t101 * t6 - t7 * t97) * t98, t143 * t9 + t6 * t40 + t7 * t41, t18 * t61, -t11 - t154, -t130 + t152, t16 * t59, -t32 + t135, -t126, -t1 * t102 + t63 * t16 + t19 * t43 + t8 * t59, t2 * t102 + t63 * t18 - t20 * t43 + t8 * t61, -t1 * t61 - t20 * t16 - t19 * t18 - t2 * t59, t1 * t19 + t2 * t20 + t8 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t161, -0.2e1 * t162, 0 (t92 ^ 2 + t94 ^ 2) * pkin(2) ^ 2, t89, t116, 0, t91, 0, 0, t80 * t167, t80 * t168, 0.2e1 * t139 * t79, t91 * t77 + t80 ^ 2 + t71, t83, -0.2e1 * t89 * t122, -0.2e1 * t98 * t119, t81, t97 * t116, t91, -0.2e1 * t40 * t102 + 0.2e1 * t150 * t89, 0.2e1 * t89 * t79 * t101 + 0.2e1 * t41 * t102 (-t101 * t40 - t41 * t97) * t168, t40 ^ 2 + t41 ^ 2 + t71, t58, -0.2e1 * t61 * t59, t61 * t167, t57, -t59 * t167, t91, -0.2e1 * t19 * t102 + 0.2e1 * t63 * t59, 0.2e1 * t20 * t102 + 0.2e1 * t63 * t61, -0.2e1 * t19 * t61 - 0.2e1 * t20 * t59, t19 ^ 2 + t20 ^ 2 + t63 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t56, 0, t70, 0, 0, 0, 0, 0, 0, t52, -t144, -t36 - t125, t14 * t102 + t15 * t98, 0, 0, 0, 0, 0, 0, -t114 - t134, -t113 - t127, -t24 + t115, -t9 * t102 + t108 * t98, 0, 0, 0, 0, 0, 0, -t32 - t135, -t130 - t152, -t11 + t154, -t1 * t59 - t8 * t102 + t2 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t107 - t132) * t98, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t102 - t19 * t59 + t20 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 + t81 + t91, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 + t57 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, -t43, t54, t14, -t15, 0, 0, t153, -t97 * t28 + t128, t146, -t129, t137, 0, -pkin(4) * t28 - pkin(10) * t146 - t148, -pkin(4) * t30 - pkin(10) * t137 + t163 (-t129 + t153) * pkin(10) + t108, -t9 * pkin(4) + pkin(10) * t108, t18 * t69, -t69 * t16 - t18 * t67, t69 * t43, t16 * t67, -t67 * t43, 0, t85 * t16 + t49 * t43 + t8 * t67, t85 * t18 - t50 * t43 + t8 * t69, -t1 * t69 - t50 * t16 - t49 * t18 - t2 * t67, t1 * t49 + t2 * t50 + t8 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, 0, t102, 0, -t143, -t132, 0, 0, t112, t84 - t149, -t121, -t112, -t119, 0, t109 * t97 - t136 * t79, t101 * t109 + t143 * t97, t107, -pkin(4) * t143 + pkin(10) * t107, t61 * t69, -t37 - t151, -t124, t59 * t67, t133, 0, -t49 * t102 + t85 * t59 + t63 * t67, t50 * t102 + t85 * t61 + t63 * t69, -t19 * t69 - t20 * t67 - t49 * t61 - t50 * t59, t19 * t49 + t20 * t50 + t63 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, -t98, 0, 0, 0, 0, 0, 0, 0, 0, t119, -t121, t84 + t149, t110 * t98 + t156, 0, 0, 0, 0, 0, 0, -t133, -t124, -t37 + t151, -t102 * t85 - t59 * t49 + t61 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t88, 0.2e1 * t122, 0, t90, 0, 0, 0.2e1 * pkin(4) * t101, -0.2e1 * pkin(4) * t97, 0.2e1 * t110, pkin(10) ^ 2 * t140 + pkin(4) ^ 2, t69 ^ 2, -0.2e1 * t69 * t67, 0, t67 ^ 2, 0, 0, t67 * t170, t69 * t170, -0.2e1 * t49 * t69 - 0.2e1 * t50 * t67, t49 ^ 2 + t50 ^ 2 + t85 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, -t28, t43, t6, -t7, 0, 0, 0, 0, t18, 0, -t16, t43, t158 * t43 + t1, -t157 + (-t4 - t164) * t96 (-t100 * t18 - t16 * t96) * pkin(5) (t1 * t100 + t2 * t96) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, -t145, -t102, t40, -t41, 0, 0, 0, 0, t61, 0, -t59, -t102, -t100 * t155 + t19, -t138 + (-t34 + t155) * t96 (-t100 * t61 - t59 * t96) * pkin(5) (t100 * t19 + t20 * t96) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, -t136, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t61, 0 (-t100 * t59 + t61 * t96) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, 0, t101, 0, -t97 * pkin(10), -t101 * pkin(10), 0, 0, 0, 0, t69, 0, -t67, 0, t49, -t50 (-t100 * t69 - t67 * t96) * pkin(5) (t100 * t49 + t50 * t96) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t158, -0.2e1 * t159, 0 (t100 ^ 2 + t96 ^ 2) * pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, t43, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, -t59, -t102, t19, -t20, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t61, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t67, 0, t49, -t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t158, -t159, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t3;
