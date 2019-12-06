% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5PRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRPRR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:57:54
% EndTime: 2019-12-05 15:57:59
% DurationCPUTime: 1.48s
% Computational Cost: add. (2737->246), mult. (7239->357), div. (0->0), fcn. (5757->10), ass. (0->136)
t163 = cos(qJ(4));
t91 = cos(pkin(10));
t119 = t163 * t91;
t89 = sin(pkin(10));
t94 = sin(qJ(4));
t145 = t94 * t89;
t102 = t119 - t145;
t90 = sin(pkin(5));
t97 = cos(qJ(2));
t147 = t90 * t97;
t100 = t102 * t147;
t144 = pkin(7) + qJ(3);
t74 = t144 * t89;
t75 = t144 * t91;
t103 = -t163 * t74 - t94 * t75;
t142 = -qJD(1) * t100 + qJD(3) * t102 + qJD(4) * t103;
t137 = qJD(1) * t90;
t95 = sin(qJ(2));
t124 = t95 * t137;
t120 = qJD(4) * t145;
t65 = -qJD(4) * t119 + t120;
t71 = t163 * t89 + t94 * t91;
t66 = t71 * qJD(4);
t174 = -pkin(4) * t66 - pkin(8) * t65 + t124;
t135 = qJD(1) * t97;
t123 = t90 * t135;
t69 = (qJD(3) + t123) * qJD(2);
t173 = t69 * t71;
t139 = pkin(7) * qJD(2);
t73 = qJD(2) * qJ(3) + t124;
t92 = cos(pkin(5));
t136 = qJD(1) * t92;
t79 = t91 * t136;
t39 = t79 + (-t73 - t139) * t89;
t50 = t89 * t136 + t91 * t73;
t40 = t91 * t139 + t50;
t105 = -t163 * t39 + t94 * t40;
t171 = t102 * t69;
t11 = -qJD(4) * t105 + t171;
t20 = t163 * t40 + t94 * t39;
t16 = qJD(4) * pkin(8) + t20;
t134 = qJD(2) * t71;
t113 = qJD(3) - t123;
t83 = -pkin(3) * t91 - pkin(2);
t133 = qJD(2) * t83;
t56 = t113 + t133;
t122 = qJD(2) * t145;
t80 = qJD(2) * t119;
t61 = -t80 + t122;
t21 = pkin(4) * t61 - pkin(8) * t134 + t56;
t93 = sin(qJ(5));
t96 = cos(qJ(5));
t111 = t16 * t93 - t21 * t96;
t77 = qJD(4) * t80;
t54 = qJD(2) * t120 - t77;
t55 = qJD(2) * t66;
t132 = qJD(2) * t95;
t121 = t90 * t132;
t76 = qJD(1) * t121;
t26 = pkin(4) * t55 + pkin(8) * t54 + t76;
t1 = -t111 * qJD(5) + t11 * t96 + t26 * t93;
t57 = qJD(5) + t61;
t114 = t111 * t57 + t1;
t8 = t16 * t96 + t21 * t93;
t2 = -qJD(5) * t8 - t11 * t93 + t96 * t26;
t172 = t57 * t8 + t2;
t110 = (-t73 * t89 + t79) * t89 - t50 * t91;
t170 = t110 * t97;
t117 = t57 * t93;
t48 = qJD(4) * t93 + t134 * t96;
t169 = t48 * t117;
t131 = qJD(5) * t48;
t28 = -t54 * t93 + t131;
t168 = t134 ^ 2;
t38 = -pkin(4) * t102 - pkin(8) * t71 + t83;
t43 = t163 * t75 - t94 * t74;
t13 = t38 * t96 - t43 * t93;
t165 = t13 * qJD(5) + t142 * t96 - t174 * t93;
t14 = t38 * t93 + t43 * t96;
t164 = t14 * qJD(5) + t142 * t93 + t174 * t96;
t148 = t90 * t95;
t59 = -t89 * t148 + t91 * t92;
t60 = t91 * t148 + t89 * t92;
t104 = t163 * t59 - t94 * t60;
t12 = t20 * qJD(4) + t173;
t162 = t12 * t104;
t161 = t12 * t103;
t128 = t96 * qJD(4);
t130 = qJD(5) * t93;
t27 = -qJD(5) * t128 + t130 * t134 + t96 * t54;
t160 = t27 * t93;
t159 = t28 * t96;
t46 = t134 * t93 - t128;
t158 = t46 * t61;
t157 = t46 * t134;
t156 = t48 * t46;
t155 = t48 * t134;
t153 = t55 * t102;
t152 = t55 * t93;
t151 = t134 * t61;
t150 = t71 * t93;
t149 = t71 * t96;
t98 = qJD(2) ^ 2;
t146 = t90 * t98;
t25 = t93 * t28;
t52 = t96 * t55;
t129 = qJD(5) * t96;
t143 = -t46 * t129 - t25;
t99 = t71 * t147;
t141 = -qJD(1) * t99 + qJD(3) * t71 + qJD(4) * t43;
t140 = t89 ^ 2 + t91 ^ 2;
t138 = qJD(2) * pkin(2);
t127 = t95 * t146;
t126 = t97 * t146;
t125 = t90 ^ 2 * t135;
t118 = t140 * t69;
t116 = t57 * t96;
t115 = -t111 * t96 + t8 * t93;
t109 = -t61 * t117 - t57 * t130 + t52;
t34 = t163 * t60 + t94 * t59;
t23 = -t96 * t147 - t34 * t93;
t108 = t93 * t147 - t34 * t96;
t107 = t71 * t129 - t65 * t93;
t106 = -t71 * t130 - t65 * t96;
t15 = -qJD(4) * pkin(4) + t105;
t101 = -pkin(8) * t55 + t57 * t15;
t72 = t113 - t138;
t58 = t61 ^ 2;
t36 = pkin(4) * t134 + pkin(8) * t61;
t18 = qJD(2) * t99 + qJD(4) * t34;
t17 = qJD(2) * t100 + qJD(4) * t104;
t10 = -t105 * t96 + t36 * t93;
t9 = t105 * t93 + t36 * t96;
t6 = qJD(5) * t108 + t96 * t121 - t17 * t93;
t5 = qJD(5) * t23 + t93 * t121 + t17 * t96;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127, -t126, 0, 0, 0, 0, 0, 0, 0, 0, -t91 * t127, t89 * t127, t140 * t126, (-t59 * t89 + t60 * t91) * t69 + (-t95 * t125 + (t72 * t95 - t170) * t90) * qJD(2), 0, 0, 0, 0, 0, 0, -qJD(4) * t18 + (t61 * t132 - t55 * t97) * t90, -qJD(4) * t17 + (t132 * t134 + t54 * t97) * t90, t104 * t54 + t134 * t18 - t17 * t61 - t34 * t55, t11 * t34 - t162 + t17 * t20 + t18 * t105 + (t56 * t90 - t125) * t132, 0, 0, 0, 0, 0, 0, -t104 * t28 + t18 * t46 + t23 * t55 + t57 * t6, t104 * t27 + t108 * t55 + t18 * t48 - t5 * t57, t108 * t28 + t23 * t27 - t46 * t5 - t48 * t6, -t1 * t108 - t111 * t6 + t15 * t18 + t2 * t23 + t5 * t8 - t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * qJD(2) * t140 + t118, -t110 * qJD(3) + qJ(3) * t118 + (t170 + (-t72 - t138) * t95) * t137, -t134 * t65 - t54 * t71, -t102 * t54 - t134 * t66 - t55 * t71 + t61 * t65, -t65 * qJD(4), t61 * t66 - t153, -t66 * qJD(4), 0, t55 * t83 + t56 * t66 - t141 * qJD(4) + (-qJD(2) * t102 - t61) * t124, -t142 * qJD(4) - t54 * t83 - t56 * t65, t102 * t11 + t103 * t54 - t105 * t65 + t12 * t71 + t134 * t141 - t142 * t61 - t20 * t66 - t43 * t55, t11 * t43 - t161 + t142 * t20 + t141 * t105 + (-t56 + t133) * t124, t106 * t48 - t149 * t27, (t46 * t96 + t48 * t93) * t65 + (t160 - t159 + (t46 * t93 - t48 * t96) * qJD(5)) * t71, t102 * t27 + t106 * t57 + t48 * t66 + t71 * t52, t107 * t46 + t71 * t25, t102 * t28 - t107 * t57 - t150 * t55 - t46 * t66, t57 * t66 - t153, -t102 * t2 - t103 * t28 + t107 * t15 - t111 * t66 + t12 * t150 + t13 * t55 + t141 * t46 - t164 * t57, t1 * t102 + t103 * t27 + t106 * t15 + t12 * t149 - t14 * t55 + t141 * t48 - t165 * t57 - t66 * t8, t13 * t27 - t14 * t28 + t115 * t65 + t164 * t48 - t165 * t46 + (-t1 * t93 - t2 * t96 + (-t111 * t93 - t8 * t96) * qJD(5)) * t71, t1 * t14 + t111 * t164 + t13 * t2 + t141 * t15 + t165 * t8 - t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t98, qJD(2) * t110 + t76, 0, 0, 0, 0, 0, 0, 0.2e1 * t134 * qJD(4), t77 + (-t61 - t122) * qJD(4), -t58 - t168, -t105 * t134 + t20 * t61 + t76, 0, 0, 0, 0, 0, 0, t109 - t157, -t57 ^ 2 * t96 - t152 - t155, (t27 - t158) * t96 + t169 + t143, t114 * t93 - t134 * t15 + t172 * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, -t58 + t168, t77 + (t61 - t122) * qJD(4), -t151, 0, 0, -t134 * t56 - t173, t56 * t61 - t171, 0, 0, t116 * t48 - t160, (-t27 - t158) * t96 - t169 + t143, t116 * t57 + t152 - t155, t117 * t46 - t159, t109 + t157, -t57 * t134, -pkin(4) * t28 - t12 * t96 - t20 * t46 + t134 * t111 + (-pkin(8) * t129 - t9) * t57 + t101 * t93, pkin(4) * t27 + t12 * t93 - t20 * t48 + t134 * t8 + (pkin(8) * t130 + t10) * t57 + t101 * t96, t10 * t46 + t48 * t9 + ((-t28 + t131) * pkin(8) + t114) * t96 + ((qJD(5) * t46 - t27) * pkin(8) - t172) * t93, -pkin(4) * t12 - t10 * t8 - t15 * t20 + t111 * t9 + (-qJD(5) * t115 + t1 * t96 - t2 * t93) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t46 ^ 2 + t48 ^ 2, t46 * t57 - t27, -t156, t48 * t57 - t28, t55, -t15 * t48 + t172, t15 * t46 - t114, 0, 0;];
tauc_reg = t3;
