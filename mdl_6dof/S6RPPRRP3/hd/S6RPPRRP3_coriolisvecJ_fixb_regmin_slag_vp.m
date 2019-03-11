% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:03:48
% EndTime: 2019-03-09 02:03:53
% DurationCPUTime: 1.40s
% Computational Cost: add. (2011->260), mult. (4004->346), div. (0->0), fcn. (2249->6), ass. (0->140)
t79 = sin(qJ(4));
t136 = qJD(1) * t79;
t67 = qJD(5) + t136;
t104 = qJD(5) * t79 + qJD(1);
t81 = cos(qJ(4));
t135 = qJD(1) * t81;
t80 = cos(qJ(5));
t115 = t80 * t135;
t125 = t80 * qJD(4);
t127 = qJD(5) * t81;
t78 = sin(qJ(5));
t112 = t78 * t127;
t87 = t79 * t125 + t112;
t33 = t87 * qJD(1) - qJD(5) * t125;
t147 = t81 * t33;
t148 = t80 * t81;
t154 = t67 * t78;
t134 = qJD(4) * t78;
t59 = t115 + t134;
t176 = qJD(4) * ((-t59 + t115) * t79 + t67 * t148) - t104 * t154 - t147;
t124 = qJD(1) * qJD(4);
t108 = t81 * t124;
t103 = pkin(5) * t108;
t128 = qJD(5) * t80;
t129 = qJD(5) * t78;
t66 = -cos(pkin(9)) * pkin(1) - pkin(2) - pkin(7);
t54 = t66 * qJD(1) + qJD(3);
t73 = t81 * qJD(2);
t36 = t79 * t54 + t73;
t29 = qJD(4) * pkin(8) + t36;
t173 = -t79 * qJD(2) + t54 * t81;
t30 = t173 * qJD(4);
t68 = sin(pkin(9)) * pkin(1) + qJ(3);
t51 = pkin(4) * t79 - pkin(8) * t81 + t68;
t39 = t51 * qJD(1);
t101 = pkin(4) * t81 + pkin(8) * t79;
t55 = t101 * qJD(4) + qJD(3);
t42 = t55 * qJD(1);
t107 = t29 * t128 + t39 * t129 + t78 * t30 - t80 * t42;
t2 = -t103 + t107;
t10 = t29 * t80 + t39 * t78;
t7 = qJ(6) * t67 + t10;
t175 = -t67 * t7 + t2;
t132 = qJD(4) * t81;
t130 = qJD(5) * t59;
t109 = t78 * t124;
t64 = t79 * t109;
t34 = -t64 + t130;
t23 = t79 * t34;
t57 = t78 * t135 - t125;
t174 = t57 * t132 + t23;
t24 = t79 * t33;
t142 = t59 * t132 - t24;
t151 = t79 * t80;
t141 = t66 * t151 + t78 * t51;
t75 = t81 ^ 2;
t65 = t80 * t75 * t124;
t172 = -t67 * t87 + t65;
t133 = qJD(4) * t79;
t31 = qJD(4) * t73 + t54 * t133;
t3 = pkin(5) * t34 + qJ(6) * t33 - qJD(6) * t59 + t31;
t9 = -t29 * t78 + t39 * t80;
t138 = qJD(6) - t9;
t6 = -pkin(5) * t67 + t138;
t98 = t6 * t78 + t7 * t80;
t171 = qJD(4) * t98 - t3;
t168 = t59 ^ 2;
t167 = 0.2e1 * qJD(3);
t166 = t3 * t78;
t165 = t3 * t80;
t28 = -qJD(4) * pkin(4) - t173;
t8 = pkin(5) * t57 - qJ(6) * t59 + t28;
t164 = t59 * t8;
t162 = t28 * t78;
t161 = t28 * t80;
t160 = t31 * t78;
t159 = t31 * t80;
t157 = t57 * t67;
t156 = t59 * t57;
t155 = t59 * t67;
t153 = t67 * t80;
t152 = t78 * t79;
t150 = t80 * t51;
t149 = t80 * t55;
t82 = qJD(4) ^ 2;
t146 = t82 * t79;
t145 = t82 * t81;
t96 = pkin(5) * t78 - qJ(6) * t80;
t144 = t78 * qJD(6) - t67 * t96 + t36;
t61 = t101 * qJD(1);
t143 = t173 * t80 + t78 * t61;
t140 = t79 ^ 2 - t75;
t83 = qJD(1) ^ 2;
t139 = -t82 - t83;
t63 = qJD(1) * t68;
t137 = qJD(1) * t75;
t131 = qJD(5) * t57;
t126 = t63 * qJD(1);
t123 = pkin(8) * t154;
t122 = pkin(8) * t153;
t121 = t67 * t152;
t120 = t67 * t151;
t118 = -t39 * t128 - t80 * t30 - t78 * t42;
t117 = t81 * t66 * t125 + t51 * t128 + t78 * t55;
t116 = pkin(8) * t132;
t45 = t57 * t133;
t114 = t59 * t133;
t113 = t59 * t127;
t111 = t67 * t128;
t110 = t66 * t78 - pkin(5);
t106 = -t33 + t131;
t105 = t81 * t111;
t102 = qJ(6) * t108;
t99 = t6 * t80 - t7 * t78;
t97 = pkin(5) * t80 + qJ(6) * t78;
t95 = -t173 * t78 + t61 * t80;
t93 = -t79 * t8 + t116;
t92 = -t66 + t96;
t91 = t28 * t79 - t116;
t90 = t10 * t67 - t107;
t89 = -t105 + t174;
t88 = t29 * t129 + t118;
t1 = qJD(6) * t67 + t102 - t88;
t86 = t99 * qJD(5) + t1 * t80 + t2 * t78;
t85 = qJD(4) * t8 + t86;
t84 = t45 + (-t34 - t64) * t81 + (-t104 * t80 - t78 * t132) * t67;
t62 = -pkin(4) - t97;
t53 = qJD(4) * t121;
t32 = t92 * t81;
t20 = pkin(5) * t59 + qJ(6) * t57;
t19 = t34 * t148;
t17 = t110 * t79 - t150;
t16 = qJ(6) * t79 + t141;
t14 = t157 - t33;
t13 = -pkin(5) * t135 - t95;
t12 = qJ(6) * t135 + t143;
t11 = (t97 * qJD(5) - qJD(6) * t80) * t81 - t92 * t133;
t5 = qJD(5) * t141 + t110 * t132 - t149;
t4 = qJ(6) * t132 + (-t66 * t129 + qJD(6)) * t79 + t117;
t15 = [0, 0, 0, 0, 0, qJD(1) * t167, t63 * t167, -0.2e1 * t79 * t108, 0.2e1 * t140 * t124, -t146, -t145, 0, t63 * t132 - t66 * t146 + (t132 * t68 + t79 * t167) * qJD(1), -t63 * t133 - t66 * t145 + (-t133 * t68 + t81 * t167) * qJD(1), -t80 * t147 - t59 * t87, -t19 + (-t113 + t45) * t80 + (t114 + (t33 + t131) * t81) * t78, t142 + t172, -t105 - t23 + t53 + (-t137 * t78 - t57 * t81) * qJD(4) (t67 + t136) * t132 (-t129 * t51 + t149) * t67 + (-t66 * t111 + (t57 * t66 - t162) * qJD(4) - t107) * t79 + (t28 * t128 + t160 - t66 * t34 + (-t66 * t154 + (-t66 * t152 + t150) * qJD(1) + t9) * qJD(4)) * t81, -t117 * t67 + ((t66 * t67 + t29) * t129 + (t59 * t66 - t161) * qJD(4) + t118) * t79 + (-t28 * t129 + t159 + t66 * t33 + (-t141 * qJD(1) - t10) * qJD(4)) * t81, t11 * t57 + t32 * t34 - t5 * t67 + (-t134 * t8 - t2) * t79 + (t8 * t128 + t166 + (-qJD(1) * t17 - t6) * qJD(4)) * t81, -t16 * t34 - t17 * t33 - t4 * t57 + t5 * t59 - t99 * t133 + (-qJD(5) * t98 - t1 * t78 + t2 * t80) * t81, -t11 * t59 + t32 * t33 + t4 * t67 + (t125 * t8 + t1) * t79 + (t8 * t129 - t165 + (qJD(1) * t16 + t7) * qJD(4)) * t81, t1 * t16 + t11 * t8 + t17 * t2 + t3 * t32 + t4 * t7 + t5 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t145, t146, 0, 0, 0, 0, 0, -t109 * t75 + t53 + t89, t142 - t172 (t67 * t79 - t137) * t134 + t89, -t19 + (t113 + t45) * t80 + (t106 * t81 - t114) * t78, -t67 * t112 + t24 + t65 + (-t59 * t81 - t120) * qJD(4), -t171 * t79 + t85 * t81; 0, 0, 0, 0, 0, -t83, -t126, 0, 0, 0, 0, 0, t139 * t79, t139 * t81, 0, 0, 0, 0, 0, t84, -t176, t84 (t104 * t59 - t174) * t80 + (t104 * t57 + t142) * t78, t176, t99 * qJD(1) + t171 * t81 + t85 * t79; 0, 0, 0, 0, 0, 0, 0, t81 * t83 * t79, -t140 * t83, 0, 0, 0, qJD(4) * t36 - t126 * t81 - t31, t79 * t126, t59 * t153 - t33 * t78 (-t33 - t157) * t80 + (-t34 - t155) * t78, t111 + (t120 + (-t59 + t134) * t81) * qJD(1), -t67 * t129 + (-t121 + (t57 + t125) * t81) * qJD(1), -t67 * t135, -pkin(4) * t34 - t159 - t95 * t67 - t36 * t57 + (-t122 + t162) * qJD(5) + (t78 * t91 - t9 * t81) * qJD(1), pkin(4) * t33 + t160 + t143 * t67 - t36 * t59 + (t123 + t161) * qJD(5) + (t10 * t81 + t80 * t91) * qJD(1), t13 * t67 - t165 + t62 * t34 - t144 * t57 + (t78 * t8 - t122) * qJD(5) + (t6 * t81 - t78 * t93) * qJD(1), t12 * t57 - t13 * t59 + (t1 + t67 * t6 + (-t34 + t130) * pkin(8)) * t80 + (pkin(8) * t106 + t175) * t78, -t12 * t67 - t166 + t62 * t33 + t144 * t59 + (-t8 * t80 - t123) * qJD(5) + (-t7 * t81 + t80 * t93) * qJD(1), t86 * pkin(8) - t7 * t12 - t6 * t13 - t144 * t8 + t3 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t156, -t57 ^ 2 + t168, t14, t155 - t34, t108, -t28 * t59 + t90, t28 * t57 + t67 * t9 + t88, -t20 * t57 + 0.2e1 * t103 - t164 + t90, pkin(5) * t33 - t34 * qJ(6) + (-t10 + t7) * t59 + (t6 - t138) * t57, 0.2e1 * t102 + t20 * t59 - t8 * t57 + (0.2e1 * qJD(6) - t9) * t67 - t88, -t2 * pkin(5) + t1 * qJ(6) - t6 * t10 + t138 * t7 - t8 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108 + t156, t14, -t67 ^ 2 - t168, t164 + t175;];
tauc_reg  = t15;
