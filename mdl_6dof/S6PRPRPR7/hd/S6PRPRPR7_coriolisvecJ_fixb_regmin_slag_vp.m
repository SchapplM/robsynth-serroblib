% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
% 
% Output:
% tauc_reg [6x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:53:47
% EndTime: 2019-03-08 19:53:52
% DurationCPUTime: 1.58s
% Computational Cost: add. (988->243), mult. (2250->335), div. (0->0), fcn. (1451->8), ass. (0->144)
t79 = cos(qJ(4));
t143 = t79 * qJD(2);
t64 = qJD(6) + t143;
t181 = qJD(6) - t64;
t74 = cos(pkin(6));
t153 = qJD(1) * t74;
t73 = sin(pkin(6));
t154 = qJD(1) * t73;
t80 = cos(qJ(2));
t126 = t80 * t154;
t105 = qJD(3) - t126;
t82 = -pkin(2) - pkin(8);
t34 = t82 * qJD(2) + t105;
t76 = sin(qJ(4));
t20 = t76 * t153 - t79 * t34;
t180 = qJD(5) + t20;
t75 = sin(qJ(6));
t145 = t75 * qJD(4);
t151 = qJD(2) * t76;
t78 = cos(qJ(6));
t179 = t78 * t151 - t145;
t139 = qJD(4) * qJ(5);
t149 = qJD(4) * t79;
t178 = pkin(4) * t149 + t76 * t139;
t101 = (t64 + t143) * t76;
t14 = -qJD(4) * pkin(4) + t180;
t21 = t79 * t153 + t76 * t34;
t15 = -t139 - t21;
t132 = pkin(5) * t143;
t141 = t132 + t180;
t138 = qJD(2) * qJD(4);
t122 = t79 * t138;
t144 = t78 * qJD(4);
t42 = t151 * t75 + t144;
t23 = qJD(6) * t42 - t78 * t122;
t152 = qJD(2) * t73;
t77 = sin(qJ(2));
t131 = t77 * t152;
t112 = t76 * t131;
t84 = qJD(2) ^ 2;
t163 = t73 * t84;
t135 = t80 * t163;
t164 = t73 * t80;
t32 = t164 * t79 + t74 * t76;
t16 = qJD(4) * t32 - t112;
t177 = qJD(4) * (-t16 + t112) - t79 * t135;
t113 = t79 * t131;
t33 = -t164 * t76 + t74 * t79;
t17 = qJD(4) * t33 - t113;
t176 = qJD(4) * (-t17 + t113) + t76 * t135;
t127 = t77 * t154;
t111 = qJD(2) * t127;
t123 = qJD(4) * t153;
t133 = -t34 * t149 + (-t111 + t123) * t76;
t6 = -qJD(4) * qJD(5) + t133;
t150 = qJD(4) * t76;
t159 = t79 * t123 + t34 * t150;
t7 = -t111 * t79 + t159;
t85 = -t6 * t76 - t7 * t79 + (t14 * t76 - t15 * t79) * qJD(4);
t81 = -pkin(4) - pkin(9);
t3 = (qJD(5) - t132) * qJD(4) - t133;
t175 = t3 * t75;
t174 = t3 * t78;
t173 = pkin(5) - t82;
t22 = t179 * qJD(6) + t75 * t122;
t172 = t22 * t78;
t171 = t179 * t64;
t170 = t179 * t76;
t169 = t42 * t64;
t140 = qJD(2) * qJ(3);
t44 = t127 + t140;
t168 = t44 * t80;
t167 = t64 * t79;
t166 = t64 * t81;
t165 = t73 * t77;
t162 = t76 * t22;
t161 = t78 * t79;
t83 = qJD(4) ^ 2;
t160 = t82 * t83;
t121 = t76 * t138;
t158 = pkin(4) * t122 + qJ(5) * t121;
t43 = pkin(4) * t143 + qJ(5) * t151;
t71 = t76 ^ 2;
t72 = t79 ^ 2;
t157 = t71 - t72;
t156 = t83 + t84;
t155 = qJD(2) * pkin(2);
t148 = qJD(6) * t78;
t13 = -pkin(5) * t151 + t21;
t10 = t13 + t139;
t147 = t10 * qJD(6);
t109 = pkin(4) * t151 + t127;
t117 = -t79 * qJ(5) + qJ(3);
t25 = qJD(2) * t117 + t109;
t146 = t25 * qJD(2);
t142 = t79 * qJD(5);
t137 = t75 * t167;
t136 = t77 * t163;
t134 = t79 * t84 * t76;
t130 = t80 * t152;
t129 = qJD(6) * t75 * t64;
t128 = t76 * t148;
t114 = t79 * t127;
t5 = (-pkin(5) * t150 - t114) * qJD(2) + t159;
t90 = qJD(3) + (qJD(4) * pkin(9) - qJD(5)) * t79;
t8 = (t90 + t126) * qJD(2) + t158;
t124 = t78 * t5 - t75 * t8;
t120 = qJD(4) * t173;
t110 = -t44 + t127;
t107 = t21 * qJD(4) - t159;
t92 = t76 * pkin(9) + t117;
t18 = qJD(2) * t92 + t109;
t9 = qJD(4) * t81 + t141;
t2 = t78 * t18 + t75 * t9;
t106 = t75 * t18 - t78 * t9;
t104 = qJD(3) + t126;
t70 = t76 * pkin(4);
t35 = t70 + t92;
t48 = t173 * t79;
t102 = t78 * t35 + t75 * t48;
t100 = -qJD(2) * t71 + t167;
t98 = t64 * (qJD(6) * t79 + qJD(2));
t94 = t165 * t78 + t32 * t75;
t93 = -t165 * t75 + t32 * t78;
t91 = t110 * qJD(2);
t46 = t117 + t70;
t89 = -qJD(2) * t46 + t127 - t25;
t88 = t110 - t140;
t11 = (t104 - t142) * qJD(2) + t158;
t30 = qJD(3) - t142 + t178;
t87 = t160 - t11 + (-t30 + t126) * qJD(2);
t38 = t104 * qJD(2);
t86 = qJD(2) * t105 - t160 + t38;
t55 = t75 * t121;
t51 = t156 * t79;
t50 = t156 * t76;
t47 = t173 * t76;
t39 = t105 - t155;
t37 = t79 * t120;
t36 = t76 * t120;
t31 = pkin(9) * t143 + t43;
t24 = t90 + t178;
t19 = t25 * t143;
t1 = [0, 0, -t136, -t135, t136, t135 (t38 * t77 + (t168 + (t39 - t126) * t77) * qJD(2)) * t73, 0, 0, 0, 0, 0, t176, -t177 (t16 * t76 + t17 * t79 + (-t32 * t76 - t33 * t79) * qJD(4)) * qJD(2), -t176, t177, t14 * t17 + t15 * t16 + t7 * t32 - t6 * t33 + (t11 * t77 + t146 * t80) * t73, 0, 0, 0, 0, 0 (-qJD(6) * t94 - t130 * t75 + t17 * t78) * t64 - t93 * t121 + t16 * t179 + t33 * t23 -(qJD(6) * t93 + t130 * t78 + t17 * t75) * t64 + t94 * t121 - t16 * t42 + t33 * t22; 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3), t38 * qJ(3) + t44 * qJD(3) + (-t168 + (-t39 - t155) * t77) * t154, -0.2e1 * t79 * t121, 0.2e1 * t157 * t138, -t83 * t76, -t83 * t79, 0, -t149 * t88 + t76 * t86, t150 * t88 + t79 * t86 (t71 + t72) * t111 - t85, t149 * t89 + t76 * t87, -t150 * t89 + t79 * t87, t11 * t46 + t25 * t30 + (-t25 * t80 + (t14 * t79 + t15 * t76) * t77) * t154 + t85 * t82, t75 * t162 + (t145 * t79 + t128) * t42 (t179 * t75 + t42 * t78) * t149 + (t172 - t23 * t75 + (t179 * t78 - t42 * t75) * qJD(6)) * t76, t64 * t128 + t22 * t79 + (t100 * t75 - t42 * t76) * qJD(4), -t76 * t129 - t23 * t79 + (t100 * t78 - t170) * qJD(4), -qJD(4) * t101, t37 * t179 - t47 * t23 + (t179 * t127 + t75 * t147 - t174 + (-(-t75 * t35 + t78 * t48) * qJD(2) + t106) * qJD(4)) * t76 + (-qJD(6) * t2 - t10 * t144 + t124) * t79 + (-t75 * t24 - t78 * t36 - t102 * qJD(6) - (-t161 * t77 - t75 * t80) * t154) * t64, -t47 * t22 - t37 * t42 + (-(qJD(6) * t9 + t8) * t79 + (-qJD(6) * t48 + t126 - t24) * t64) * t78 + (-(-qJD(6) * t35 - t36) * t64 + (t10 * qJD(4) + qJD(6) * t18 - t127 * t64 - t5) * t79) * t75 + (-t42 * t127 + t78 * t147 + t175 + (qJD(2) * t102 + t2) * qJD(4)) * t76; 0, 0, 0, 0, 0, -t84, t91, 0, 0, 0, 0, 0, -t50, -t51, 0, t50, t51, t85 - t146, 0, 0, 0, 0, 0, t76 * t23 + t75 * t98 + (t101 * t78 - t179 * t79) * qJD(4), t162 + t78 * t98 + (-t101 * t75 + t42 * t79) * qJD(4); 0, 0, 0, 0, 0, 0, 0, t134, -t157 * t84, 0, 0, 0, t79 * t91 + t107, -t20 * qJD(4) + t151 * t44 + t133, 0, t19 + (t43 * t76 - t114) * qJD(2) - t107 (0.2e1 * qJD(5) + t20) * qJD(4) + (-t25 * t76 + t43 * t79) * qJD(2) - t133, -t7 * pkin(4) - t6 * qJ(5) - t14 * t21 - t180 * t15 - t25 * t43, -t169 * t75 + t172 (-t23 - t169) * t78 + (-t22 - t171) * t75, -t129 + (-t137 + (t42 - t144) * t76) * qJD(2), -t64 * t148 + t55 + (-t161 * t64 + t170) * qJD(2), t64 * t151, qJ(5) * t23 + t175 - (t78 * t13 - t75 * t31) * t64 - t141 * t179 + (t10 * t78 - t166 * t75) * qJD(6) + (t10 * t161 + (-t144 * t81 - t106) * t76) * qJD(2), qJ(5) * t22 + t174 + (t75 * t13 + t78 * t31) * t64 + t141 * t42 + (-t10 * t75 - t166 * t78) * qJD(6) + (-t2 * t76 + (-t10 * t79 + t150 * t81) * t75) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t134, -t72 * t84 - t83, t15 * qJD(4) + t19 + t7, 0, 0, 0, 0, 0, -t129 + qJD(4) * t179 + (-t144 * t76 - t137) * qJD(2), -t64 ^ 2 * t78 - qJD(4) * t42 + t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t42 * t179, -t179 ^ 2 + t42 ^ 2, t22 - t171, t169 - t23, -t121, -t10 * t42 - t181 * t2 + t124, -t10 * t179 + t181 * t106 - t75 * t5 - t78 * t8;];
tauc_reg  = t1;
