% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRRRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:09
% EndTime: 2019-03-08 18:58:16
% DurationCPUTime: 2.40s
% Computational Cost: add. (3190->312), mult. (8447->449), div. (0->0), fcn. (7128->12), ass. (0->166)
t104 = sin(pkin(12));
t106 = sin(pkin(6));
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t107 = cos(pkin(12));
t108 = cos(pkin(7));
t184 = t107 * t108;
t123 = (t104 * t115 + t112 * t184) * t106;
t105 = sin(pkin(7));
t187 = t105 * t112;
t109 = cos(pkin(6));
t96 = qJD(1) * t109 + qJD(2);
t57 = qJD(1) * t123 + t96 * t187;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t141 = pkin(4) * t111 - pkin(10) * t114;
t91 = t141 * qJD(4);
t222 = t57 - t91;
t186 = t105 * t115;
t214 = (-t104 * t112 + t115 * t184) * t106;
t221 = t109 * t186 + t214;
t169 = t114 * qJD(3);
t97 = -qJD(5) + t169;
t220 = qJD(1) * t214 + t96 * t186;
t168 = qJD(3) * qJD(4);
t148 = t111 * t168;
t146 = pkin(5) * t148;
t110 = sin(qJ(5));
t113 = cos(qJ(5));
t54 = qJD(3) * pkin(9) + t57;
t180 = qJD(1) * t106;
t155 = t107 * t180;
t73 = -t105 * t155 + t108 * t96;
t216 = -t111 * t54 + t114 * t73;
t50 = t220 * qJD(3);
t12 = qJD(4) * t216 + t114 * t50;
t173 = qJD(5) * t113;
t174 = qJD(5) * t110;
t33 = t111 * t73 + t114 * t54;
t31 = qJD(4) * pkin(10) + t33;
t145 = t108 * t155;
t178 = qJD(3) * t112;
t154 = t105 * t178;
t156 = t104 * t180;
t177 = qJD(3) * t115;
t51 = t145 * t178 + t96 * t154 + t156 * t177;
t39 = qJD(3) * t91 + t51;
t56 = -t112 * t156 + t115 * (t105 * t96 + t145);
t94 = -t114 * pkin(4) - t111 * pkin(10) - pkin(3);
t42 = qJD(3) * t94 - t56;
t147 = -t110 * t12 + t113 * t39 - t31 * t173 - t42 * t174;
t2 = -t146 - t147;
t9 = t110 * t42 + t113 * t31;
t7 = -t97 * qJ(6) + t9;
t219 = t7 * t97 + t2;
t183 = t110 * t114;
t218 = t222 * t113 + t94 * t174 - t56 * t183;
t182 = t113 * t114;
t217 = -t222 * t110 + t94 * t173 - t56 * t182;
t215 = qJD(3) * t57 - t51;
t151 = t111 * t173;
t167 = qJD(4) * qJD(5);
t171 = t110 * qJD(4);
t72 = qJD(3) * (t114 * t171 + t151) + t110 * t167;
t179 = qJD(3) * t111;
t87 = t113 * t179 + t171;
t213 = t87 ^ 2;
t30 = -qJD(4) * pkin(4) - t216;
t150 = t110 * t179;
t170 = t113 * qJD(4);
t85 = t150 - t170;
t15 = t85 * pkin(5) - t87 * qJ(6) + t30;
t211 = t15 * t87;
t175 = qJD(4) * t114;
t176 = qJD(4) * t111;
t13 = t111 * t50 + t54 * t175 + t73 * t176;
t152 = t114 * t170;
t71 = -qJD(3) * t152 + qJD(5) * t150 - t113 * t167;
t5 = pkin(5) * t72 + qJ(6) * t71 - qJD(6) * t87 + t13;
t210 = t5 * t110;
t209 = t5 * t113;
t208 = t56 * t85;
t207 = t56 * t87;
t206 = t85 * t97;
t205 = t87 * t85;
t204 = t87 * t97;
t138 = pkin(5) * t110 - qJ(6) * t113;
t203 = t110 * qJD(6) + t97 * t138 + t33;
t172 = qJD(5) * t114;
t202 = -pkin(5) * t176 + (-t111 * t171 + t113 * t172) * pkin(9) + t218;
t201 = -qJ(6) * t176 + t114 * qJD(6) - (-t110 * t172 - t111 * t170) * pkin(9) - t217;
t90 = t141 * qJD(3);
t200 = t110 * t90 + t113 * t216;
t198 = pkin(9) * t182 + t110 * t94;
t197 = qJD(3) * pkin(3);
t196 = t110 * t97;
t195 = t113 * t94;
t194 = t113 * t97;
t192 = t13 * t110;
t191 = t13 * t113;
t190 = t71 * t110;
t8 = -t110 * t31 + t113 * t42;
t189 = qJD(6) - t8;
t117 = qJD(3) ^ 2;
t185 = t105 * t117;
t102 = t111 ^ 2;
t181 = -t114 ^ 2 + t102;
t166 = pkin(10) * t196;
t165 = pkin(10) * t194;
t164 = pkin(10) * t176;
t163 = pkin(10) * t170;
t161 = t97 * t174;
t160 = t97 * t173;
t158 = t111 * t187;
t157 = t114 * t187;
t153 = t105 * t177;
t144 = t114 * t153;
t143 = t111 * t153;
t142 = qJ(6) * t148;
t6 = t97 * pkin(5) + t189;
t140 = -t110 * t7 + t113 * t6;
t139 = t113 * pkin(5) + t110 * qJ(6);
t137 = -t110 * t216 + t113 * t90;
t61 = t109 * t187 + t123;
t76 = -t105 * t106 * t107 + t108 * t109;
t41 = t76 * t111 + t61 * t114;
t21 = -t110 * t221 + t41 * t113;
t20 = t41 * t110 + t113 * t221;
t40 = t61 * t111 - t76 * t114;
t135 = qJD(3) * t102 - t114 * t97;
t134 = pkin(9) + t138;
t133 = -t9 * t97 + t147;
t116 = qJD(4) ^ 2;
t132 = pkin(9) * t116 - t215;
t53 = -t56 - t197;
t131 = qJD(4) * (t53 + t56 - t197);
t78 = t111 * t108 + t157;
t65 = t110 * t78 + t113 * t186;
t66 = -t110 * t186 + t113 * t78;
t77 = -t114 * t108 + t158;
t128 = -t110 * t39 - t113 * t12 - t42 * t173 + t31 * t174;
t58 = t221 * qJD(3);
t16 = qJD(4) * t41 + t58 * t111;
t17 = -qJD(4) * t40 + t58 * t114;
t59 = t61 * qJD(3);
t3 = qJD(5) * t21 + t17 * t110 - t59 * t113;
t122 = -t148 * t20 + t16 * t85 + t3 * t97 + t40 * t72;
t63 = -qJD(4) * t77 + t144;
t29 = qJD(5) * t66 + t110 * t63 - t113 * t154;
t64 = qJD(4) * t78 + t143;
t121 = -t148 * t65 + t29 * t97 + t64 * t85 + t77 * t72;
t4 = -qJD(5) * t20 + t59 * t110 + t17 * t113;
t119 = t148 * t21 - t16 * t87 - t4 * t97 + t40 * t71;
t28 = -qJD(5) * t65 + t110 * t154 + t113 * t63;
t118 = t148 * t66 - t28 * t97 - t64 * t87 + t77 * t71;
t93 = -pkin(4) - t139;
t74 = t134 * t111;
t70 = -t195 + (pkin(9) * t110 + pkin(5)) * t114;
t69 = -t114 * qJ(6) + t198;
t62 = pkin(5) * t87 + qJ(6) * t85;
t49 = -t71 - t206;
t44 = (qJD(5) * t139 - qJD(6) * t113) * t111 + t134 * t175;
t19 = -pkin(5) * t179 - t137;
t18 = qJ(6) * t179 + t200;
t1 = -t97 * qJD(6) - t128 + t142;
t10 = [0, 0, 0, -t59 * qJD(3), -t58 * qJD(3), 0, 0, 0, 0, 0, -t16 * qJD(4) + (-t114 * t59 - t176 * t221) * qJD(3), -t17 * qJD(4) + (t111 * t59 - t175 * t221) * qJD(3), 0, 0, 0, 0, 0, t122, -t119, t122, -t20 * t71 - t21 * t72 + t3 * t87 - t4 * t85, t119, t1 * t21 + t15 * t16 + t2 * t20 + t3 * t6 + t4 * t7 + t40 * t5; 0, 0, 0, -t112 * t185, -t115 * t185, 0, 0, 0, 0, 0, -t117 * t157 + (-t64 - t143) * qJD(4), t117 * t158 + (-t63 - t144) * qJD(4), 0, 0, 0, 0, 0, t121, -t118, t121, -t28 * t85 + t29 * t87 - t65 * t71 - t66 * t72, t118, t1 * t66 + t15 * t64 + t2 * t65 + t28 * t7 + t29 * t6 + t5 * t77; 0, 0, 0, t215 (-t220 + t56) * qJD(3), 0.2e1 * t114 * t148, -0.2e1 * t181 * t168, t116 * t114, -t116 * t111, 0, t111 * t131 - t114 * t132, t111 * t132 + t114 * t131, t87 * t152 + (-t71 * t113 - t87 * t174) * t111 (-t110 * t87 - t113 * t85) * t175 + (t190 - t113 * t72 + (t110 * t85 - t113 * t87) * qJD(5)) * t111, t111 * t161 + t71 * t114 + (t111 * t87 + t113 * t135) * qJD(4), t97 * t151 + t72 * t114 + (-t110 * t135 - t111 * t85) * qJD(4) (-t97 - t169) * t176, t218 * t97 + (t30 * t171 + (qJD(4) * t85 + t160) * pkin(9) - t147) * t114 + (t30 * t173 + pkin(9) * t72 + t192 - t208 + (-pkin(9) * t196 + (-pkin(9) * t183 + t195) * qJD(3) + t8) * qJD(4)) * t111, t217 * t97 + (t30 * t170 + (qJD(4) * t87 - t161) * pkin(9) - t128) * t114 + (-t30 * t174 - pkin(9) * t71 + t191 - t207 + (-pkin(9) * t194 - qJD(3) * t198 - t9) * qJD(4)) * t111, t44 * t85 + t74 * t72 + t202 * t97 + (t15 * t171 + t2) * t114 + (t15 * t173 + t210 - t208 + (-qJD(3) * t70 - t6) * qJD(4)) * t111, -t69 * t72 - t70 * t71 + t202 * t87 + t201 * t85 + t140 * t175 + (-t1 * t110 + t113 * t2 + (-t110 * t6 - t113 * t7) * qJD(5)) * t111, -t44 * t87 + t74 * t71 + t201 * t97 + (-t15 * t170 - t1) * t114 + (t15 * t174 - t209 + t207 + (qJD(3) * t69 + t7) * qJD(4)) * t111, t1 * t69 + t2 * t70 + t5 * t74 - t201 * t7 + t202 * t6 + (-t111 * t56 + t44) * t15; 0, 0, 0, 0, 0, -t111 * t117 * t114, t181 * t117, 0, 0, 0, t33 * qJD(4) - t53 * t179 - t13 (-qJD(3) * t53 - t50) * t114, -t87 * t194 - t190 (-t71 + t206) * t113 + (-t72 + t204) * t110, -t160 + (t97 * t182 + (-t87 + t171) * t111) * qJD(3), t161 + (-t97 * t183 + (t85 + t170) * t111) * qJD(3), t97 * t179, -pkin(4) * t72 - t191 + t137 * t97 - t33 * t85 + (t30 * t110 + t165) * qJD(5) + (-t8 * t111 + (-t114 * t30 - t164) * t110) * qJD(3), pkin(4) * t71 + t192 - t200 * t97 - t33 * t87 + (t30 * t113 - t166) * qJD(5) + (-t30 * t182 + (t9 - t163) * t111) * qJD(3), -t209 - t19 * t97 + t93 * t72 - t203 * t85 + (t110 * t15 + t165) * qJD(5) + (t111 * t6 + (-t114 * t15 - t164) * t110) * qJD(3), t18 * t85 - t19 * t87 + (t1 - t97 * t6 + (qJD(5) * t87 - t72) * pkin(10)) * t113 + ((qJD(5) * t85 - t71) * pkin(10) + t219) * t110, -t210 + t18 * t97 + t93 * t71 + t203 * t87 + (-t113 * t15 + t166) * qJD(5) + (t15 * t182 + (-t7 + t163) * t111) * qJD(3), -t7 * t18 - t6 * t19 + t5 * t93 - t203 * t15 + (qJD(5) * t140 + t1 * t113 + t2 * t110) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t205, -t85 ^ 2 + t213, t49, -t204 - t72, t148, -t30 * t87 + t133, t30 * t85 - t8 * t97 + t128, -t62 * t85 + t133 + 0.2e1 * t146 - t211, pkin(5) * t71 - qJ(6) * t72 + (t7 - t9) * t87 + (t6 - t189) * t85, 0.2e1 * t142 - t15 * t85 + t62 * t87 + (-0.2e1 * qJD(6) + t8) * t97 - t128, -pkin(5) * t2 + qJ(6) * t1 - t15 * t62 + t189 * t7 - t6 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148 + t205, t49, -t97 ^ 2 - t213, t211 + t219;];
tauc_reg  = t10;
