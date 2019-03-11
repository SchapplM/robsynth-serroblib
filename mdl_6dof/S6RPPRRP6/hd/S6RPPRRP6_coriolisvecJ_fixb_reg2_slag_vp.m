% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:11:18
% EndTime: 2019-03-09 02:11:24
% DurationCPUTime: 2.38s
% Computational Cost: add. (2570->322), mult. (5014->417), div. (0->0), fcn. (2689->4), ass. (0->174)
t92 = sin(qJ(4));
t175 = qJD(1) * t92;
t221 = -t175 - qJD(5);
t93 = cos(qJ(5));
t118 = t221 * t93;
t165 = t93 * qJD(4);
t94 = cos(qJ(4));
t174 = qJD(1) * t94;
t91 = sin(qJ(5));
t63 = t91 * t174 - t165;
t111 = t63 * t118;
t153 = t93 * t174;
t172 = qJD(4) * t91;
t65 = t153 + t172;
t191 = t65 * t221;
t166 = qJD(5) * t94;
t150 = t91 * t166;
t107 = t92 * t165 + t150;
t32 = qJD(1) * t107 - qJD(5) * t165;
t31 = t93 * t32;
t161 = qJD(1) * qJD(4);
t145 = t92 * t161;
t74 = t91 * t145;
t33 = qJD(5) * t65 - t74;
t224 = (t33 - t191) * t91 - t111 + t31;
t216 = t33 + t191;
t167 = qJD(5) * t93;
t169 = qJD(5) * t91;
t162 = qJD(1) * qJD(2);
t170 = qJD(4) * t94;
t82 = qJD(1) * qJ(2) + qJD(3);
t76 = -pkin(7) * qJD(1) + t82;
t47 = t92 * t162 + t76 * t170;
t90 = pkin(1) + qJ(3);
t70 = pkin(4) * t92 - pkin(8) * t94 + t90;
t48 = qJD(1) * t70 - qJD(2);
t136 = pkin(4) * t94 + pkin(8) * t92;
t61 = t136 * qJD(4) + qJD(3);
t49 = t61 * qJD(1);
t68 = t92 * t76;
t55 = qJD(4) * pkin(8) + t68;
t144 = t55 * t167 + t48 * t169 + t91 * t47 - t93 * t49;
t20 = t48 * t91 + t55 * t93;
t112 = -t20 * t221 - t144;
t171 = qJD(4) * t92;
t190 = t65 * t93;
t198 = t33 * t93;
t200 = t32 * t91;
t220 = ((t63 * t91 - t190) * qJD(5) - t198 + t200) * t94 + (t63 * t93 + t65 * t91) * t171;
t88 = t94 ^ 2;
t121 = qJD(1) * t88 + t221 * t92;
t149 = t93 * t166;
t219 = (t121 * t91 + t63 * t94) * qJD(4) - t221 * t149 + t92 * t33;
t168 = qJD(5) * t92;
t140 = qJD(1) + t168;
t180 = t94 * t32;
t181 = t93 * t94;
t188 = t221 * t91;
t218 = ((-t65 + t153) * t92 - t221 * t181) * qJD(4) + t140 * t188 - t180;
t16 = -qJ(6) * t221 + t20;
t80 = t94 * t161;
t139 = pkin(5) * t80;
t2 = -t139 + t144;
t217 = t16 * t221 + t2;
t119 = t63 * t221;
t17 = -t32 - t119;
t215 = qJD(1) * t90;
t87 = t92 ^ 2;
t177 = t87 + t88;
t213 = t177 * qJD(2);
t133 = -t167 * t221 + t91 * t80;
t184 = t92 * t93;
t157 = t221 * t184;
t189 = t65 * t94;
t209 = (-t157 + t189) * qJD(1) + t133;
t89 = -pkin(7) + qJ(2);
t151 = t89 * t168;
t173 = qJD(2) * t92;
t12 = -(qJD(5) * t70 + t89 * t170 + t173) * t91 - (-t61 + t151) * t93;
t207 = t65 ^ 2;
t206 = pkin(8) * t65;
t46 = -t94 * t162 + t76 * t171;
t5 = pkin(5) * t33 + qJ(6) * t32 - qJD(6) * t65 + t46;
t205 = t5 * t91;
t204 = t5 * t93;
t56 = -qJD(4) * pkin(4) - t76 * t94;
t18 = pkin(5) * t63 - qJ(6) * t65 + t56;
t202 = t18 * t65;
t199 = t32 * t92;
t197 = t46 * t91;
t196 = t46 * t93;
t195 = t56 * t91;
t194 = t56 * t93;
t192 = t65 * t63;
t96 = qJD(1) ^ 2;
t187 = t87 * t96;
t186 = t91 * t92;
t185 = t91 * t94;
t67 = t136 * qJD(1);
t183 = t93 * t67;
t182 = t93 * t70;
t131 = pkin(5) * t91 - qJ(6) * t93;
t179 = t91 * qJD(6) + t221 * t131 + t68;
t24 = t76 * t181 + t91 * t67;
t158 = t221 * t186;
t69 = t221 * t169;
t178 = qJD(1) * t158 + t69;
t37 = t89 * t184 + t91 * t70;
t95 = qJD(4) ^ 2;
t176 = -t95 - t96;
t77 = -qJD(2) + t215;
t164 = qJD(2) - t77;
t19 = t48 * t93 - t55 * t91;
t163 = qJD(6) - t19;
t160 = pkin(8) * t188;
t159 = pkin(8) * t118;
t155 = t94 * t96 * t92;
t154 = pkin(8) * t170;
t152 = t94 * t165;
t148 = t221 * t174;
t85 = 0.2e1 * t162;
t147 = 0.2e1 * qJD(3) * qJD(1);
t146 = t63 ^ 2 - t207;
t143 = t89 * t152 + t70 * t167 + t93 * t173 + t91 * t61;
t142 = t77 + t215;
t141 = t164 * qJD(1);
t138 = t92 * t80;
t137 = qJ(6) * t80;
t132 = pkin(5) * t93 + qJ(6) * t91;
t15 = pkin(5) * t221 + t163;
t130 = t15 * t93 - t16 * t91;
t129 = t15 * t91 + t16 * t93;
t128 = t19 * t93 + t20 * t91;
t127 = t19 * t91 - t20 * t93;
t123 = qJD(2) + t142;
t122 = (qJD(5) * t63 - t32) * pkin(8);
t117 = -t89 * t95 + t147;
t115 = t131 - t89;
t114 = -t18 * t92 + t154;
t113 = t56 * t92 - t154;
t110 = -t48 * t167 + t55 * t169 - t93 * t47 - t91 * t49;
t109 = -t93 * t80 - t178;
t106 = t216 * t91;
t104 = -t91 * t119 - t198;
t1 = -qJD(6) * t221 - t110 + t137;
t103 = -t129 * qJD(5) - t1 * t91 + t2 * t93;
t102 = t130 * qJD(5) + t1 * t93 + t2 * t91;
t101 = t127 * qJD(5) + t110 * t91 + t144 * t93;
t100 = -t128 * qJD(5) - t110 * t93 + t144 * t91;
t99 = t33 * t185 + (-t91 * t171 + t149) * t63;
t98 = t63 * t171 + (-t33 - t74) * t94 - (-t140 * t93 - t91 * t170) * t221;
t97 = -t33 * t184 - t63 * t152 + t140 * t190 + (t140 * t63 + t65 * t170 - t199) * t91;
t84 = t88 * t96;
t71 = -pkin(4) - t132;
t50 = t63 * t174;
t43 = t115 * t94;
t42 = (-t221 + t175) * t170;
t36 = -t89 * t186 + t182;
t30 = -t182 + (t89 * t91 - pkin(5)) * t92;
t29 = pkin(8) * t198;
t28 = qJ(6) * t92 + t37;
t26 = pkin(5) * t65 + qJ(6) * t63;
t23 = -t76 * t185 + t183;
t22 = -t183 + (-pkin(5) * qJD(1) + t76 * t91) * t94;
t21 = qJ(6) * t174 + t24;
t14 = (-t157 - t189) * qJD(1) + t133;
t13 = -t115 * t171 + (t132 * qJD(5) - qJD(6) * t93 - qJD(2)) * t94;
t11 = -t91 * t151 + t143;
t10 = -t65 * t118 - t200;
t9 = -pkin(5) * t170 - t12;
t8 = -t107 * t65 - t93 * t180;
t7 = qJ(6) * t170 + (-t89 * t169 + qJD(6)) * t92 + t143;
t6 = t221 * t150 - t199 + (t121 * t93 + t189) * qJD(4);
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, qJ(2) * t85, 0, 0, 0, 0, 0, 0, 0, t85, t147, t82 * qJD(2) + t77 * qJD(3) + (qJ(2) * qJD(2) + qJD(3) * t90) * qJD(1), -0.2e1 * t138, 0.2e1 * (t87 - t88) * t161, -t95 * t92, 0.2e1 * t138, -t95 * t94, 0, t117 * t92 + t123 * t170, t117 * t94 - t123 * t171, -t177 * t85, t142 * qJD(3) + (qJD(1) * t89 + t76) * t213, t8, t220, t6, t99, -t219, t42, -t12 * t221 + (-t144 + (t63 * t89 - t195) * qJD(4)) * t92 + (t56 * t167 - qJD(2) * t63 - t33 * t89 + t197 + (qJD(1) * t36 + t19) * qJD(4)) * t94, t11 * t221 + (t110 + (t65 * t89 - t194) * qJD(4)) * t92 + (-t56 * t169 - qJD(2) * t65 + t32 * t89 + t196 + (-qJD(1) * t37 - t20) * qJD(4)) * t94, t101 * t94 - t11 * t63 - t12 * t65 + t128 * t171 + t32 * t36 - t33 * t37, -t46 * t89 * t94 + t11 * t20 + t12 * t19 - t110 * t37 - t36 * t144 + (-qJD(2) * t94 + t89 * t171) * t56, t8, t6, -t220, t42, t219, t99, t13 * t63 + t33 * t43 + t221 * t9 + (-t18 * t172 - t2) * t92 + (t18 * t167 + t205 + (-qJD(1) * t30 - t15) * qJD(4)) * t94, t103 * t94 - t130 * t171 - t28 * t33 - t30 * t32 - t63 * t7 + t65 * t9, -t13 * t65 + t32 * t43 - t7 * t221 + (t18 * t165 + t1) * t92 + (t18 * t169 - t204 + (qJD(1) * t28 + t16) * qJD(4)) * t94, t1 * t28 + t13 * t18 + t15 * t9 + t16 * t7 + t2 * t30 + t43 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t96 * qJ(2), 0, 0, 0, 0, 0, 0, 0, -t96, 0 (-qJD(3) - t82) * qJD(1), 0, 0, 0, 0, 0, 0, -0.2e1 * t80, 0.2e1 * t145, t84 + t187 (-t177 * t76 - qJD(3)) * qJD(1), 0, 0, 0, 0, 0, 0, t109 + t50, t209, t17 * t93 + t106 (t127 * t92 + t56 * t94) * qJD(1) + t101, 0, 0, 0, 0, 0, 0, -t69 + t50 + (-t152 - t158) * qJD(1), -t31 - t111 + t106, -t209 (-t129 * t92 + t18 * t94) * qJD(1) + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, t141, 0, 0, 0, 0, 0, 0, t176 * t92, t176 * t94, 0 (-t77 + t213) * qJD(1), 0, 0, 0, 0, 0, 0, t98, -t218, t97, -t128 * qJD(1) + (-t127 * qJD(4) - t46) * t94 + (qJD(4) * t56 + t100) * t92, 0, 0, 0, 0, 0, 0, t98, t97, t218, t130 * qJD(1) + (qJD(4) * t129 - t5) * t94 + (qJD(4) * t18 + t102) * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t84 - t187, 0, -t155, 0, 0, t94 * t141, -t164 * t175, 0, 0, t10, -t224, t14, t104 (t63 + t165) * t174 + t178, t148, -t63 * t68 - pkin(4) * t33 + t23 * t221 - t196 + (t159 + t195) * qJD(5) + (t113 * t91 - t19 * t94) * qJD(1), -t65 * t68 + pkin(4) * t32 - t24 * t221 + t197 + (-t160 + t194) * qJD(5) + (t113 * t93 + t20 * t94) * qJD(1), t23 * t65 + t24 * t63 - t29 + (-t19 * t175 - t110 + (-t19 + t206) * qJD(5)) * t93 + (t122 - t112) * t91, -pkin(4) * t46 + pkin(8) * t100 - t19 * t23 - t20 * t24 - t56 * t68, t10, t14, t224, t148, t109 - t50, t104, -t22 * t221 + t33 * t71 - t204 - t179 * t63 + (t18 * t91 + t159) * qJD(5) + (-t114 * t91 + t15 * t94) * qJD(1), t21 * t63 - t22 * t65 - t29 + (t15 * t175 + t1 + (t15 + t206) * qJD(5)) * t93 + (t122 + t217) * t91, t21 * t221 + t32 * t71 - t205 + t179 * t65 + (-t18 * t93 + t160) * qJD(5) + (t114 * t93 - t16 * t94) * qJD(1), t102 * pkin(8) - t15 * t22 - t16 * t21 - t179 * t18 + t5 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, -t146, t17, -t192, -t216, t80, -t56 * t65 + t112, -t19 * t221 + t56 * t63 + t110, 0, 0, t192, t17, t146, t80, t216, -t192, -t26 * t63 + t112 + 0.2e1 * t139 - t202, pkin(5) * t32 - qJ(6) * t33 + (t16 - t20) * t65 + (t15 - t163) * t63, 0.2e1 * t137 - t18 * t63 + t26 * t65 - (0.2e1 * qJD(6) - t19) * t221 - t110, -pkin(5) * t2 + qJ(6) * t1 - t15 * t20 + t16 * t163 - t18 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80 + t192, t17, -t221 ^ 2 - t207, t202 + t217;];
tauc_reg  = t3;
