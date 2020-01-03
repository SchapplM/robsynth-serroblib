% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR15_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:08
% EndTime: 2019-12-31 20:43:15
% DurationCPUTime: 2.01s
% Computational Cost: add. (1722->286), mult. (4076->418), div. (0->0), fcn. (2596->6), ass. (0->172)
t213 = pkin(3) + pkin(6);
t125 = sin(qJ(2));
t123 = sin(qJ(5));
t124 = sin(qJ(4));
t126 = cos(qJ(5));
t127 = cos(qJ(4));
t79 = t123 * t127 + t124 * t126;
t137 = t79 * t125;
t215 = qJD(4) + qJD(5);
t211 = -qJD(1) * t137 - t215 * t79;
t185 = qJD(2) * t124;
t128 = cos(qJ(2));
t186 = qJD(1) * t128;
t76 = t127 * t186 + t185;
t165 = t124 * t186;
t183 = qJD(2) * t127;
t78 = -t165 + t183;
t148 = t123 * t76 - t126 * t78;
t25 = t123 * t78 + t126 * t76;
t223 = t148 * t25;
t222 = t148 ^ 2 - t25 ^ 2;
t187 = qJD(1) * t125;
t109 = qJD(4) + t187;
t102 = qJD(5) + t109;
t177 = qJD(5) * t126;
t178 = qJD(5) * t123;
t175 = qJD(1) * qJD(2);
t162 = t125 * t175;
t39 = -qJD(4) * t76 + t124 * t162;
t180 = qJD(4) * t127;
t40 = qJD(2) * t180 - qJD(4) * t165 - t127 * t162;
t8 = -t123 * t40 + t126 * t39 - t76 * t177 - t78 * t178;
t221 = t102 * t25 + t8;
t129 = -pkin(2) - pkin(7);
t161 = -t125 * qJ(3) - pkin(1);
t72 = t129 * t128 + t161;
t46 = t72 * qJD(1);
t112 = pkin(6) * t187;
t216 = qJD(3) + t112;
t189 = pkin(3) * t187 + t216;
t51 = t129 * qJD(2) + t189;
t18 = t124 * t51 + t127 * t46;
t13 = -pkin(8) * t76 + t18;
t11 = t13 * t178;
t120 = qJD(2) * qJ(3);
t113 = pkin(6) * t186;
t86 = pkin(3) * t186 + t113;
t65 = t120 + t86;
t32 = pkin(4) * t76 + t65;
t220 = t25 * t32 + t11;
t111 = t128 * t175;
t108 = pkin(2) * t162;
t150 = pkin(7) * t125 - qJ(3) * t128;
t176 = t125 * qJD(3);
t136 = t150 * qJD(2) - t176;
t34 = qJD(1) * t136 + t108;
t107 = pkin(6) * t111;
t71 = pkin(3) * t111 + t107;
t158 = -t124 * t34 + t127 * t71;
t134 = -t18 * qJD(4) + t158;
t4 = pkin(4) * t111 - pkin(8) * t39 + t134;
t174 = -t124 * t71 - t127 * t34 - t51 * t180;
t181 = qJD(4) * t124;
t140 = -t46 * t181 - t174;
t5 = -pkin(8) * t40 + t140;
t170 = -t123 * t5 + t126 * t4;
t17 = -t124 * t46 + t127 * t51;
t12 = -pkin(8) * t78 + t17;
t10 = pkin(4) * t109 + t12;
t205 = t126 * t13;
t2 = t10 * t123 + t205;
t219 = -t2 * qJD(5) + t32 * t148 + t170;
t133 = t148 * qJD(5) - t123 * t39 - t126 * t40;
t218 = -t102 * t148 + t133;
t217 = -0.2e1 * t175;
t214 = t215 * t128;
t212 = pkin(8) - t129;
t169 = t127 * t187;
t194 = t126 * t127;
t196 = t123 * t124;
t210 = -t123 * t181 - t124 * t178 + t126 * t169 - t187 * t196 + t215 * t194;
t116 = pkin(2) * t187;
t56 = t150 * qJD(1) + t116;
t209 = t124 * t86 + t127 * t56;
t95 = t213 * t125;
t81 = t124 * t95;
t208 = t127 * t72 + t81;
t207 = qJD(2) * pkin(2);
t206 = t109 * t78;
t204 = t128 * t78;
t203 = t39 * t127;
t119 = qJD(2) * qJD(3);
t184 = qJD(2) * t125;
t85 = t213 * t184;
t53 = -qJD(1) * t85 + t119;
t202 = t53 * t124;
t201 = t53 * t127;
t200 = t76 * t109;
t171 = -pkin(4) * t127 - pkin(3);
t199 = pkin(4) * t180 - t171 * t187 + t216;
t92 = -pkin(2) * t128 + t161;
t66 = qJD(1) * t92;
t198 = t109 * t125;
t197 = t109 * t129;
t195 = t124 * t125;
t193 = t127 * t128;
t131 = qJD(1) ^ 2;
t192 = t128 * t131;
t130 = qJD(2) ^ 2;
t191 = t130 * t125;
t190 = t130 * t128;
t96 = t213 * t128;
t121 = t125 ^ 2;
t122 = t128 ^ 2;
t188 = t121 - t122;
t182 = qJD(2) * t128;
t179 = qJD(4) * t128;
t173 = t127 * t198;
t172 = t125 * t192;
t168 = t124 * t179;
t167 = t109 * t180;
t166 = t127 * t179;
t164 = pkin(8) * t128 - t72;
t91 = t212 * t127;
t163 = t210 * t102;
t160 = qJD(5) * t10 + t5;
t115 = pkin(2) * t184;
t42 = t115 + t136;
t87 = t213 * t182;
t157 = -t124 * t42 + t127 * t87;
t156 = -t124 * t56 + t127 * t86;
t155 = pkin(1) * t217;
t154 = qJD(3) - t207;
t146 = -t194 + t196;
t153 = t211 * t102 - t146 * t111;
t143 = pkin(4) * t128 - pkin(8) * t195;
t90 = t212 * t124;
t152 = t143 * qJD(1) - qJD(5) * t90 - t212 * t181 + t156;
t151 = pkin(8) * t169 + t215 * t91 + t209;
t82 = t127 * t95;
t21 = t125 * pkin(4) + t164 * t124 + t82;
t23 = -pkin(8) * t193 + t208;
t149 = t123 * t21 + t126 * t23;
t147 = -0.2e1 * qJD(2) * t66;
t145 = -qJD(1) * t122 + t198;
t144 = t109 * t124;
t138 = -qJ(3) * t182 - t176;
t44 = qJD(1) * t138 + t108;
t61 = t115 + t138;
t142 = pkin(6) * t130 + qJD(1) * t61 + t44;
t141 = t125 * t65 + t129 * t182;
t139 = t124 * t87 + t127 * t42 + t95 * t180 - t72 * t181;
t88 = pkin(6) * t162 - t119;
t89 = t112 + t154;
t94 = -t113 - t120;
t132 = -t88 * t128 + (t128 * t89 + (t94 + t113) * t125) * qJD(2);
t110 = pkin(4) * t124 + qJ(3);
t99 = t127 * t111;
t98 = t125 * t111;
t83 = -qJ(3) * t186 + t116;
t64 = pkin(4) * t193 + t96;
t58 = t79 * t128;
t57 = t146 * t128;
t50 = t66 * t187;
t33 = -pkin(4) * t168 + (-pkin(6) + t171) * t184;
t20 = t40 * pkin(4) + t53;
t15 = -t146 * t184 + t79 * t214;
t14 = qJD(2) * t137 + t146 * t214;
t7 = (t125 * t183 + t168) * pkin(8) + t139;
t6 = t143 * qJD(2) + (t164 * t127 - t81) * qJD(4) + t157;
t1 = t10 * t126 - t123 * t13;
t3 = [0, 0, 0, 0.2e1 * t98, t188 * t217, t190, -t191, 0, -pkin(6) * t190 + t125 * t155, pkin(6) * t191 + t128 * t155, t132, t125 * t147 + t142 * t128, -t142 * t125 + t128 * t147, pkin(6) * t132 + t44 * t92 + t61 * t66, -t78 * t166 + (-t128 * t39 + t78 * t184) * t124, (-t124 * t76 + t127 * t78) * t184 + (t124 * t40 - t203 + (t124 * t78 + t127 * t76) * qJD(4)) * t128, -t109 * t166 + t39 * t125 + (t124 * t145 + t204) * qJD(2), t109 * t168 - t40 * t125 + (t127 * t145 - t128 * t76) * qJD(2), t109 * t182 + t98, t157 * t109 - t85 * t76 + t96 * t40 + (-t65 * t183 + t158) * t125 + (-t208 * t109 - t18 * t125) * qJD(4) + (-t65 * t181 + t201 + ((-t124 * t72 + t82) * qJD(1) + t17) * qJD(2)) * t128, -t139 * t109 - t85 * t78 + t96 * t39 + ((qJD(2) * t65 + qJD(4) * t46) * t124 + t174) * t125 + (-t65 * t180 - t202 + (-t208 * qJD(1) - t18) * qJD(2)) * t128, -t14 * t148 - t58 * t8, -t133 * t58 - t14 * t25 - t148 * t15 + t57 * t8, t14 * t102 + t8 * t125 + (-qJD(1) * t58 - t148) * t182, t15 * t102 + t133 * t125 + (qJD(1) * t57 - t25) * t182, t102 * t182 + t98, (-t123 * t7 + t126 * t6) * t102 + t170 * t125 + t33 * t25 - t64 * t133 - t20 * t57 - t32 * t15 + (-t102 * t149 - t125 * t2) * qJD(5) + ((-t123 * t23 + t126 * t21) * qJD(1) + t1) * t182, t11 * t125 + t32 * t14 - t20 * t58 - t33 * t148 + t64 * t8 + (-(-qJD(5) * t23 + t6) * t102 - t4 * t125) * t123 + (-(qJD(5) * t21 + t7) * t102 - t160 * t125) * t126 + (-qJD(1) * t149 - t2) * t182; 0, 0, 0, -t172, t188 * t131, 0, 0, 0, t131 * pkin(1) * t125, pkin(1) * t192, ((-t94 - t120) * t125 + (t154 - t89) * t128) * qJD(1), -t83 * t186 + t50, 0.2e1 * t119 + (t125 * t83 + t128 * t66) * qJD(1), -qJ(3) * t88 - qJD(3) * t94 - t66 * t83 + (-t125 * t94 + (-t89 - t207) * t128) * qJD(1) * pkin(6), -t144 * t78 + t203, (-t40 - t206) * t127 + (-t39 + t200) * t124, -t109 * t181 + t99 + (-t109 * t195 - t204) * qJD(1), -t167 + (-t173 + (t76 - t185) * t128) * qJD(1), -t109 * t186, qJ(3) * t40 + t202 - t156 * t109 + t189 * t76 + (-t124 * t197 + t127 * t65) * qJD(4) + (t127 * t141 - t17 * t128) * qJD(1), qJ(3) * t39 + t201 + t209 * t109 + t189 * t78 + (-t124 * t65 - t127 * t197) * qJD(4) + (-t124 * t141 + t18 * t128) * qJD(1), -t146 * t8 - t148 * t211, -t133 * t146 + t148 * t210 - t211 * t25 - t79 * t8, t148 * t186 + t153, -t163 + (-qJD(2) * t79 + t25) * t186, -t102 * t186, -t110 * t133 + t20 * t79 + t210 * t32 + t199 * t25 + (t123 * t151 - t126 * t152) * t102 + ((t123 * t90 - t126 * t91) * qJD(2) - t1) * t186, t110 * t8 - t20 * t146 + t211 * t32 - t199 * t148 + (t123 * t152 + t126 * t151) * t102 + (-(-t123 * t91 - t126 * t90) * qJD(2) + t2) * t186; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t121 * t131 - t130, qJD(2) * t94 + t107 + t50, 0, 0, 0, 0, 0, -qJD(2) * t76 - t109 * t144 + t99, -t167 - qJD(2) * t78 + (-t124 * t182 - t173) * qJD(1), 0, 0, 0, 0, 0, -qJD(2) * t25 + t153, -t163 + (-t186 * t79 + t148) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t39 + t200, -t40 + t206, t111, t109 * t18 - t65 * t78 + t134, t109 * t17 + t65 * t76 - t140, -t223, t222, t221, t218, t111, -(-t12 * t123 - t205) * t102 + (-t102 * t178 + t111 * t126 - t25 * t78) * pkin(4) + t219, (-t102 * t13 - t4) * t123 + (t102 * t12 - t160) * t126 + (-t102 * t177 - t111 * t123 + t148 * t78) * pkin(4) + t220; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t222, t221, t218, t111, t102 * t2 + t219, t1 * t102 - t123 * t4 - t126 * t160 + t220;];
tauc_reg = t3;
