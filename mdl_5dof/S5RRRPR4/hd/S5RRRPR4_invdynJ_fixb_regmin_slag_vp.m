% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:42
% EndTime: 2019-12-31 21:11:48
% DurationCPUTime: 2.23s
% Computational Cost: add. (1666->287), mult. (2459->349), div. (0->0), fcn. (1494->10), ass. (0->185)
t123 = sin(qJ(5));
t128 = cos(qJ(3));
t124 = sin(qJ(3));
t127 = cos(qJ(5));
t202 = t124 * t127;
t229 = -t123 * t128 + t202;
t191 = qJD(3) * t128;
t192 = qJD(3) * t124;
t22 = t229 * qJD(5) + t123 * t191 - t127 * t192;
t122 = qJ(1) + qJ(2);
t108 = cos(t122);
t221 = g(2) * t108;
t109 = t124 * qJ(4);
t111 = t128 * pkin(3);
t199 = t111 + t109;
t187 = pkin(2) + t199;
t115 = qJD(3) - qJD(5);
t116 = qJD(1) + qJD(2);
t63 = t123 * t124 + t127 * t128;
t44 = t63 * t116;
t114 = qJDD(1) + qJDD(2);
t143 = t63 * qJD(5);
t177 = t116 * t191;
t201 = t128 * t114;
t235 = t116 * t192 - t201;
t9 = t114 * t202 - t116 * t143 + t235 * t123 + t127 * t177;
t237 = -t115 * t44 + t9;
t107 = sin(t122);
t216 = g(1) * t108 + g(2) * t107;
t117 = qJDD(3) * qJ(4);
t118 = qJD(3) * qJD(4);
t125 = sin(qJ(2));
t190 = qJDD(1) * t125;
t129 = cos(qJ(2));
t194 = qJD(2) * t129;
t48 = pkin(7) * t114 + (qJD(1) * t194 + t190) * pkin(1);
t37 = t128 * t48;
t215 = pkin(1) * qJD(1);
t183 = t125 * t215;
t68 = pkin(7) * t116 + t183;
t17 = -t68 * t192 + t117 + t118 + t37;
t36 = t124 * t48;
t54 = t68 * t191;
t180 = qJDD(4) + t36 + t54;
t209 = qJDD(3) * pkin(3);
t18 = t180 - t209;
t234 = t18 * t124 + t17 * t128;
t10 = t63 * t114 + t22 * t116;
t46 = t229 * t116;
t233 = t115 * t46 + t10;
t120 = t124 ^ 2;
t121 = t128 ^ 2;
t197 = t120 + t121;
t231 = t116 * t197;
t57 = t124 * t68;
t230 = qJD(4) + t57;
t228 = qJD(5) + t115;
t227 = pkin(3) + pkin(4);
t226 = pkin(7) - pkin(8);
t225 = pkin(1) * t129;
t224 = pkin(2) * t114;
t223 = pkin(2) * t116;
t132 = qJD(3) ^ 2;
t222 = pkin(7) * t132;
t98 = g(1) * t107;
t220 = t46 * t44;
t101 = pkin(1) * t125 + pkin(7);
t219 = -pkin(8) + t101;
t206 = t108 * t124;
t207 = t107 * t124;
t217 = g(1) * t207 - g(2) * t206;
t58 = t128 * t68;
t212 = -qJD(2) * t183 + qJDD(1) * t225;
t211 = pkin(7) * qJDD(3);
t210 = qJ(4) * t128;
t208 = t101 * t132;
t205 = t115 * t129;
t204 = t116 * t124;
t92 = t124 * t114;
t200 = -pkin(8) * t204 + t230;
t198 = t120 - t121;
t196 = qJD(1) * t129;
t195 = qJD(2) * t125;
t119 = qJD(3) * qJ(4);
t193 = qJD(3) * t116;
t106 = t124 * qJD(4);
t189 = qJDD(3) * t101;
t168 = t116 * t183;
t186 = pkin(1) * t196;
t88 = t128 * t98;
t188 = t128 * t168 + t186 * t192 + t88;
t185 = pkin(1) * t195;
t184 = pkin(1) * t194;
t182 = t44 ^ 2 - t46 ^ 2;
t82 = t226 * t128;
t112 = t116 ^ 2;
t181 = t124 * t112 * t128;
t52 = pkin(3) * t192 - qJ(4) * t191 - t106;
t102 = -pkin(2) - t225;
t179 = t116 * t195;
t47 = -t212 - t224;
t174 = -t47 - t221;
t62 = t219 * t128;
t173 = pkin(2) + t109;
t69 = -t186 - t223;
t172 = t47 * t124 + t69 * t191 - t217;
t171 = t197 * t114;
t170 = -qJD(3) * pkin(3) + qJD(4);
t169 = t115 ^ 2;
t33 = -pkin(8) * t116 * t128 + t58;
t167 = t222 - t224;
t126 = sin(qJ(1));
t130 = cos(qJ(1));
t165 = g(1) * t126 - g(2) * t130;
t164 = pkin(3) * t124 - t210;
t21 = -t227 * qJD(3) + t200;
t27 = t119 + t33;
t163 = t123 * t27 - t127 * t21;
t162 = -t123 * t21 - t127 * t27;
t61 = t219 * t124;
t161 = -t123 * t62 + t127 * t61;
t160 = t123 * t61 + t127 * t62;
t81 = t226 * t124;
t159 = -t123 * t82 + t127 * t81;
t158 = t123 * t81 + t127 * t82;
t43 = t170 + t57;
t51 = t58 + t119;
t157 = t124 * t43 + t128 * t51;
t156 = g(1) * t206 + g(2) * t207 - g(3) * t128 - t36;
t59 = t102 - t199;
t155 = t212 + t98 - t221;
t154 = qJ(4) * t127 - t123 * t227;
t153 = -qJ(4) * t123 - t127 * t227;
t151 = -t173 - t111;
t12 = (t164 * qJD(3) - t106) * t116 + t151 * t114 - t212;
t150 = t114 * t187 - t12 - t222;
t142 = t227 * t128 + t173;
t19 = t142 * t116 + t186;
t40 = t63 * t107;
t42 = t63 * t108;
t146 = -t227 * t124 + t210;
t5 = (t146 * qJD(3) + t106) * t116 + t142 * t114 + t212;
t149 = g(1) * t40 - g(2) * t42 + t19 * t22 + t5 * t63;
t23 = t63 * qJD(3) - t143;
t39 = t229 * t107;
t41 = t229 * t108;
t148 = g(1) * t39 - g(2) * t41 + t19 * t23 + t229 * t5;
t147 = t116 * t59 - t184;
t144 = -qJDD(4) + t156;
t34 = -pkin(4) * t192 - t52;
t140 = t43 * t191 - t51 * t192 - t216 + t234;
t35 = t52 + t185;
t139 = -t114 * t59 - t116 * t35 - t12 - t208;
t138 = pkin(1) * t179 + t102 * t114 + t208;
t137 = -t189 + (t102 * t116 - t184) * qJD(3);
t136 = (-t124 * t51 + t128 * t43) * qJD(3) + t234;
t11 = t235 * pkin(8) + t17;
t8 = -t227 * qJDD(3) + (-t177 - t92) * pkin(8) + t180;
t135 = -g(1) * t41 - g(2) * t39 + g(3) * t63 - t123 * t11 + t127 * t8 - t19 * t46;
t134 = g(1) * t42 + g(2) * t40 + g(3) * t229 - t127 * t11 - t123 * t8 + t19 * t44;
t133 = -t216 * pkin(7) - t151 * t98 - t187 * t221;
t113 = qJDD(3) - qJDD(5);
t110 = t128 * pkin(4);
t104 = pkin(8) * t192;
t77 = qJDD(3) * t128 - t124 * t132;
t76 = qJDD(3) * t124 + t128 * t132;
t66 = qJD(3) * t82;
t65 = -pkin(7) * t192 + t104;
t60 = t110 + t187;
t55 = t69 * t192;
t53 = t164 * t116;
t50 = t114 * t120 + 0.2e1 * t124 * t177;
t49 = t110 - t59;
t31 = t146 * t116;
t29 = qJD(3) * t62 + t124 * t184;
t28 = -t101 * t192 + t128 * t184 + t104;
t26 = t151 * t116 - t186;
t25 = t34 - t185;
t24 = 0.2e1 * t124 * t201 - 0.2e1 * t198 * t193;
t20 = t26 * t192;
t14 = -t113 * t229 - t115 * t23;
t13 = t113 * t63 + t115 * t22;
t2 = t229 * t9 + t23 * t46;
t1 = -t10 * t229 - t22 * t46 - t23 * t44 - t63 * t9;
t3 = [qJDD(1), t165, g(1) * t130 + g(2) * t126, t114, (t114 * t129 - t179) * pkin(1) + t155, ((-qJDD(1) - t114) * t125 + (-qJD(1) - t116) * t194) * pkin(1) + t216, t50, t24, t76, t77, 0, t55 + t88 + t137 * t124 + (-t138 + t174) * t128, t138 * t124 + t137 * t128 + t172, t20 + t88 + (t147 * qJD(3) - t189) * t124 + (t139 - t221) * t128, t101 * t171 + t184 * t231 + t140, (t189 + (-t147 - t26) * qJD(3)) * t128 + t139 * t124 + t217, t12 * t59 + t26 * t35 + (t157 * t194 + t165) * pkin(1) + t136 * t101 + t133, t2, t1, t14, t13, 0, t25 * t44 + t49 * t10 - (-qJD(5) * t160 - t123 * t28 + t127 * t29) * t115 - t161 * t113 + t149, t25 * t46 + t49 * t9 + (qJD(5) * t161 + t123 * t29 + t127 * t28) * t115 + t160 * t113 + t148; 0, 0, 0, t114, t155 + t168, (-t190 + (-qJD(2) + t116) * t196) * pkin(1) + t216, t50, t24, t76, t77, 0, t55 + (-pkin(2) * t193 - t211) * t124 + (-t167 + t174) * t128 + t188, (-t211 + (t186 - t223) * qJD(3)) * t128 + (t167 - t168) * t124 + t172, t20 + (-t187 * t193 - t211) * t124 + (-t116 * t52 + t150 - t221) * t128 + t188, pkin(7) * t171 - t186 * t231 + t140, (t211 + (t116 * t187 - t186 - t26) * qJD(3)) * t128 + ((-t52 + t183) * t116 + t150) * t124 + t217, -t12 * t187 + t26 * t52 + (-t125 * t26 - t129 * t157) * t215 + t136 * pkin(7) + t133, t2, t1, t14, t13, 0, t34 * t44 + t60 * t10 - (-qJD(5) * t158 - t123 * t65 + t127 * t66) * t115 - t159 * t113 + (t125 * t44 + t205 * t229) * t215 + t149, t34 * t46 + t60 * t9 + (qJD(5) * t159 + t123 * t66 + t127 * t65) * t115 + t158 * t113 + (t125 * t46 - t205 * t63) * t215 + t148; 0, 0, 0, 0, 0, 0, -t181, t198 * t112, t92, t201, qJDD(3), -t69 * t204 + t156, g(3) * t124 - t37 + (-t116 * t69 + t216) * t128, 0.2e1 * t209 + (-t124 * t26 + t128 * t53) * t116 + t144, -t164 * t114 + ((t51 - t119) * t124 + (t170 - t43) * t128) * t116, 0.2e1 * t117 + 0.2e1 * t118 + t37 + (t116 * t53 - g(3)) * t124 + (t116 * t26 - t216) * t128, -t18 * pkin(3) - g(3) * t199 + t17 * qJ(4) + t216 * t164 + t230 * t51 - t26 * t53 - t43 * t58, -t220, t182, -t237, t233, t113, -t153 * t113 - t31 * t44 + (t123 * t200 + t127 * t33) * t115 + (t115 * t154 - t162) * qJD(5) - t135, t154 * t113 - t31 * t46 + (-t123 * t33 + t127 * t200) * t115 + (t115 * t153 - t163) * qJD(5) - t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(3) - t181, t92, -t112 * t120 - t132, -qJD(3) * t51 + t204 * t26 - t144 - t209 + t54, 0, 0, 0, 0, 0, -t127 * t113 - t123 * t169 - t204 * t44, t123 * t113 - t127 * t169 - t204 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t220, -t182, t237, -t233, -t113, t228 * t162 + t135, t228 * t163 + t134;];
tau_reg = t3;
