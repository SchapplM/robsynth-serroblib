% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:04:00
% EndTime: 2020-01-03 12:04:04
% DurationCPUTime: 1.50s
% Computational Cost: add. (2455->249), mult. (3712->312), div. (0->0), fcn. (2849->16), ass. (0->165)
t141 = qJ(1) + qJ(2);
t130 = sin(t141);
t131 = cos(t141);
t174 = g(2) * t131 + g(3) * t130;
t146 = sin(qJ(2));
t225 = t146 * pkin(1);
t198 = qJD(1) * t225;
t150 = cos(qJ(2));
t224 = t150 * pkin(1);
t208 = -qJD(2) * t198 + qJDD(1) * t224;
t189 = qJDD(3) - t208;
t135 = qJDD(1) + qJDD(2);
t226 = t135 * pkin(2);
t70 = t189 - t226;
t243 = t70 + t174;
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t140 = qJD(1) + qJD(2);
t143 = cos(pkin(9));
t149 = cos(qJ(4));
t211 = t149 * t143;
t193 = t140 * t211;
t142 = sin(pkin(9));
t145 = sin(qJ(4));
t214 = t145 * t142;
t194 = t140 * t214;
t71 = -t193 + t194;
t92 = t149 * t142 + t145 * t143;
t73 = t92 * t140;
t168 = t144 * t71 - t148 * t73;
t196 = qJD(4) * t193 + t92 * t135;
t37 = -qJD(4) * t194 + t196;
t91 = -t211 + t214;
t172 = t91 * t135;
t82 = t92 * qJD(4);
t38 = t140 * t82 + t172;
t156 = t168 * qJD(5) - t144 * t37 - t148 * t38;
t139 = qJD(4) + qJD(5);
t219 = t168 * t139;
t242 = t156 - t219;
t206 = qJD(1) * t150;
t197 = pkin(1) * t206;
t175 = qJD(3) - t197;
t30 = t144 * t73 + t148 * t71;
t218 = t30 * t139;
t201 = qJD(5) * t148;
t202 = qJD(5) * t144;
t7 = -t144 * t38 + t148 * t37 - t71 * t201 - t73 * t202;
t241 = t7 + t218;
t240 = t168 * t30;
t207 = t142 ^ 2 + t143 ^ 2;
t200 = qJDD(1) * t146;
t204 = qJD(2) * t150;
t59 = t135 * qJ(3) + t140 * qJD(3) + (qJD(1) * t204 + t200) * pkin(1);
t185 = t207 * t59;
t123 = g(3) * t131;
t234 = g(2) * t130 - t123;
t239 = t168 ^ 2 - t30 ^ 2;
t138 = pkin(9) + qJ(4);
t129 = qJ(5) + t138;
t117 = sin(t129);
t118 = cos(t129);
t98 = t140 * qJ(3) + t198;
t190 = pkin(7) * t140 + t98;
t60 = t190 * t142;
t61 = t190 * t143;
t166 = t145 * t60 - t149 * t61;
t21 = -t71 * pkin(8) - t166;
t120 = -t143 * pkin(3) - pkin(2);
t69 = t120 * t140 + t175;
t41 = t71 * pkin(4) + t69;
t191 = pkin(7) * t135 + t59;
t45 = t191 * t142;
t46 = t191 * t143;
t187 = -t145 * t46 - t149 * t45;
t5 = qJDD(4) * pkin(4) - t37 * pkin(8) + t166 * qJD(4) + t187;
t238 = t41 * t30 + g(1) * t117 + t21 * t202 + (-t21 * t139 - t5) * t144 + t234 * t118;
t236 = -t145 * t61 - t149 * t60;
t20 = -t73 * pkin(8) + t236;
t19 = qJD(4) * pkin(4) + t20;
t220 = t148 * t21;
t171 = -t144 * t19 - t220;
t215 = t117 * t131;
t216 = t117 * t130;
t167 = -t145 * t45 + t149 * t46;
t6 = -t38 * pkin(8) + t236 * qJD(4) + t167;
t237 = -g(1) * t118 + g(2) * t216 - g(3) * t215 + t171 * qJD(5) - t144 * t6 + t148 * t5 + t41 * t168;
t103 = (-pkin(7) - qJ(3)) * t142;
t132 = t143 * pkin(7);
t104 = t143 * qJ(3) + t132;
t222 = t145 * t103 + t149 * t104;
t235 = -t222 * qJD(4) - t175 * t92;
t203 = qJD(4) * t149;
t233 = (qJD(3) * t142 + qJD(4) * t104) * t145 - t91 * t197 - qJD(3) * t211 - t103 * t203;
t81 = t91 * qJD(4);
t232 = t81 * pkin(8);
t231 = t82 * pkin(4);
t230 = t91 * pkin(4);
t229 = t92 * pkin(8);
t119 = qJ(3) + t225;
t83 = (-pkin(7) - t119) * t142;
t84 = t143 * t119 + t132;
t223 = t145 * t83 + t149 * t84;
t210 = t131 * pkin(2) + t130 * qJ(3);
t205 = qJD(2) * t146;
t199 = pkin(1) * t205;
t195 = t243 * t142;
t192 = t140 * t205;
t186 = -t145 * t84 + t149 * t83;
t51 = t144 * t92 + t148 * t91;
t22 = -t51 * qJD(5) - t144 * t82 - t148 * t81;
t53 = t120 * t135 + t189;
t24 = t38 * pkin(4) + t53;
t52 = -t144 * t91 + t148 * t92;
t184 = g(2) * t215 + g(3) * t216 + t41 * t22 + t24 * t52;
t183 = t149 * t103 - t145 * t104;
t111 = pkin(1) * t204 + qJD(3);
t182 = t111 * t207;
t181 = t207 * t135;
t127 = sin(t138);
t180 = t174 * t127 + t53 * t92 - t69 * t81;
t179 = -t234 + t185;
t178 = t140 * t198;
t42 = t183 - t229;
t79 = t82 * pkin(8);
t177 = -qJD(5) * t42 + t233 + t79;
t86 = t91 * pkin(8);
t43 = -t86 + t222;
t176 = qJD(5) * t43 - t232 - t235;
t27 = t186 - t229;
t28 = -t86 + t223;
t170 = -t144 * t28 + t148 * t27;
t169 = t144 * t27 + t148 * t28;
t102 = t120 - t224;
t164 = -t198 + t231;
t163 = t174 - t208;
t125 = -pkin(2) - t224;
t162 = pkin(1) * t192 + t125 * t135;
t23 = t52 * qJD(5) - t144 * t81 + t148 * t82;
t161 = -t174 * t118 + t41 * t23 + t24 * t51;
t128 = cos(t138);
t160 = -t174 * t128 + t53 * t91 + t69 * t82;
t159 = t83 * t203 + t111 * t211 + (-qJD(4) * t84 - t111 * t142) * t145;
t158 = -t174 + t178;
t157 = t175 * t207;
t154 = -t223 * qJD(4) - t92 * t111;
t151 = cos(qJ(1));
t147 = sin(qJ(1));
t134 = qJDD(4) + qJDD(5);
t121 = t130 * pkin(2);
t93 = -t140 * pkin(2) + t175;
t68 = t120 + t230;
t63 = t199 + t231;
t57 = t102 + t230;
t50 = -t82 * qJD(4) - t91 * qJDD(4);
t49 = -t81 * qJD(4) + t92 * qJDD(4);
t17 = t154 + t232;
t16 = -t79 + t159;
t15 = t37 * t92 - t73 * t81;
t10 = -t51 * t134 - t23 * t139;
t9 = t52 * t134 + t22 * t139;
t4 = -t37 * t91 - t92 * t38 + t81 * t71 - t73 * t82;
t2 = -t168 * t22 + t7 * t52;
t1 = t156 * t52 + t168 * t23 - t22 * t30 - t7 * t51;
t3 = [qJDD(1), -g(2) * t151 - g(3) * t147, g(2) * t147 - g(3) * t151, t135, (t135 * t150 - t192) * pkin(1) - t163, ((-qJDD(1) - t135) * t146 + (-qJD(1) - t140) * t204) * pkin(1) + t234, (-t162 - t243) * t143, t162 * t142 + t195, t119 * t181 + t140 * t182 + t179, t70 * t125 + t93 * t199 - g(2) * (t151 * pkin(1) + t210) - g(3) * (t147 * pkin(1) - t131 * qJ(3) + t121) + t119 * t185 + t98 * t182, t15, t4, t49, t50, 0, t154 * qJD(4) + t186 * qJDD(4) + t102 * t38 + t71 * t199 + t160, -t159 * qJD(4) - t223 * qJDD(4) + t102 * t37 + t73 * t199 + t180, t2, t1, t9, t10, 0, t63 * t30 - t57 * t156 + (-t169 * qJD(5) - t144 * t16 + t148 * t17) * t139 + t170 * t134 + t161, -t63 * t168 + t57 * t7 - (qJD(5) * t170 + t144 * t17 + t148 * t16) * t139 - t169 * t134 + t184; 0, 0, 0, t135, t158 + t208, (-t200 + (-qJD(2) + t140) * t206) * pkin(1) + t234, (t158 - t70 + t226) * t143, (-t178 - t226) * t142 + t195, qJ(3) * t181 + t157 * t140 + t179, -t70 * pkin(2) - t93 * t198 - g(2) * t210 - g(3) * t121 + (t185 + t123) * qJ(3) + t157 * t98, t15, t4, t49, t50, 0, t235 * qJD(4) + t183 * qJDD(4) + t120 * t38 - t71 * t198 + t160, t233 * qJD(4) - t222 * qJDD(4) + t120 * t37 - t73 * t198 + t180, t2, t1, t9, t10, 0, -t68 * t156 + (-t144 * t43 + t148 * t42) * t134 + t164 * t30 + (t177 * t144 - t176 * t148) * t139 + t161, t68 * t7 - (t144 * t42 + t148 * t43) * t134 - t164 * t168 + (t144 * t176 + t148 * t177) * t139 + t184; 0, 0, 0, 0, 0, 0, -t143 * t135, t142 * t135, -t207 * t140 ^ 2, -t207 * t98 * t140 + qJDD(3) + t163 - t226, 0, 0, 0, 0, 0, 0.2e1 * t73 * qJD(4) + t172, (-t71 - t194) * qJD(4) + t196, 0, 0, 0, 0, 0, -t156 - t219, t7 - t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, (t71 - t194) * qJD(4) + t196, -t172, qJDD(4), -g(1) * t128 + t127 * t234 - t69 * t73 + t187, g(1) * t127 + t128 * t234 + t69 * t71 - t167, -t240, t239, t241, t242, t134, -(-t144 * t20 - t220) * t139 + (t148 * t134 - t139 * t202 - t73 * t30) * pkin(4) + t237, (-qJD(5) * t19 + t20 * t139 - t6) * t148 + (-t144 * t134 - t139 * t201 + t168 * t73) * pkin(4) + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t239, t241, t242, t134, -t139 * t171 + t237, (-t6 + (-qJD(5) + t139) * t19) * t148 + t238;];
tau_reg = t3;
