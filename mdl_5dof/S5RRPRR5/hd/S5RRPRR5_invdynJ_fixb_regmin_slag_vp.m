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
% Datum: 2019-12-05 18:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:34:07
% EndTime: 2019-12-05 18:34:11
% DurationCPUTime: 1.59s
% Computational Cost: add. (2455->248), mult. (3712->312), div. (0->0), fcn. (2849->16), ass. (0->165)
t144 = sin(qJ(5));
t148 = cos(qJ(5));
t140 = qJD(1) + qJD(2);
t143 = cos(pkin(9));
t149 = cos(qJ(4));
t210 = t149 * t143;
t193 = t140 * t210;
t142 = sin(pkin(9));
t145 = sin(qJ(4));
t213 = t145 * t142;
t194 = t140 * t213;
t71 = -t193 + t194;
t92 = t149 * t142 + t145 * t143;
t73 = t92 * t140;
t167 = t144 * t71 - t148 * t73;
t135 = qJDD(1) + qJDD(2);
t195 = qJD(4) * t193 + t92 * t135;
t37 = -qJD(4) * t194 + t195;
t91 = -t210 + t213;
t171 = t91 * t135;
t82 = t92 * qJD(4);
t38 = t140 * t82 + t171;
t156 = qJD(5) * t167 - t144 * t37 - t148 * t38;
t139 = qJD(4) + qJD(5);
t218 = t167 * t139;
t240 = t156 - t218;
t150 = cos(qJ(2));
t205 = qJD(1) * t150;
t196 = pkin(1) * t205;
t174 = qJD(3) - t196;
t30 = t144 * t73 + t148 * t71;
t217 = t30 * t139;
t200 = qJD(5) * t148;
t201 = qJD(5) * t144;
t7 = -t144 * t38 + t148 * t37 - t71 * t200 - t201 * t73;
t239 = t7 + t217;
t238 = t167 * t30;
t206 = t142 ^ 2 + t143 ^ 2;
t146 = sin(qJ(2));
t199 = qJDD(1) * t146;
t203 = qJD(2) * t150;
t59 = t135 * qJ(3) + t140 * qJD(3) + (qJD(1) * t203 + t199) * pkin(1);
t185 = t206 * t59;
t141 = qJ(1) + qJ(2);
t130 = sin(t141);
t131 = cos(t141);
t173 = g(2) * t131 + g(3) * t130;
t123 = g(2) * t130;
t231 = g(3) * t131 - t123;
t237 = t167 ^ 2 - t30 ^ 2;
t138 = pkin(9) + qJ(4);
t129 = qJ(5) + t138;
t117 = sin(t129);
t224 = t146 * pkin(1);
t197 = qJD(1) * t224;
t98 = t140 * qJ(3) + t197;
t190 = pkin(7) * t140 + t98;
t60 = t190 * t142;
t61 = t190 * t143;
t165 = t145 * t60 - t149 * t61;
t21 = -t71 * pkin(8) - t165;
t118 = cos(t129);
t214 = t118 * t131;
t215 = t118 * t130;
t120 = -t143 * pkin(3) - pkin(2);
t69 = t120 * t140 + t174;
t41 = t71 * pkin(4) + t69;
t191 = pkin(7) * t135 + t59;
t45 = t191 * t142;
t46 = t191 * t143;
t187 = -t145 * t46 - t149 * t45;
t5 = qJDD(4) * pkin(4) - t37 * pkin(8) + qJD(4) * t165 + t187;
t236 = t41 * t30 + g(1) * t117 + t21 * t201 + g(3) * t214 + (-t21 * t139 - t5) * t144 - g(2) * t215;
t233 = -t145 * t61 - t149 * t60;
t20 = -t73 * pkin(8) + t233;
t19 = qJD(4) * pkin(4) + t20;
t219 = t148 * t21;
t170 = -t144 * t19 - t219;
t166 = -t145 * t45 + t149 * t46;
t6 = -t38 * pkin(8) + t233 * qJD(4) + t166;
t235 = -g(1) * t118 + qJD(5) * t170 + t231 * t117 - t144 * t6 + t148 * t5 + t41 * t167;
t223 = t150 * pkin(1);
t207 = -qJD(2) * t197 + qJDD(1) * t223;
t189 = qJDD(3) - t207;
t225 = t135 * pkin(2);
t70 = t189 - t225;
t234 = t173 - t70;
t103 = (-pkin(7) - qJ(3)) * t142;
t132 = t143 * pkin(7);
t104 = t143 * qJ(3) + t132;
t221 = t145 * t103 + t149 * t104;
t232 = -t221 * qJD(4) - t174 * t92;
t202 = qJD(4) * t149;
t230 = (qJD(3) * t142 + qJD(4) * t104) * t145 - t91 * t196 - qJD(3) * t210 - t103 * t202;
t81 = t91 * qJD(4);
t229 = t81 * pkin(8);
t228 = t82 * pkin(4);
t227 = t91 * pkin(4);
t226 = t92 * pkin(8);
t119 = qJ(3) + t224;
t83 = (-pkin(7) - t119) * t142;
t84 = t143 * t119 + t132;
t222 = t145 * t83 + t149 * t84;
t209 = t173 * t143;
t204 = qJD(2) * t146;
t198 = pkin(1) * t204;
t192 = t140 * t204;
t186 = -t145 * t84 + t149 * t83;
t52 = -t144 * t91 + t148 * t92;
t23 = qJD(5) * t52 - t144 * t81 + t148 * t82;
t53 = t120 * t135 + t189;
t24 = t38 * pkin(4) + t53;
t51 = t144 * t92 + t148 * t91;
t184 = g(2) * t214 + g(3) * t215 + t41 * t23 + t24 * t51;
t183 = t149 * t103 - t145 * t104;
t111 = pkin(1) * t203 + qJD(3);
t182 = t111 * t206;
t181 = t206 * t135;
t128 = cos(t138);
t180 = t173 * t128 + t53 * t91 + t69 * t82;
t179 = -t231 + t185;
t178 = t140 * t197;
t177 = -t207 - t173;
t42 = t183 - t226;
t79 = t82 * pkin(8);
t176 = -qJD(5) * t42 + t230 + t79;
t86 = t91 * pkin(8);
t43 = -t86 + t221;
t175 = qJD(5) * t43 - t229 - t232;
t27 = t186 - t226;
t28 = -t86 + t222;
t169 = -t144 * t28 + t148 * t27;
t168 = t144 * t27 + t148 * t28;
t102 = t120 - t223;
t163 = -t197 + t228;
t162 = t178 + t225;
t125 = -pkin(2) - t223;
t161 = -pkin(1) * t192 - t125 * t135;
t22 = -qJD(5) * t51 - t144 * t82 - t148 * t81;
t160 = -t117 * t173 + t41 * t22 + t24 * t52;
t127 = sin(t138);
t159 = -t127 * t173 + t53 * t92 - t69 * t81;
t158 = t83 * t202 + t111 * t210 + (-qJD(4) * t84 - t111 * t142) * t145;
t157 = t174 * t206;
t154 = -t222 * qJD(4) - t92 * t111;
t151 = cos(qJ(1));
t147 = sin(qJ(1));
t134 = qJDD(4) + qJDD(5);
t115 = t131 * qJ(3);
t93 = -t140 * pkin(2) + t174;
t68 = t120 + t227;
t63 = t198 + t228;
t62 = t70 * t142;
t57 = t102 + t227;
t50 = -t82 * qJD(4) - t91 * qJDD(4);
t49 = -t81 * qJD(4) + t92 * qJDD(4);
t17 = t154 + t229;
t16 = -t79 + t158;
t15 = t37 * t92 - t73 * t81;
t10 = -t51 * t134 - t23 * t139;
t9 = t52 * t134 + t22 * t139;
t4 = -t37 * t91 - t92 * t38 + t81 * t71 - t73 * t82;
t2 = -t167 * t22 + t7 * t52;
t1 = t156 * t52 + t167 * t23 - t22 * t30 - t7 * t51;
t3 = [qJDD(1), g(2) * t151 + g(3) * t147, -g(2) * t147 + g(3) * t151, t135, (t135 * t150 - t192) * pkin(1) - t177, ((-qJDD(1) - t135) * t146 + (-qJD(1) - t140) * t203) * pkin(1) + t231, (t161 - t70) * t143 + t209, t62 + (-t161 - t173) * t142, t119 * t181 + t140 * t182 + t179, t70 * t125 + t93 * t198 - g(2) * (-t151 * pkin(1) - t131 * pkin(2) - t130 * qJ(3)) - g(3) * (-t147 * pkin(1) - t130 * pkin(2) + t115) + t119 * t185 + t98 * t182, t15, t4, t49, t50, 0, qJD(4) * t154 + qJDD(4) * t186 + t102 * t38 + t198 * t71 + t180, -t158 * qJD(4) - t222 * qJDD(4) + t102 * t37 + t73 * t198 + t159, t2, t1, t9, t10, 0, t63 * t30 - t57 * t156 + (-qJD(5) * t168 - t144 * t16 + t148 * t17) * t139 + t169 * t134 + t184, -t63 * t167 + t57 * t7 - (qJD(5) * t169 + t144 * t17 + t148 * t16) * t139 - t168 * t134 + t160; 0, 0, 0, t135, -t177 + t178, (-t199 + (-qJD(2) + t140) * t205) * pkin(1) + t231, (t162 - t70) * t143 + t209, t62 + (-t162 - t173) * t142, qJ(3) * t181 + t140 * t157 + t179, -t93 * t197 - g(3) * t115 + t234 * pkin(2) + (t185 + t123) * qJ(3) + t157 * t98, t15, t4, t49, t50, 0, t232 * qJD(4) + t183 * qJDD(4) + t120 * t38 - t71 * t197 + t180, t230 * qJD(4) - t221 * qJDD(4) + t120 * t37 - t73 * t197 + t159, t2, t1, t9, t10, 0, -t68 * t156 + (-t144 * t43 + t148 * t42) * t134 + t163 * t30 + (t144 * t176 - t148 * t175) * t139 + t184, t68 * t7 - (t144 * t42 + t148 * t43) * t134 - t163 * t167 + (t144 * t175 + t148 * t176) * t139 + t160; 0, 0, 0, 0, 0, 0, -t143 * t135, t142 * t135, -t206 * t140 ^ 2, -t140 * t206 * t98 - t234, 0, 0, 0, 0, 0, 0.2e1 * t73 * qJD(4) + t171, (-t71 - t194) * qJD(4) + t195, 0, 0, 0, 0, 0, -t156 - t218, t7 - t217; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73 * t71, -t71 ^ 2 + t73 ^ 2, (t71 - t194) * qJD(4) + t195, -t171, qJDD(4), -g(1) * t128 + t127 * t231 - t69 * t73 + t187, g(1) * t127 + t128 * t231 + t69 * t71 - t166, -t238, t237, t239, t240, t134, -(-t144 * t20 - t219) * t139 + (t148 * t134 - t139 * t201 - t73 * t30) * pkin(4) + t235, (-qJD(5) * t19 + t20 * t139 - t6) * t148 + (-t144 * t134 - t139 * t200 + t167 * t73) * pkin(4) + t236; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t237, t239, t240, t134, -t139 * t170 + t235, (-t6 + (-qJD(5) + t139) * t19) * t148 + t236;];
tau_reg = t3;
