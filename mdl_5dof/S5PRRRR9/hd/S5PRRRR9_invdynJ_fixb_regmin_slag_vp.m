% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:14
% EndTime: 2019-12-05 17:21:24
% DurationCPUTime: 3.29s
% Computational Cost: add. (2334->375), mult. (5448->566), div. (0->0), fcn. (4346->14), ass. (0->194)
t146 = cos(qJ(4));
t142 = sin(qJ(4));
t223 = qJD(3) * t142;
t143 = sin(qJ(3));
t225 = qJD(2) * t143;
t100 = t146 * t225 + t223;
t141 = sin(qJ(5));
t145 = cos(qJ(5));
t214 = t146 * qJD(3);
t98 = t142 * t225 - t214;
t174 = t145 * t100 - t141 * t98;
t40 = t100 * t141 + t145 * t98;
t279 = t174 * t40;
t218 = qJD(4) * t143;
t278 = -qJD(2) * t218 + qJDD(3);
t277 = t174 ^ 2 - t40 ^ 2;
t147 = cos(qJ(3));
t224 = qJD(2) * t147;
t123 = -qJD(4) + t224;
t119 = -qJD(5) + t123;
t215 = qJD(5) * t145;
t216 = qJD(5) * t141;
t212 = qJD(2) * qJD(3);
t193 = t147 * t212;
t210 = t143 * qJDD(2);
t37 = qJD(4) * t214 + (t193 + t210) * t146 + t278 * t142;
t38 = ((qJD(4) + t224) * qJD(3) + t210) * t142 - t278 * t146;
t6 = -t100 * t216 - t141 * t38 + t145 * t37 - t98 * t215;
t276 = -t119 * t40 + t6;
t144 = sin(qJ(2));
t139 = sin(pkin(5));
t229 = qJD(1) * t139;
t108 = qJD(2) * pkin(7) + t144 * t229;
t140 = cos(pkin(5));
t228 = qJD(1) * t140;
t121 = t143 * t228;
t68 = t147 * t108 + t121;
t62 = qJD(3) * pkin(8) + t68;
t110 = -pkin(3) * t147 - pkin(8) * t143 - pkin(2);
t148 = cos(qJ(2));
t227 = qJD(1) * t148;
t197 = t139 * t227;
t70 = t110 * qJD(2) - t197;
t28 = t142 * t70 + t146 * t62;
t14 = -pkin(9) * t98 + t28;
t10 = t14 * t216;
t137 = qJ(4) + qJ(5);
t133 = sin(t137);
t134 = cos(t137);
t236 = t139 * t148;
t271 = -t143 * t108 + t147 * t228;
t61 = -qJD(3) * pkin(3) - t271;
t33 = pkin(4) * t98 + t61;
t247 = cos(pkin(10));
t188 = t139 * t247;
t187 = t247 * t144;
t138 = sin(pkin(10));
t239 = t138 * t148;
t79 = t140 * t187 + t239;
t51 = -t143 * t188 + t79 * t147;
t186 = t247 * t148;
t240 = t138 * t144;
t81 = -t140 * t240 + t186;
t53 = t138 * t139 * t143 + t147 * t81;
t78 = -t140 * t186 + t240;
t80 = t140 * t239 + t187;
t237 = t139 * t147;
t86 = t140 * t143 + t144 * t237;
t275 = t33 * t40 - g(1) * (-t133 * t80 - t134 * t53) - g(2) * (-t133 * t78 - t134 * t51) - g(3) * (t133 * t236 - t134 * t86) + t10;
t181 = pkin(3) * t143 - pkin(8) * t147;
t106 = t181 * qJD(3);
t222 = qJD(3) * t143;
t232 = t147 * t148;
t261 = pkin(7) * t142;
t273 = (-t142 * t232 + t144 * t146) * t229 - t146 * t106 - t222 * t261;
t217 = qJD(4) * t146;
t272 = -(t142 * t144 + t146 * t232) * t229 + t142 * t106 + t110 * t217;
t102 = t141 * t146 + t142 * t145;
t75 = t102 * t143;
t270 = -t142 * t218 + t147 * t214;
t269 = qJD(4) + qJD(5);
t211 = qJDD(1) * t140;
t190 = t143 * t211;
t213 = qJD(1) * qJD(2);
t72 = qJDD(2) * pkin(7) + (qJDD(1) * t144 + t148 * t213) * t139;
t18 = qJDD(3) * pkin(8) + qJD(3) * t271 + t147 * t72 + t190;
t194 = t144 * t213;
t178 = -qJDD(1) * t236 + t139 * t194;
t32 = qJD(2) * t106 + t110 * qJDD(2) + t178;
t31 = t146 * t32;
t155 = -t28 * qJD(4) - t142 * t18 + t31;
t131 = t147 * qJDD(2);
t95 = t143 * t212 + qJDD(4) - t131;
t2 = pkin(4) * t95 - pkin(9) * t37 + t155;
t209 = t142 * t32 + t146 * t18 + t70 * t217;
t219 = qJD(4) * t142;
t166 = t62 * t219 - t209;
t3 = -pkin(9) * t38 - t166;
t208 = -t141 * t3 + t145 * t2;
t251 = t14 * t145;
t27 = -t142 * t62 + t146 * t70;
t13 = -pkin(9) * t100 + t27;
t9 = -pkin(4) * t123 + t13;
t5 = t141 * t9 + t251;
t268 = -t33 * t174 - g(1) * (-t133 * t53 + t134 * t80) - g(2) * (-t133 * t51 + t134 * t78) - g(3) * (-t133 * t86 - t134 * t236) - t5 * qJD(5) + t208;
t7 = t174 * qJD(5) + t141 * t37 + t145 * t38;
t267 = -t119 * t174 - t7;
t149 = qJD(3) ^ 2;
t183 = g(1) * t80 + g(2) * t78;
t266 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t149 + t139 * (-g(3) * t148 + t194) - t178 + t183;
t195 = qJD(5) * t9 + t3;
t265 = t141 * t2 + t145 * t195;
t263 = -g(3) * t236 + (pkin(7) * t123 + t62) * qJD(4) + t183;
t262 = pkin(8) + pkin(9);
t233 = t146 * t147;
t124 = pkin(7) * t233;
t172 = pkin(4) * t143 - pkin(9) * t233;
t259 = -t172 * qJD(3) - (-t124 + (pkin(9) * t143 - t110) * t142) * qJD(4) + t273;
t221 = qJD(3) * t147;
t200 = t142 * t221;
t160 = t143 * t217 + t200;
t258 = -t160 * pkin(9) + (-t143 * t214 - t147 * t219) * pkin(7) + t272;
t101 = t141 * t142 - t145 * t146;
t162 = t101 * t147;
t257 = qJD(2) * t162 - t269 * t101;
t256 = (-t224 + t269) * t102;
t103 = t181 * qJD(2);
t255 = t142 * t103 + t146 * t271;
t253 = qJD(2) * pkin(2);
t252 = t123 * t98;
t250 = t37 * t142;
t248 = t142 * t110 + t124;
t246 = qJD(3) * t98;
t244 = t100 * t123;
t243 = t100 * t146;
t242 = t133 * t147;
t241 = t134 * t147;
t238 = t139 * t144;
t235 = t142 * t143;
t234 = t143 * t146;
t231 = qJDD(1) - g(3);
t135 = t143 ^ 2;
t230 = -t147 ^ 2 + t135;
t226 = qJD(2) * t139;
t220 = qJD(4) * t123;
t207 = qJD(4) * t262;
t205 = t143 * t227;
t204 = t144 * t226;
t203 = t148 * t226;
t202 = t142 * t224;
t201 = t123 * t214;
t192 = t148 * t212;
t184 = -t68 + (-t202 + t219) * pkin(4);
t182 = g(1) * t81 + g(2) * t79;
t114 = t262 * t142;
t180 = pkin(9) * t202 - qJD(5) * t114 - t142 * t207 - t255;
t115 = t262 * t146;
t89 = t146 * t103;
t179 = t172 * qJD(2) + qJD(5) * t115 - t142 * t271 + t146 * t207 + t89;
t97 = t146 * t110;
t46 = -pkin(9) * t234 + t97 + (-pkin(4) - t261) * t147;
t58 = -pkin(9) * t235 + t248;
t177 = t141 * t46 + t145 * t58;
t170 = t142 * t236 - t146 * t86;
t56 = -t142 * t86 - t146 * t236;
t176 = t141 * t170 + t145 * t56;
t175 = t141 * t56 - t145 * t170;
t150 = qJD(2) ^ 2;
t173 = qJDD(2) * t148 - t144 * t150;
t168 = -t123 * t217 + t142 * t95;
t167 = t123 * t219 + t146 * t95;
t85 = -t140 * t147 + t143 * t238;
t164 = g(1) * (t138 * t237 - t143 * t81) + g(2) * (-t79 * t143 - t147 * t188) - g(3) * t85;
t163 = qJD(3) * t121 + t108 * t221 + t143 * t72 - t147 * t211;
t159 = -g(3) * t238 - t182;
t157 = -pkin(8) * t95 - t123 * t61;
t19 = -qJDD(3) * pkin(3) + t163;
t153 = pkin(8) * t220 - t164 - t19;
t109 = -t197 - t253;
t152 = -pkin(7) * qJDD(3) + (t109 + t197 - t253) * qJD(3);
t128 = -pkin(4) * t146 - pkin(3);
t107 = (pkin(4) * t142 + pkin(7)) * t143;
t91 = qJDD(5) + t95;
t76 = t101 * t143;
t69 = t160 * pkin(4) + pkin(7) * t221;
t55 = t86 * qJD(3) + t143 * t203;
t54 = -t85 * qJD(3) + t147 * t203;
t21 = -t216 * t235 + (t269 * t234 + t200) * t145 + t270 * t141;
t20 = -qJD(3) * t162 - t269 * t75;
t12 = t56 * qJD(4) + t142 * t204 + t54 * t146;
t11 = t170 * qJD(4) - t54 * t142 + t146 * t204;
t8 = pkin(4) * t38 + t19;
t4 = -t14 * t141 + t145 * t9;
t1 = [t231, 0, t173 * t139, (-qJDD(2) * t144 - t148 * t150) * t139, 0, 0, 0, 0, 0, -t55 * qJD(3) - t85 * qJDD(3) + (-t143 * t192 + t147 * t173) * t139, -t54 * qJD(3) - t86 * qJDD(3) + (-t143 * t173 - t147 * t192) * t139, 0, 0, 0, 0, 0, -t11 * t123 + t38 * t85 + t55 * t98 + t56 * t95, t100 * t55 + t12 * t123 + t170 * t95 + t37 * t85, 0, 0, 0, 0, 0, t55 * t40 + t85 * t7 - (-qJD(5) * t175 + t145 * t11 - t141 * t12) * t119 + t176 * t91, t55 * t174 + t85 * t6 + (qJD(5) * t176 + t141 * t11 + t145 * t12) * t119 - t175 * t91; 0, qJDD(2), t231 * t236 + t183, -t231 * t238 + t182, qJDD(2) * t135 + 0.2e1 * t143 * t193, 0.2e1 * t143 * t131 - 0.2e1 * t230 * t212, qJDD(3) * t143 + t147 * t149, qJDD(3) * t147 - t143 * t149, 0, t152 * t143 + t266 * t147, -t266 * t143 + t152 * t147, t270 * t100 + t37 * t234, (-t100 * t142 - t146 * t98) * t221 + (-t250 - t146 * t38 + (t142 * t98 - t243) * qJD(4)) * t143, (-t37 - t201) * t147 + (qJD(3) * t100 + t167) * t143, (t123 * t223 + t38) * t147 + (-t168 - t246) * t143, -t123 * t222 - t147 * t95, t97 * t95 + t273 * t123 + (t110 * t220 + t159) * t142 + (pkin(7) * t246 - t31 + (-pkin(7) * t95 + qJD(3) * t61 + qJD(4) * t70 + t18) * t142 + t263 * t146) * t147 + (pkin(7) * t38 + qJD(3) * t27 + t19 * t142 - t197 * t98 + t217 * t61) * t143, -t248 * t95 + t272 * t123 + t159 * t146 + ((pkin(7) * t100 + t146 * t61) * qJD(3) - t263 * t142 + t209) * t147 + (-t100 * t197 - t61 * t219 - t28 * qJD(3) + t19 * t146 + (t37 - t201) * pkin(7)) * t143, t174 * t20 - t6 * t76, -t174 * t21 - t20 * t40 - t6 * t75 + t7 * t76, -t119 * t20 - t147 * t6 + t174 * t222 - t76 * t91, t119 * t21 + t147 * t7 - t40 * t222 - t75 * t91, -t119 * t222 - t147 * t91, t69 * t40 + t107 * t7 + t8 * t75 + t33 * t21 + (-t141 * t58 + t145 * t46) * t91 - t208 * t147 + t4 * t222 - g(1) * (t133 * t81 - t80 * t241) - g(2) * (t133 * t79 - t78 * t241) + (t258 * t141 + t259 * t145) * t119 + (t119 * t177 + t147 * t5) * qJD(5) + (-t40 * t205 - g(3) * (t133 * t144 + t134 * t232)) * t139, t69 * t174 + t107 * t6 - t8 * t76 + t33 * t20 - t177 * t91 + (-t10 + t265) * t147 - t5 * t222 - g(1) * (t134 * t81 + t80 * t242) - g(2) * (t134 * t79 + t78 * t242) + ((qJD(5) * t46 + t258) * t145 + (-qJD(5) * t58 - t259) * t141) * t119 + (-t174 * t205 - g(3) * (-t133 * t232 + t134 * t144)) * t139; 0, 0, 0, 0, -t143 * t150 * t147, t230 * t150, t210, t131, qJDD(3), qJD(3) * t68 - t109 * t225 - t163 - t164, -t190 + g(1) * t53 + g(2) * t51 + g(3) * t86 + (-qJD(2) * t109 - t72) * t147, -t123 * t243 + t250, (t37 + t252) * t146 + (-t38 + t244) * t142, (-t100 * t143 + t123 * t233) * qJD(2) + t168, (-t123 * t142 * t147 + t143 * t98) * qJD(2) + t167, t123 * t225, -t27 * t225 - pkin(3) * t38 + t89 * t123 - t68 * t98 + (-t123 * t271 + t157) * t142 + t153 * t146, -pkin(3) * t37 - t68 * t100 - t255 * t123 - t153 * t142 + t157 * t146 + t28 * t225, t6 * t102 + t174 * t257, -t101 * t6 - t102 * t7 - t174 * t256 - t257 * t40, t102 * t91 - t257 * t119 - t174 * t225, -t101 * t91 + t256 * t119 + t40 * t225, t119 * t225, t128 * t7 + t8 * t101 + (-t114 * t145 - t115 * t141) * t91 - t4 * t225 + t184 * t40 + t256 * t33 + (t141 * t180 + t145 * t179) * t119 - t164 * t134, t128 * t6 + t8 * t102 - (-t114 * t141 + t115 * t145) * t91 + t5 * t225 + t184 * t174 + t257 * t33 + (-t141 * t179 + t145 * t180) * t119 + t164 * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t98, t100 ^ 2 - t98 ^ 2, t37 - t252, -t244 - t38, t95, -t28 * t123 - t61 * t100 - g(1) * (-t142 * t53 + t146 * t80) - g(2) * (-t142 * t51 + t146 * t78) - g(3) * t56 + t155, -t27 * t123 + t61 * t98 - g(1) * (-t142 * t80 - t146 * t53) - g(2) * (-t142 * t78 - t146 * t51) - g(3) * t170 + t166, t279, t277, t276, t267, t91, (-t13 * t141 - t251) * t119 + (-t100 * t40 + t119 * t216 + t145 * t91) * pkin(4) + t268, (t119 * t14 - t2) * t141 + (-t119 * t13 - t195) * t145 + (-t100 * t174 + t119 * t215 - t141 * t91) * pkin(4) + t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t277, t276, t267, t91, -t119 * t5 + t268, -t119 * t4 - t265 + t275;];
tau_reg = t1;
