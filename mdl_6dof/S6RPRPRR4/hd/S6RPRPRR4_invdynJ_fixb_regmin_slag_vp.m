% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:46:30
% EndTime: 2019-03-09 03:46:39
% DurationCPUTime: 3.52s
% Computational Cost: add. (2984->393), mult. (5918->520), div. (0->0), fcn. (3878->14), ass. (0->232)
t169 = sin(pkin(10));
t143 = t169 * pkin(1) + pkin(7);
t126 = t143 * qJDD(1);
t319 = qJD(2) * qJD(3) + t126;
t173 = sin(qJ(3));
t250 = qJD(1) * qJD(3);
t233 = t173 * t250;
t177 = cos(qJ(3));
t245 = t177 * qJDD(1);
t317 = t233 - t245;
t129 = t143 * qJD(1);
t267 = t177 * qJD(2) - t173 * t129;
t315 = qJD(4) - t267;
t298 = pkin(4) + t143;
t171 = sin(qJ(6));
t172 = sin(qJ(5));
t175 = cos(qJ(6));
t176 = cos(qJ(5));
t119 = t171 * t176 + t175 * t172;
t201 = t119 * t173;
t308 = qJD(5) + qJD(6);
t294 = -qJD(1) * t201 - t308 * t119;
t254 = t172 * qJD(3);
t264 = qJD(1) * t177;
t116 = t176 * t264 + t254;
t236 = t172 * t264;
t251 = t176 * qJD(3);
t118 = -t236 + t251;
t215 = t171 * t116 - t175 * t118;
t48 = t175 * t116 + t171 * t118;
t318 = t215 * t48;
t69 = -qJD(3) * pkin(3) + t315;
t258 = qJD(5) * t177;
t237 = t172 * t258;
t316 = t173 * t251 + t237;
t314 = t215 ^ 2 - t48 ^ 2;
t253 = t173 * qJD(1);
t142 = qJD(5) + t253;
t136 = qJD(6) + t142;
t255 = qJD(6) * t175;
t256 = qJD(6) * t171;
t41 = -qJD(5) * t116 + t176 * qJDD(3) + t317 * t172;
t42 = -qJD(5) * t236 + t172 * qJDD(3) + (qJD(3) * qJD(5) - t317) * t176;
t7 = -t116 * t255 - t118 * t256 - t171 * t42 + t175 * t41;
t313 = t136 * t48 + t7;
t179 = -pkin(3) - pkin(8);
t269 = pkin(4) * t253 + t315;
t53 = t179 * qJD(3) + t269;
t170 = cos(pkin(10));
t299 = t170 * pkin(1);
t144 = -pkin(2) - t299;
t275 = t173 * qJ(4);
t84 = t179 * t177 + t144 - t275;
t60 = t84 * qJD(1);
t18 = t172 * t53 + t176 * t60;
t15 = -t116 * pkin(9) + t18;
t13 = t15 * t256;
t168 = qJ(5) + qJ(6);
t158 = sin(t168);
t161 = g(3) * t177;
t165 = qJD(3) * qJ(4);
t87 = t173 * qJD(2) + t177 * t129;
t64 = pkin(4) * t264 + t87;
t57 = t165 + t64;
t34 = t116 * pkin(5) + t57;
t162 = qJ(1) + pkin(10);
t152 = sin(t162);
t153 = cos(t162);
t159 = cos(t168);
t281 = t158 * t173;
t66 = t152 * t159 + t153 * t281;
t68 = -t152 * t281 + t153 * t159;
t312 = g(1) * t66 - g(2) * t68 - t158 * t161 + t34 * t48 + t13;
t232 = t177 * t250;
t246 = t173 * qJDD(1);
t200 = t232 + t246;
t111 = qJDD(5) + t200;
t244 = t177 * qJDD(2);
t261 = qJD(3) * t177;
t307 = -t129 * t261 - t319 * t173;
t206 = t244 + t307;
t198 = qJDD(4) - t206;
t24 = t200 * pkin(4) + t179 * qJDD(3) + t198;
t141 = pkin(3) * t233;
t286 = qJ(4) * t177;
t218 = pkin(8) * t173 - t286;
t252 = t173 * qJD(4);
t191 = t218 * qJD(3) - t252;
t28 = qJD(1) * t191 + qJDD(1) * t84 + t141;
t229 = -t172 * t28 + t176 * t24;
t189 = -t18 * qJD(5) + t229;
t2 = t111 * pkin(5) - t41 * pkin(9) + t189;
t259 = qJD(5) * t176;
t243 = -t172 * t24 - t176 * t28 - t53 * t259;
t260 = qJD(5) * t172;
t3 = -t42 * pkin(9) - t60 * t260 - t243;
t239 = -t171 * t3 + t175 * t2;
t17 = -t172 * t60 + t176 * t53;
t14 = -t118 * pkin(9) + t17;
t12 = t142 * pkin(5) + t14;
t289 = t175 * t15;
t5 = t171 * t12 + t289;
t280 = t159 * t173;
t65 = -t152 * t158 + t153 * t280;
t67 = t152 * t280 + t153 * t158;
t311 = -g(1) * t65 - g(2) * t67 - t5 * qJD(6) + t159 * t161 + t34 * t215 + t239;
t188 = t215 * qJD(6) - t171 * t41 - t175 * t42;
t310 = -t136 * t215 + t188;
t72 = -t165 - t87;
t263 = qJD(3) * t116;
t98 = t176 * t111;
t309 = t263 - t98;
t224 = g(1) * t153 + g(2) * t152;
t306 = t179 * t111 + t142 * t57;
t196 = t173 * t254 - t176 * t258;
t276 = t172 * t177;
t305 = t111 * t276 - t142 * t196;
t304 = t177 * t308;
t303 = t7 * t173 - t215 * t261;
t300 = g(3) * t173;
t297 = pkin(9) - t179;
t101 = qJDD(6) + t111;
t273 = t175 * t176;
t279 = t171 * t172;
t214 = -t273 + t279;
t262 = qJD(3) * t173;
t23 = t119 * t304 - t214 * t262;
t89 = t214 * t177;
t296 = t89 * t101 + t23 * t136;
t295 = t118 * t261 + t41 * t173;
t238 = t176 * t253;
t293 = -t171 * t260 - t172 * t256 + t175 * t238 - t253 * t279 + t308 * t273;
t103 = t298 * t173;
t91 = t172 * t103;
t292 = t176 * t84 + t91;
t151 = pkin(3) * t253;
t88 = t218 * qJD(1) + t151;
t291 = t172 * t64 + t176 * t88;
t290 = t173 * t42;
t288 = t41 * t176;
t240 = -pkin(5) * t176 - pkin(4);
t287 = pkin(5) * t259 - t240 * t253 + t315;
t285 = qJD(5) * t60;
t284 = qJDD(3) * pkin(3);
t283 = t116 * t142;
t282 = t118 * t142;
t278 = t172 * t111;
t277 = t172 * t173;
t274 = t173 * t176;
t272 = t176 * t177;
t270 = t87 * qJD(3);
t268 = t316 * t142;
t104 = t298 * t177;
t166 = t173 ^ 2;
t167 = t177 ^ 2;
t266 = t166 - t167;
t220 = t177 * pkin(3) + t275;
t210 = pkin(2) + t220;
t102 = -t210 - t299;
t73 = qJD(1) * t102;
t121 = -qJ(4) * t264 + t151;
t265 = qJD(1) * t121;
t130 = qJD(1) * t144;
t257 = qJD(5) * t179;
t248 = qJDD(1) * t102;
t247 = qJDD(3) * t143;
t181 = qJD(1) ^ 2;
t241 = t177 * t181 * t173;
t234 = pkin(9) * t177 - t84;
t123 = t297 * t176;
t231 = qJD(6) * t12 + t3;
t150 = pkin(3) * t262;
t70 = t150 + t191;
t96 = t298 * t261;
t228 = -t172 * t70 + t176 * t96;
t226 = t173 * qJDD(2) - t129 * t262 + t319 * t177;
t225 = -t214 * t101 + t294 * t136;
t174 = sin(qJ(1));
t178 = cos(qJ(1));
t223 = g(1) * t174 - g(2) * t178;
t122 = t297 * t172;
t211 = pkin(5) * t177 - pkin(9) * t277;
t59 = t176 * t64;
t222 = t211 * qJD(1) - qJD(6) * t122 - t172 * t88 - t297 * t260 + t59;
t221 = pkin(9) * t238 + t308 * t123 + t291;
t219 = pkin(3) * t173 - t286;
t22 = qJD(3) * t201 + t214 * t304;
t90 = t119 * t177;
t217 = t101 * t90 - t136 * t22;
t213 = t142 * t172;
t208 = t223 * pkin(1);
t207 = t173 * t188 - t48 * t261;
t163 = qJDD(3) * qJ(4);
t164 = qJD(3) * qJD(4);
t31 = -t163 - t164 - t226;
t205 = t103 * t259 + t172 * t96 + t176 * t70 - t84 * t260;
t204 = -t119 * t101 - t293 * t136;
t203 = -qJ(4) * t261 - t252;
t202 = qJD(3) * t267 - t226;
t180 = qJD(3) ^ 2;
t199 = g(1) * t152 - g(2) * t153 - t143 * t180;
t197 = -qJD(1) * t130 + t224;
t194 = -0.2e1 * t73 * qJD(3) + t247;
t193 = 0.2e1 * t130 * qJD(3) - t247;
t192 = -t224 * t177 - t300;
t25 = -t317 * pkin(4) - t31;
t190 = t25 + t192;
t187 = -0.2e1 * qJDD(1) * t144 + t199;
t32 = t198 - t284;
t186 = t32 * t173 - t31 * t177 + (t173 * t72 + t177 * t69) * qJD(3);
t185 = -t224 * t173 + t73 * t253 + qJDD(4) + t161 - t307;
t33 = qJD(1) * t203 + t141 + t248;
t97 = t150 + t203;
t184 = -qJD(1) * t97 + t199 - t248 - t33;
t145 = t172 * pkin(5) + qJ(4);
t125 = qJDD(3) * t177 - t180 * t173;
t124 = qJDD(3) * t173 + t180 * t177;
t95 = t298 * t262;
t92 = t176 * t103;
t83 = -t152 * t277 + t153 * t176;
t82 = t152 * t274 + t153 * t172;
t81 = t152 * t176 + t153 * t277;
t80 = -t152 * t172 + t153 * t274;
t74 = pkin(5) * t272 + t104;
t45 = -pkin(5) * t237 + (-t143 + t240) * t262;
t30 = -pkin(9) * t272 + t292;
t29 = t173 * pkin(5) + t234 * t172 + t92;
t11 = t42 * pkin(5) + t25;
t10 = t316 * pkin(9) + t205;
t9 = t211 * qJD(3) + (t234 * t176 - t91) * qJD(5) + t228;
t4 = t175 * t12 - t171 * t15;
t1 = [qJDD(1), t223, g(1) * t178 + g(2) * t174 (t169 ^ 2 + t170 ^ 2) * pkin(1) ^ 2 * qJDD(1) + t208, t166 * qJDD(1) + 0.2e1 * t173 * t232, 0.2e1 * t173 * t245 - 0.2e1 * t266 * t250, t124, t125, 0, t173 * t193 + t177 * t187, -t173 * t187 + t177 * t193 (t166 + t167) * t126 + t186 - t224, t173 * t194 - t177 * t184, t173 * t184 + t177 * t194, t33 * t102 + t73 * t97 + t208 + (-g(1) * pkin(7) - g(2) * t210) * t153 + (-g(2) * pkin(7) + g(1) * t210) * t152 + t186 * t143, t118 * t196 - t41 * t276 (-t116 * t172 + t118 * t176) * t262 + (t172 * t42 - t288 + (t116 * t176 + t118 * t172) * qJD(5)) * t177, t295 - t305, -t290 + (-t263 - t98) * t177 + t268, t111 * t173 + t142 * t261, t228 * t142 + (-t172 * t84 + t92) * t111 + t229 * t173 - t95 * t116 + t104 * t42 + t25 * t272 - g(1) * t83 - g(2) * t81 + (t17 * t177 - t57 * t274) * qJD(3) + (-t292 * t142 - t18 * t173 - t57 * t276) * qJD(5), -t205 * t142 - t292 * t111 - t95 * t118 + t104 * t41 + g(1) * t82 - g(2) * t80 + ((qJD(3) * t57 + t285) * t172 + t243) * t173 + (-qJD(3) * t18 - t25 * t172 - t57 * t259) * t177, -t215 * t22 - t7 * t90, -t188 * t90 - t215 * t23 - t22 * t48 + t7 * t89, -t217 + t303, t207 + t296, t101 * t173 + t136 * t261 (-t171 * t10 + t175 * t9) * t136 + (-t171 * t30 + t175 * t29) * t101 + t239 * t173 + t4 * t261 + t45 * t48 - t74 * t188 - t11 * t89 - t34 * t23 - g(1) * t68 - g(2) * t66 + ((-t171 * t29 - t175 * t30) * t136 - t5 * t173) * qJD(6), -t5 * t261 + g(1) * t67 - g(2) * t65 - t11 * t90 + t13 * t173 + t34 * t22 - t45 * t215 + t74 * t7 + (-(-qJD(6) * t30 + t9) * t136 - t29 * t101 - t2 * t173) * t171 + (-(qJD(6) * t29 + t10) * t136 - t30 * t101 - t231 * t173) * t175; 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, t125, -t124, 0, -t125, t124, -t31 * t173 - t32 * t177 - g(3) + (t173 * t69 - t177 * t72) * qJD(3), 0, 0, 0, 0, 0, t309 * t177 + t268 + t290, t295 + t305, 0, 0, 0, 0, 0, -t207 + t296, t217 + t303; 0, 0, 0, 0, -t241, t266 * t181, t246, t245, qJDD(3), t173 * t197 - t161 + t206 + t270, t177 * t197 + t202 + t300, -t219 * qJDD(1), -0.2e1 * t284 - t270 + (-qJDD(2) - t265) * t177 + t185, 0.2e1 * t163 + 0.2e1 * t164 + (-g(3) + t265) * t173 + (qJD(1) * t73 - t224) * t177 - t202, -t32 * pkin(3) - g(3) * t220 - t31 * qJ(4) - t73 * t121 + t224 * t219 - t315 * t72 - t69 * t87, -t118 * t213 + t288 (-t42 - t282) * t176 + (-t41 + t283) * t172, -t142 * t260 + t98 + (-t118 * t177 - t142 * t277) * qJD(1), -t142 * t259 - t278 + (t116 * t177 - t142 * t274) * qJD(1), -t142 * t264, -t17 * t264 + qJ(4) * t42 - t59 * t142 + t269 * t116 + t306 * t176 + ((t88 - t257) * t142 + t190) * t172, qJ(4) * t41 + t291 * t142 + t18 * t264 + t269 * t118 - t306 * t172 + (-t142 * t257 + t190) * t176, -t214 * t7 - t215 * t294, -t7 * t119 - t188 * t214 + t215 * t293 - t294 * t48, t215 * t264 + t225, t264 * t48 + t204, -t136 * t264 (t122 * t171 - t123 * t175) * t101 - t145 * t188 + t11 * t119 - t4 * t264 + t287 * t48 + t293 * t34 + (t171 * t221 - t175 * t222) * t136 + t192 * t158 -(-t122 * t175 - t123 * t171) * t101 + t145 * t7 - t11 * t214 + t5 * t264 - t287 * t215 + t294 * t34 + (t171 * t222 + t175 * t221) * t136 + t192 * t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, qJDD(3) + t241, -t166 * t181 - t180, t72 * qJD(3) + t185 - t244 - t284, 0, 0, 0, 0, 0, -t142 * t213 - t309, -t142 ^ 2 * t176 - qJD(3) * t118 - t278, 0, 0, 0, 0, 0, -qJD(3) * t48 + t225, qJD(3) * t215 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118 * t116, -t116 ^ 2 + t118 ^ 2, t41 + t283, t282 - t42, t111, -g(1) * t80 - g(2) * t82 + g(3) * t272 - t57 * t118 + t18 * t142 + t189, g(1) * t81 - g(2) * t83 + t57 * t116 + t17 * t142 + (t285 - t161) * t172 + t243, -t318, t314, t313, t310, t101 -(-t14 * t171 - t289) * t136 + (t101 * t175 - t118 * t48 - t136 * t256) * pkin(5) + t311 (-t136 * t15 - t2) * t171 + (t136 * t14 - t231) * t175 + (-t101 * t171 + t118 * t215 - t136 * t255) * pkin(5) + t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t314, t313, t310, t101, t5 * t136 + t311, t4 * t136 - t171 * t2 - t175 * t231 + t312;];
tau_reg  = t1;
