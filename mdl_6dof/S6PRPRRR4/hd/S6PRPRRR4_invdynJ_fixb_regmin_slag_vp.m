% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR4
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x29]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:39:04
% EndTime: 2019-03-08 20:39:15
% DurationCPUTime: 5.13s
% Computational Cost: add. (4743->456), mult. (11247->641), div. (0->0), fcn. (9580->18), ass. (0->241)
t190 = cos(pkin(12));
t323 = cos(qJ(4));
t258 = t323 * t190;
t170 = qJD(2) * t258;
t187 = sin(pkin(12));
t194 = sin(qJ(4));
t284 = t194 * t187;
t255 = qJD(2) * t284;
t141 = -t170 + t255;
t330 = qJD(5) + qJD(6);
t344 = t141 + t330;
t217 = t258 - t284;
t189 = sin(pkin(6));
t198 = cos(qJ(2));
t287 = t189 * t198;
t206 = t217 * t287;
t120 = qJD(1) * t206;
t320 = pkin(8) + qJ(3);
t155 = t320 * t187;
t156 = t320 * t190;
t218 = -t323 * t155 - t194 * t156;
t64 = t217 * qJD(3) + t218 * qJD(4);
t343 = t120 - t64;
t144 = t217 * qJD(4);
t148 = t323 * t187 + t194 * t190;
t145 = t148 * qJD(4);
t195 = sin(qJ(2));
t279 = qJD(1) * t189;
t257 = t195 * t279;
t342 = -pkin(4) * t145 + pkin(9) * t144 + t257;
t143 = t148 * qJD(2);
t193 = sin(qJ(5));
t197 = cos(qJ(5));
t272 = t197 * qJD(4);
t121 = t143 * t193 - t272;
t123 = qJD(4) * t193 + t143 * t197;
t192 = sin(qJ(6));
t196 = cos(qJ(6));
t231 = t121 * t192 - t196 * t123;
t53 = t196 * t121 + t123 * t192;
t341 = t231 * t53;
t151 = t192 * t197 + t193 * t196;
t308 = t344 * t151;
t154 = qJD(2) * qJ(3) + t257;
t191 = cos(pkin(6));
t278 = qJD(1) * t191;
t166 = t190 * t278;
t317 = pkin(8) * qJD(2);
t110 = t166 + (-t154 - t317) * t187;
t126 = t190 * t154 + t187 * t278;
t111 = t190 * t317 + t126;
t45 = t194 * t110 + t323 * t111;
t340 = t45 * qJD(4);
t276 = qJD(5) * t193;
t299 = t141 * t193;
t339 = t276 + t299;
t338 = t231 ^ 2 - t53 ^ 2;
t134 = qJD(5) + t141;
t132 = qJD(6) + t134;
t273 = qJD(6) * t196;
t274 = qJD(6) * t192;
t249 = qJDD(2) * t323;
t267 = t190 * qJDD(2);
t262 = qJD(4) * t170 + t187 * t249 + t194 * t267;
t98 = -qJD(4) * t255 + t262;
t38 = qJD(5) * t272 + t193 * qJDD(4) - t143 * t276 + t197 * t98;
t39 = t123 * qJD(5) - t197 * qJDD(4) + t193 * t98;
t8 = -t121 * t273 - t123 * t274 - t192 * t39 + t196 * t38;
t337 = t132 * t53 + t8;
t305 = cos(pkin(11));
t247 = t305 * t195;
t188 = sin(pkin(11));
t289 = t188 * t198;
t138 = t191 * t247 + t289;
t185 = pkin(12) + qJ(4);
t176 = sin(t185);
t177 = cos(t185);
t248 = t189 * t305;
t101 = t138 * t177 - t176 * t248;
t246 = t305 * t198;
t290 = t188 * t195;
t140 = -t191 * t290 + t246;
t291 = t188 * t189;
t103 = t140 * t177 + t176 * t291;
t288 = t189 * t195;
t131 = t176 * t191 + t177 * t288;
t137 = -t191 * t246 + t290;
t139 = t191 * t289 + t247;
t41 = qJD(4) * pkin(9) + t45;
t172 = -pkin(3) * t190 - pkin(2);
t256 = t198 * t279;
t236 = qJD(3) - t256;
t133 = t172 * qJD(2) + t236;
t57 = pkin(4) * t141 - pkin(9) * t143 + t133;
t22 = t193 * t57 + t197 * t41;
t17 = -pkin(10) * t121 + t22;
t15 = t17 * t274;
t186 = qJ(5) + qJ(6);
t181 = sin(t186);
t182 = cos(t186);
t331 = t323 * t110 - t194 * t111;
t40 = -qJD(4) * pkin(4) - t331;
t30 = t121 * pkin(5) + t40;
t336 = t30 * t53 - g(1) * (-t103 * t182 - t139 * t181) - g(2) * (-t101 * t182 - t137 * t181) - g(3) * (-t131 * t182 + t181 * t287) + t15;
t269 = qJDD(2) * qJ(3);
t271 = qJDD(1) * t189;
t124 = t195 * t271 + t269 + (qJD(3) + t256) * qJD(2);
t270 = qJDD(1) * t191;
t164 = t190 * t270;
t69 = t164 + (-pkin(8) * qJDD(2) - t124) * t187;
t88 = t190 * t124 + t187 * t270;
t70 = pkin(8) * t267 + t88;
t227 = t194 * t69 + t323 * t70;
t13 = qJDD(4) * pkin(9) + qJD(4) * t331 + t227;
t277 = qJD(2) * t195;
t252 = qJD(1) * t277;
t266 = t189 * t252 + qJDD(3);
t223 = -t198 * t271 + t266;
t113 = t172 * qJDD(2) + t223;
t268 = t187 * qJDD(2);
t235 = -t190 * t249 + t194 * t268;
t99 = qJD(2) * t145 + t235;
t29 = pkin(4) * t99 - pkin(9) * t98 + t113;
t27 = t197 * t29;
t203 = -t22 * qJD(5) - t193 * t13 + t27;
t92 = qJDD(5) + t99;
t2 = pkin(5) * t92 - pkin(10) * t38 + t203;
t275 = qJD(5) * t197;
t222 = -t197 * t13 - t193 * t29 - t57 * t275 + t41 * t276;
t3 = -pkin(10) * t39 - t222;
t261 = -t192 * t3 + t196 * t2;
t21 = -t193 * t41 + t197 * t57;
t16 = -pkin(10) * t123 + t21;
t11 = pkin(5) * t134 + t16;
t313 = t17 * t196;
t5 = t11 * t192 + t313;
t335 = t30 * t231 - g(1) * (-t103 * t181 + t139 * t182) - g(2) * (-t101 * t181 + t137 * t182) - g(3) * (-t131 * t181 - t182 * t287) - t5 * qJD(6) + t261;
t202 = t231 * qJD(6) - t192 * t38 - t196 * t39;
t334 = -t132 * t231 + t202;
t333 = t120 * t193 - t197 * t342;
t118 = -t194 * t155 + t323 * t156;
t95 = -pkin(4) * t217 - pkin(9) * t148 + t172;
t332 = t118 * t276 + t342 * t193 + t197 * t343 - t95 * t275;
t207 = t148 * t287;
t306 = -qJD(1) * t207 + t148 * qJD(3) + t118 * qJD(4);
t81 = t151 * t148;
t199 = qJD(2) ^ 2;
t212 = (qJDD(2) * t198 - t195 * t199) * t189;
t283 = t197 * t144;
t215 = -t148 * t276 + t283;
t150 = t192 * t193 - t196 * t197;
t309 = t344 * t150;
t86 = qJDD(6) + t92;
t329 = t132 * t309 - t151 * t86;
t304 = qJDD(2) * pkin(2);
t127 = t223 - t304;
t240 = g(1) * t139 + g(2) * t137;
t211 = t240 + t304;
t328 = (-g(3) * t198 + t252) * t189 - t127 + t211;
t251 = qJD(6) * t11 + t3;
t327 = t192 * t2 + t196 * t251;
t230 = (-t154 * t187 + t166) * t187 - t126 * t190;
t326 = t198 * t230 - (-qJD(2) * pkin(2) + t236) * t195;
t325 = pkin(9) + pkin(10);
t107 = t197 * t118;
t324 = -pkin(10) * t283 + pkin(5) * t145 - t193 * t64 + (-t107 + (pkin(10) * t148 - t95) * t193) * qJD(5) + t333;
t322 = g(3) * t189;
t286 = t193 * t144;
t216 = t148 * t275 + t286;
t319 = t216 * pkin(10) + t332;
t93 = pkin(4) * t143 + pkin(9) * t141;
t318 = t193 * t93 + t197 * t331;
t316 = t143 * t53;
t315 = t143 * t231;
t312 = t193 * t38;
t311 = t193 * t92;
t310 = t193 * t95 + t107;
t307 = t216 * pkin(5) + t306;
t303 = t121 * t134;
t302 = t121 * t143;
t301 = t123 * t134;
t300 = t123 * t143;
t298 = t148 * t193;
t297 = t148 * t197;
t295 = t177 * t181;
t294 = t177 * t182;
t293 = t177 * t193;
t292 = t177 * t198;
t285 = t193 * t198;
t282 = t198 * t199;
t281 = qJDD(1) - g(3);
t280 = t187 ^ 2 + t190 ^ 2;
t265 = g(3) * t288;
t264 = t189 * t285;
t263 = t197 * t287;
t260 = qJD(5) * t325;
t254 = t189 * t277;
t244 = t189 * t281;
t243 = t134 * t197;
t242 = -t132 * t308 - t150 * t86;
t241 = pkin(5) * t339 - t45;
t239 = g(1) * t140 + g(2) * t138;
t159 = t325 * t193;
t238 = pkin(10) * t299 + qJD(6) * t159 + t193 * t260 + t318;
t160 = t325 * t197;
t79 = t197 * t93;
t237 = pkin(5) * t143 + qJD(6) * t160 - t193 * t331 + t79 + (pkin(10) * t141 + t260) * t197;
t84 = t197 * t95;
t28 = -pkin(5) * t217 - pkin(10) * t297 - t118 * t193 + t84;
t31 = -pkin(10) * t298 + t310;
t234 = t192 * t28 + t196 * t31;
t135 = -t187 * t288 + t190 * t191;
t136 = t187 * t191 + t190 * t288;
t74 = t194 * t135 + t323 * t136;
t226 = -t197 * t74 + t264;
t58 = -t193 * t74 - t263;
t233 = t192 * t226 + t196 * t58;
t232 = t192 * t58 - t196 * t226;
t228 = -t134 * t339 + t197 * t92;
t225 = t194 * t70 - t323 * t69 + t340;
t221 = -pkin(9) * t92 + t134 * t40;
t219 = t323 * t135 - t194 * t136;
t213 = g(1) * (-t140 * t176 + t177 * t291) + g(2) * (-t138 * t176 - t177 * t248) + g(3) * (-t176 * t288 + t177 * t191);
t209 = g(3) * t287 - t240;
t14 = -qJDD(4) * pkin(4) + t225;
t205 = t209 * t177;
t87 = -t124 * t187 + t164;
t204 = -t187 * t87 + t190 * t88 - t239;
t201 = pkin(9) * qJD(5) * t134 + t14 + t213;
t175 = -pkin(5) * t197 - pkin(4);
t82 = t150 * t148;
t66 = pkin(5) * t298 - t218;
t43 = qJD(2) * t207 + t74 * qJD(4);
t42 = qJD(2) * t206 + t219 * qJD(4);
t25 = -t274 * t298 + (t297 * t330 + t286) * t196 + t215 * t192;
t24 = -t150 * t144 - t330 * t81;
t20 = t226 * qJD(5) - t193 * t42 + t197 * t254;
t19 = t58 * qJD(5) + t193 * t254 + t197 * t42;
t6 = t39 * pkin(5) + t14;
t4 = t11 * t196 - t17 * t192;
t1 = [t281, 0, t212 (-qJDD(2) * t195 - t282) * t189, t190 * t212, -t187 * t212, t280 * t189 * t282 + (-t135 * t187 + t136 * t190) * qJDD(2), t135 * t87 + t136 * t88 - g(3) + (-qJD(2) * t326 - t127 * t198) * t189, 0, 0, 0, 0, 0, -qJD(4) * t43 + qJDD(4) * t219 + (t141 * t277 - t198 * t99) * t189, -qJD(4) * t42 - qJDD(4) * t74 + (t143 * t277 - t198 * t98) * t189, 0, 0, 0, 0, 0, t121 * t43 + t134 * t20 - t219 * t39 + t58 * t92, t123 * t43 - t134 * t19 - t219 * t38 + t226 * t92, 0, 0, 0, 0, 0 (-qJD(6) * t232 - t19 * t192 + t196 * t20) * t132 + t233 * t86 + t43 * t53 + t219 * t202 -(qJD(6) * t233 + t19 * t196 + t192 * t20) * t132 - t232 * t86 - t43 * t231 - t219 * t8; 0, qJDD(2), t281 * t287 + t240, -t195 * t244 + t239, t328 * t190, -t328 * t187, -t265 + t204 + (qJD(2) * t236 + t269) * t280, -t230 * qJD(3) + (-t127 + t240) * pkin(2) + t204 * qJ(3) + (-g(3) * (pkin(2) * t198 + qJ(3) * t195) + t326 * qJD(1)) * t189, t143 * t144 + t148 * t98, -t141 * t144 - t143 * t145 - t148 * t99 + t217 * t98, qJD(4) * t144 + qJDD(4) * t148, -qJD(4) * t145 + qJDD(4) * t217, 0, -t306 * qJD(4) + qJDD(4) * t218 - t113 * t217 + t133 * t145 - t141 * t257 + t172 * t99 - t205, qJD(4) * t343 - qJDD(4) * t118 + t113 * t148 + t133 * t144 - t143 * t257 + t172 * t98 + t209 * t176, t123 * t215 + t38 * t297 (-t121 * t197 - t123 * t193) * t144 + (-t312 - t197 * t39 + (t121 * t193 - t123 * t197) * qJD(5)) * t148, t123 * t145 + t134 * t215 - t217 * t38 + t92 * t297, -t121 * t145 - t134 * t216 + t217 * t39 - t92 * t298, t134 * t145 - t217 * t92, -t218 * t39 + t21 * t145 - t27 * t217 + t84 * t92 + t333 * t134 + t306 * t121 + (-t205 + (-t118 * t134 + t148 * t40 + t217 * t41) * qJD(5)) * t197 + ((-qJD(5) * t95 - t64) * t134 - t118 * t92 - (-qJD(5) * t57 - t13) * t217 + t14 * t148 + t40 * t144 - t265 - t239) * t193, -t310 * t92 - t222 * t217 - t22 * t145 - t218 * t38 + t40 * t283 - g(1) * (t139 * t293 + t140 * t197) - g(2) * (t137 * t293 + t138 * t197) - (-t177 * t285 + t195 * t197) * t322 + (t14 * t197 - t276 * t40) * t148 + t332 * t134 + t306 * t123, -t231 * t24 - t8 * t82, -t202 * t82 + t231 * t25 - t24 * t53 - t8 * t81, t132 * t24 - t145 * t231 - t217 * t8 - t82 * t86, -t132 * t25 - t145 * t53 - t202 * t217 - t81 * t86, t132 * t145 - t217 * t86 (-t192 * t31 + t196 * t28) * t86 - t261 * t217 + t4 * t145 - t66 * t202 + t6 * t81 + t30 * t25 - g(1) * (-t139 * t294 + t140 * t181) - g(2) * (-t137 * t294 + t138 * t181) + t307 * t53 - (t181 * t195 + t182 * t292) * t322 + (t319 * t192 + t324 * t196) * t132 + (-t132 * t234 + t217 * t5) * qJD(6), -t234 * t86 + (-t15 + t327) * t217 - t5 * t145 + t66 * t8 - t6 * t82 + t30 * t24 - g(1) * (t139 * t295 + t140 * t182) - g(2) * (t137 * t295 + t138 * t182) - t307 * t231 - (-t181 * t292 + t182 * t195) * t322 + ((-qJD(6) * t28 + t319) * t196 + (qJD(6) * t31 - t324) * t192) * t132; 0, 0, 0, 0, -t267, t268, -t280 * t199, t230 * qJD(2) - t198 * t244 - t211 + t266, 0, 0, 0, 0, 0, 0.2e1 * t143 * qJD(4) + t235 (-t141 - t255) * qJD(4) + t262, 0, 0, 0, 0, 0, t228 - t302, -t134 ^ 2 * t197 - t300 - t311, 0, 0, 0, 0, 0, t242 - t316, t315 + t329; 0, 0, 0, 0, 0, 0, 0, 0, t143 * t141, -t141 ^ 2 + t143 ^ 2 (t141 - t255) * qJD(4) + t262, -t235, qJDD(4), -t133 * t143 - t213 - t225 + t340, g(1) * t103 + g(2) * t101 + g(3) * t131 + t133 * t141 - t227, t123 * t243 + t312 (t38 - t303) * t197 + (-t39 - t301) * t193, t134 * t243 - t300 + t311, t228 + t302, -t134 * t143, -pkin(4) * t39 - t45 * t121 - t79 * t134 - t21 * t143 + (t134 * t331 + t221) * t193 - t201 * t197, -pkin(4) * t38 - t45 * t123 + t318 * t134 + t22 * t143 + t201 * t193 + t221 * t197, t151 * t8 + t231 * t309, -t150 * t8 + t151 * t202 + t231 * t308 + t309 * t53, t315 - t329, t242 + t316, -t132 * t143 (-t159 * t196 - t160 * t192) * t86 - t175 * t202 + t6 * t150 - t4 * t143 + t241 * t53 + t308 * t30 + (t192 * t238 - t196 * t237) * t132 - t213 * t182 -(-t159 * t192 + t160 * t196) * t86 + t175 * t8 + t6 * t151 + t5 * t143 - t241 * t231 - t309 * t30 + (t192 * t237 + t196 * t238) * t132 + t213 * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t121, -t121 ^ 2 + t123 ^ 2, t38 + t303, t301 - t39, t92, t22 * t134 - t40 * t123 - g(1) * (-t103 * t193 + t139 * t197) - g(2) * (-t101 * t193 + t137 * t197) - g(3) * (-t131 * t193 - t263) + t203, t21 * t134 + t40 * t121 - g(1) * (-t103 * t197 - t139 * t193) - g(2) * (-t101 * t197 - t137 * t193) - g(3) * (-t131 * t197 + t264) + t222, -t341, t338, t337, t334, t86 -(-t16 * t192 - t313) * t132 + (-t123 * t53 - t132 * t274 + t196 * t86) * pkin(5) + t335 (-t132 * t17 - t2) * t192 + (t132 * t16 - t251) * t196 + (t123 * t231 - t132 * t273 - t192 * t86) * pkin(5) + t336; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t341, t338, t337, t334, t86, t132 * t5 + t335, t132 * t4 - t327 + t336;];
tau_reg  = t1;
