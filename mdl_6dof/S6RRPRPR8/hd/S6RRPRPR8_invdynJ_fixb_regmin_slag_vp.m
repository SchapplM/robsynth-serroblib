% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPR8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPR8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:53:25
% EndTime: 2019-03-09 10:53:39
% DurationCPUTime: 7.09s
% Computational Cost: add. (7121->578), mult. (16008->735), div. (0->0), fcn. (11912->12), ass. (0->270)
t241 = sin(pkin(10));
t383 = pkin(8) + qJ(3);
t197 = t383 * t241;
t242 = cos(pkin(10));
t198 = t383 * t242;
t245 = sin(qJ(4));
t393 = cos(qJ(4));
t134 = -t245 * t197 + t393 * t198;
t249 = cos(qJ(2));
t231 = t249 * qJDD(1);
t246 = sin(qJ(2));
t340 = qJD(1) * qJD(2);
t275 = t246 * t340 - t231;
t185 = qJDD(4) + t275;
t238 = pkin(10) + qJ(4);
t229 = sin(t238);
t237 = g(3) * t249;
t247 = sin(qJ(1));
t250 = cos(qJ(1));
t312 = g(1) * t250 + g(2) * t247;
t284 = t312 * t246;
t413 = t284 - t237;
t422 = t134 * t185 + t229 * t413;
t244 = sin(qJ(6));
t248 = cos(qJ(6));
t350 = qJD(1) * t246;
t328 = t241 * t350;
t347 = qJD(2) * t242;
t180 = -t328 + t347;
t331 = t242 * t350;
t348 = qJD(2) * t241;
t181 = t331 + t348;
t325 = t249 * t340;
t338 = t246 * qJDD(1);
t276 = t325 + t338;
t339 = qJDD(2) * t241;
t259 = t242 * t276 + t339;
t354 = t276 * t241;
t300 = qJDD(2) * t242 - t354;
t327 = qJD(4) * t393;
t344 = qJD(4) * t245;
t47 = -t180 * t327 + t181 * t344 - t245 * t300 - t393 * t259;
t281 = -t245 * t180 - t181 * t393;
t48 = -qJD(4) * t281 + t245 * t259 - t393 * t300;
t116 = -t393 * t180 + t181 * t245;
t57 = t116 * t244 - t248 * t281;
t12 = qJD(6) * t57 - t244 * t47 - t248 * t48;
t349 = qJD(1) * t249;
t218 = -qJD(4) + t349;
t341 = -qJD(6) - t218;
t418 = t341 * t57;
t421 = t12 + t418;
t332 = t393 * t242;
t280 = -t245 * t241 + t332;
t404 = -t241 * t344 + t242 * t327;
t356 = -t280 * t349 + t404;
t322 = qJD(2) * pkin(2) - qJD(3);
t310 = -pkin(7) * t350 + t322;
t272 = pkin(3) * t180 + t310;
t262 = -qJ(5) * t281 + t272;
t396 = pkin(4) + pkin(5);
t26 = -t116 * t396 + t262;
t230 = cos(t238);
t293 = t229 * t248 - t230 * t244;
t303 = pkin(2) * t246 - qJ(3) * t249;
t165 = qJD(2) * t303 - qJD(3) * t246;
t304 = pkin(2) * t249 + qJ(3) * t246;
t194 = -pkin(1) - t304;
t113 = qJD(1) * t165 + qJDD(1) * t194;
t153 = -pkin(7) * t275 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t65 = t242 * t113 - t241 * t153;
t38 = pkin(3) * t275 - pkin(8) * t259 + t65;
t66 = t241 * t113 + t242 * t153;
t49 = pkin(8) * t300 + t66;
t172 = t194 * qJD(1);
t227 = pkin(7) * t349;
t199 = qJD(2) * qJ(3) + t227;
t126 = t242 * t172 - t199 * t241;
t336 = pkin(3) * t349;
t74 = -pkin(8) * t181 + t126 - t336;
t127 = t241 * t172 + t242 * t199;
t83 = pkin(8) * t180 + t127;
t323 = t245 * t49 + t83 * t327 + t74 * t344 - t393 * t38;
t305 = qJDD(5) + t323;
t3 = pkin(9) * t47 - t185 * t396 + t305;
t173 = t185 * qJ(5);
t202 = t218 * qJD(5);
t282 = t245 * t38 + t74 * t327 - t344 * t83 + t393 * t49;
t7 = t173 - t202 + t282;
t5 = pkin(9) * t48 + t7;
t333 = -t244 * t5 + t248 * t3;
t384 = g(3) * t246;
t362 = t247 * t249;
t154 = t229 * t362 + t230 * t250;
t360 = t250 * t229;
t155 = t230 * t362 - t360;
t403 = t154 * t248 - t155 * t244;
t363 = t247 * t230;
t156 = t249 * t360 - t363;
t361 = t249 * t250;
t157 = t229 * t247 + t230 * t361;
t85 = t156 * t248 - t157 * t244;
t420 = -g(1) * t85 - g(2) * t403 - t26 * t57 - t293 * t384 + t333;
t419 = t116 ^ 2;
t397 = t281 ^ 2;
t189 = t303 * qJD(1);
t139 = pkin(7) * t328 + t242 * t189;
t365 = t242 * t249;
t289 = pkin(3) * t246 - pkin(8) * t365;
t103 = qJD(1) * t289 + t139;
t170 = t241 * t189;
t366 = t242 * t246;
t367 = t241 * t249;
t278 = -pkin(7) * t366 - pkin(8) * t367;
t123 = qJD(1) * t278 + t170;
t417 = qJD(3) * t332 - t393 * t123 - t197 * t327 + (-qJD(3) * t241 - qJD(4) * t198 - t103) * t245;
t416 = t116 * t218;
t415 = t218 * t281;
t297 = -t248 * t116 - t244 * t281;
t414 = t297 * t341;
t187 = t241 * t393 + t245 * t242;
t169 = t187 * qJD(4);
t269 = t249 * t187;
t355 = -qJD(1) * t269 + t169;
t412 = -t297 ^ 2 + t57 ^ 2;
t378 = -qJ(5) * t350 + t417;
t409 = t57 * t297;
t408 = t187 * qJD(3) + qJD(4) * t134 + t103 * t393 - t245 * t123;
t33 = -t245 * t83 + t393 * t74;
t358 = qJD(5) - t33;
t175 = t241 * t336 + t227;
t405 = -qJ(5) * t356 - qJD(5) * t187 - t175;
t402 = -pkin(7) * t325 - qJDD(3);
t342 = qJD(6) * t248;
t343 = qJD(6) * t244;
t319 = -t116 * t342 - t244 * t48 + t248 * t47 - t281 * t343;
t401 = t319 + t414;
t292 = t229 * t244 + t230 * t248;
t296 = t154 * t244 + t155 * t248;
t86 = t156 * t244 + t157 * t248;
t399 = -g(1) * t86 - g(2) * t296 - t26 * t297 - t292 * t384;
t398 = -0.2e1 * pkin(1);
t392 = pkin(3) * t241;
t391 = pkin(4) * t185;
t390 = pkin(7) * t180;
t389 = g(1) * t247;
t386 = g(2) * t250;
t294 = -t187 * t244 - t248 * t280;
t382 = qJD(6) * t294 + t244 * t355 + t248 * t356;
t122 = t187 * t248 - t244 * t280;
t381 = qJD(6) * t122 + t244 * t356 - t248 * t355;
t380 = -t355 * t396 - t405;
t379 = pkin(4) * t355 + t405;
t34 = t245 * t74 + t393 * t83;
t377 = pkin(4) * t350 + t408;
t376 = t218 * t34;
t203 = t218 * qJ(5);
t23 = pkin(9) * t116 + t34;
t21 = -t203 + t23;
t375 = t244 * t21;
t373 = qJ(5) * t116;
t372 = qJDD(2) * pkin(2);
t371 = t281 * t116;
t368 = t241 * t246;
t364 = t246 * t250;
t359 = -pkin(9) * t281 - t358;
t179 = t242 * t194;
t125 = -pkin(8) * t366 + t179 + (-pkin(7) * t241 - pkin(3)) * t249;
t147 = pkin(7) * t365 + t241 * t194;
t132 = -pkin(8) * t368 + t147;
t357 = t245 * t125 + t393 * t132;
t346 = qJD(2) * t246;
t335 = pkin(7) * t346;
t130 = t242 * t165 + t241 * t335;
t225 = pkin(7) * t338;
t353 = -t225 - t237;
t345 = qJD(2) * t249;
t330 = t241 * t345;
t176 = pkin(3) * t330 + pkin(7) * t345;
t190 = pkin(3) * t368 + t246 * pkin(7);
t352 = t250 * pkin(1) + t247 * pkin(7);
t239 = t246 ^ 2;
t351 = -t249 ^ 2 + t239;
t20 = t218 * t396 - t359;
t337 = t20 * t342 + t244 * t3 + t248 * t5;
t334 = t246 * t396;
t222 = t242 * pkin(3) + pkin(2);
t320 = t341 ^ 2;
t220 = t246 * t389;
t317 = -g(2) * t364 + t220;
t96 = -pkin(9) * t280 + t134;
t316 = pkin(9) * t356 - qJD(1) * t334 + qJD(6) * t96 - t408;
t133 = t197 * t393 + t245 * t198;
t95 = -t187 * pkin(9) + t133;
t315 = -pkin(9) * t355 - qJD(6) * t95 - t378;
t314 = -g(1) * t154 + g(2) * t156;
t313 = g(1) * t155 - g(2) * t157;
t311 = -t386 + t389;
t60 = -qJ(5) * t249 + t357;
t162 = t280 * t246;
t308 = qJ(5) * t162 - t190;
t306 = t125 * t393 - t245 * t132;
t10 = t244 * t20 + t248 * t21;
t61 = t249 * pkin(4) - t306;
t39 = t249 * pkin(5) - t162 * pkin(9) + t61;
t161 = t187 * t246;
t41 = pkin(9) * t161 + t60;
t302 = -t244 * t41 + t248 * t39;
t301 = t244 * t39 + t248 * t41;
t299 = qJ(5) * t248 - t244 * t396;
t298 = qJ(5) * t244 + t248 * t396;
t295 = t248 * t161 - t162 * t244;
t91 = t161 * t244 + t162 * t248;
t288 = qJ(5) * t187 + t222;
t287 = pkin(4) * t230 + qJ(5) * t229 + t222;
t286 = -t225 + t372 + t402;
t285 = t21 * t343 - t337;
t283 = -pkin(7) * qJDD(2) + t340 * t398;
t151 = t241 * t165;
t106 = qJD(2) * t278 + t151;
t92 = qJD(2) * t289 + t130;
t279 = t245 * t106 + t125 * t344 + t132 * t327 - t393 * t92;
t277 = t393 * t106 + t125 * t327 - t132 * t344 + t245 * t92;
t274 = t242 * t338 + t339;
t253 = qJD(1) ^ 2;
t273 = pkin(1) * t253 + t312;
t252 = qJD(2) ^ 2;
t271 = pkin(7) * t252 + qJDD(1) * t398 + t386;
t97 = t169 * t246 + t245 * t330 - t332 * t345;
t270 = -qJ(5) * t97 + qJD(5) * t162 - t176;
t268 = g(2) * t246 * t363 - t133 * t185 + (g(1) * t364 - t237) * t230;
t266 = -t249 * t312 - t384;
t265 = g(1) * t156 + g(2) * t154 + t229 * t384 - t323;
t264 = -t284 - t372;
t263 = t47 - t416;
t261 = t286 + t413;
t16 = qJ(5) * t346 - qJD(5) * t249 + t277;
t93 = -pkin(3) * t300 - t286;
t40 = pkin(4) * t116 - t262;
t258 = -t281 * t40 + qJDD(5) - t265;
t256 = g(1) * t157 + g(2) * t155 - t218 * t33 + t230 * t384 - t282;
t255 = -t47 * qJ(5) - qJD(5) * t281 - t93;
t13 = t48 * pkin(4) - t255;
t254 = t48 + t415;
t235 = t250 * pkin(7);
t174 = -qJDD(6) + t185;
t146 = -pkin(7) * t367 + t179;
t140 = -pkin(7) * t331 + t170;
t131 = -t242 * t335 + t151;
t114 = -pkin(4) * t280 - t288;
t98 = qJD(2) * t269 + t246 * t404;
t80 = t280 * t396 + t288;
t75 = pkin(4) * t161 - t308;
t63 = -pkin(4) * t281 + t373;
t62 = -t161 * t396 + t308;
t35 = t281 * t396 - t373;
t30 = -t203 + t34;
t29 = pkin(4) * t218 + t358;
t28 = pkin(4) * t98 - t270;
t27 = -t47 - t416;
t25 = qJD(6) * t91 - t244 * t97 - t248 * t98;
t24 = qJD(6) * t295 + t244 * t98 - t248 * t97;
t19 = -t396 * t98 + t270;
t18 = -pkin(4) * t346 + t279;
t15 = pkin(9) * t98 + t16;
t14 = t97 * pkin(9) - qJD(2) * t334 + t279;
t9 = t20 * t248 - t375;
t8 = t305 - t391;
t6 = -t396 * t48 + t255;
t1 = [qJDD(1), t311, t312, qJDD(1) * t239 + 0.2e1 * t246 * t325, 0.2e1 * t231 * t246 - 0.2e1 * t340 * t351, qJDD(2) * t246 + t249 * t252, qJDD(2) * t249 - t246 * t252, 0, t283 * t246 + (-t271 + t389) * t249, t246 * t271 + t249 * t283 - t220, -t312 * t241 + (-pkin(7) * t300 - t286 * t241 + (qJD(1) * t146 + t126) * qJD(2)) * t246 + (-qJD(1) * t130 - qJDD(1) * t146 - t65 + t311 * t242 + (-t241 * t310 - t390) * qJD(2)) * t249, -t312 * t242 + (-t286 * t242 + (-qJD(1) * t147 - t127) * qJD(2) + t274 * pkin(7)) * t246 + (t131 * qJD(1) + t147 * qJDD(1) + t66 - t311 * t241 + (-t310 * t242 + (t181 + t331) * pkin(7)) * qJD(2)) * t249, t131 * t180 - t147 * t354 - t130 * t181 + (-qJDD(2) * t146 - t127 * t345 - t246 * t66) * t241 + (t147 * qJDD(2) - t126 * t345 - t146 * t276 - t65 * t246) * t242 + t317, t66 * t147 + t127 * t131 + t65 * t146 + t126 * t130 - g(1) * t235 - g(2) * (t250 * t304 + t352) - t194 * t389 + (-t246 * t286 - t310 * t345) * pkin(7), -t162 * t47 + t281 * t97, t116 * t97 + t161 * t47 - t162 * t48 + t281 * t98, t162 * t185 + t218 * t97 + t249 * t47 - t281 * t346, -t116 * t346 - t161 * t185 + t218 * t98 + t249 * t48, -t185 * t249 - t218 * t346, t176 * t116 + t93 * t161 + t185 * t306 + t190 * t48 + t218 * t279 + t249 * t323 - t272 * t98 + t33 * t346 + t313, t93 * t162 - t176 * t281 - t185 * t357 - t190 * t47 + t218 * t277 + t249 * t282 + t272 * t97 - t34 * t346 + t314, t116 * t28 + t13 * t161 + t18 * t218 - t185 * t61 + t249 * t8 - t29 * t346 + t40 * t98 + t48 * t75 + t313, -t116 * t16 - t161 * t7 + t162 * t8 - t18 * t281 - t29 * t97 - t30 * t98 - t47 * t61 - t48 * t60 + t317, -t13 * t162 - t16 * t218 + t185 * t60 - t249 * t7 + t28 * t281 + t30 * t346 + t40 * t97 + t47 * t75 - t314, t7 * t60 + t30 * t16 + t13 * t75 + t40 * t28 + t8 * t61 + t29 * t18 - g(1) * (-pkin(4) * t155 - qJ(5) * t154 + t250 * t392 + t235) - g(2) * (pkin(4) * t157 + qJ(5) * t156 + t222 * t361 + t364 * t383 + t352) + (-g(1) * (-t222 * t249 - t246 * t383 - pkin(1)) - g(2) * t392) * t247, t24 * t57 - t319 * t91, -t12 * t91 - t24 * t297 - t25 * t57 - t295 * t319, -t174 * t91 - t24 * t341 - t249 * t319 - t346 * t57, -t12 * t249 - t174 * t295 + t25 * t341 + t297 * t346, -t174 * t249 + t341 * t346 -(t14 * t248 - t15 * t244) * t341 - t302 * t174 + t333 * t249 - t9 * t346 + t19 * t297 + t62 * t12 - t6 * t295 + t26 * t25 + g(1) * t296 - g(2) * t86 + (-t10 * t249 + t301 * t341) * qJD(6) (qJD(6) * t302 + t14 * t244 + t15 * t248) * t341 + t301 * t174 + t285 * t249 + t10 * t346 + t19 * t57 - t62 * t319 + t6 * t91 + t26 * t24 + g(1) * t403 - g(2) * t85; 0, 0, 0, -t246 * t253 * t249, t351 * t253, t338, t231, qJDD(2), t246 * t273 + t353, t384 + (-pkin(7) * qJDD(1) + t273) * t249, t241 * qJ(3) * t231 - pkin(2) * t354 + (t261 + t372) * t242 + ((-qJ(3) * t348 - t126) * t246 + (t390 + t139 + (qJD(3) + t310) * t241) * t249) * qJD(1), -t303 * t242 * qJDD(1) + (-t286 + t264 + t237) * t241 + ((-qJ(3) * t347 + t127) * t246 + (-pkin(7) * t181 - t140 + (t310 - t322) * t242) * t249) * qJD(1), t139 * t181 - t140 * t180 + (qJ(3) * t300 + qJD(3) * t180 + t126 * t349 + t66) * t242 + (qJ(3) * t259 + qJD(3) * t181 + t127 * t349 - t65) * t241 + t266, t310 * t227 - t126 * t139 - t127 * t140 + (-t126 * t241 + t127 * t242) * qJD(3) + t261 * pkin(2) + (-t65 * t241 + t66 * t242 + t266) * qJ(3), -t187 * t47 - t281 * t356, -t116 * t356 - t187 * t48 - t280 * t47 + t281 * t355, t185 * t187 - t218 * t356 + t281 * t350, t116 * t350 + t185 * t280 + t218 * t355, t218 * t350, -t175 * t116 + t218 * t408 - t222 * t48 - t272 * t355 - t280 * t93 - t33 * t350 + t268, t175 * t281 + t93 * t187 + t218 * t417 + t222 * t47 - t356 * t272 + t34 * t350 - t422, t114 * t48 + t116 * t379 - t13 * t280 + t218 * t377 + t29 * t350 + t355 * t40 + t268, -t116 * t378 - t133 * t47 - t134 * t48 + t187 * t8 + t280 * t7 - t281 * t377 + t29 * t356 - t30 * t355 + t266, t114 * t47 - t13 * t187 - t218 * t378 + t281 * t379 - t30 * t350 - t356 * t40 + t422, t13 * t114 + t8 * t133 + t7 * t134 + t379 * t40 + t378 * t30 + t377 * t29 + (-g(3) * t287 - t312 * t383) * t249 + (-g(3) * t383 + t287 * t312) * t246, -t122 * t319 + t382 * t57, -t12 * t122 - t294 * t319 - t297 * t382 - t381 * t57, -t122 * t174 - t341 * t382 + t350 * t57, -t174 * t294 - t297 * t350 + t341 * t381, -t341 * t350 -(-t244 * t96 + t248 * t95) * t174 + t80 * t12 - t6 * t294 + t380 * t297 + t381 * t26 - (t244 * t315 - t248 * t316) * t341 + t9 * t350 + t413 * t292 (t244 * t95 + t248 * t96) * t174 - t80 * t319 + t6 * t122 + t380 * t57 + t382 * t26 - (t244 * t316 + t248 * t315) * t341 - t10 * t350 + t413 * t293; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181 * t349 - t300 (-t180 + t347) * t349 + t274, -t180 ^ 2 - t181 ^ 2, t126 * t181 - t127 * t180 + t264 - t353 - t402, 0, 0, 0, 0, 0, t254, -t263, t254, -t397 - t419, t263, t116 * t30 + t281 * t29 + t13 - t413, 0, 0, 0, 0, 0, -t12 + t418, t319 - t414; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t371, t397 - t419, t27, -t48 + t415, t185, -t272 * t281 + t265 - t376, -t116 * t272 + t256, -t116 * t63 - t258 - t376 + 0.2e1 * t391, pkin(4) * t47 - qJ(5) * t48 - (t30 - t34) * t281 + (t29 - t358) * t116, -t116 * t40 - t281 * t63 + 0.2e1 * t173 - 0.2e1 * t202 - t256, t7 * qJ(5) - t8 * pkin(4) - t40 * t63 - t29 * t34 - g(1) * (-pkin(4) * t156 + qJ(5) * t157) - g(2) * (-pkin(4) * t154 + qJ(5) * t155) + t358 * t30 - (-pkin(4) * t229 + qJ(5) * t230) * t384, -t409, -t412, t401, t421, t174, t298 * t174 - t35 * t297 - (-t23 * t248 + t244 * t359) * t341 + (t299 * t341 + t10) * qJD(6) - t420, t299 * t174 - t35 * t57 - (t23 * t244 + t248 * t359) * t341 + (-t298 * t341 - t375) * qJD(6) + t337 + t399; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 - t371, t27, -t218 ^ 2 - t397, t218 * t30 + t258 - t391, 0, 0, 0, 0, 0, -t248 * t174 - t244 * t320 + t281 * t297, t244 * t174 - t248 * t320 + t281 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t409, t412, -t401, -t421, -t174 (-t341 - qJD(6)) * t10 + t420, -t341 * t9 + t285 - t399;];
tau_reg  = t1;
